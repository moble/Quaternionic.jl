module QuaternionicReverseDiffExt

# ReverseDiff differentiates quaternions whose components are `TrackedReal`s natively, by
# recording the scalar operations on their components, so this extension provides no
# derivative rules.  It only makes the source take the same branches as it does for plain
# numbers, fixes broadcasting of a `TrackedReal` with an array of quaternions, and keeps
# replayed tapes of `from_rotation_matrix` and `align` correct.  It is loaded only as a
# package extension (Julia 1.9 and later), never through Requires.

import Quaternionic
import Quaternionic: AbstractQuaternion
import ReverseDiff
import ReverseDiff: TrackedReal, TrackedArray, SpecialInstruction
using Base.Broadcast: Broadcasted
import LinearAlgebra
using LinearAlgebra: Symmetric, ⋅
using StaticArrays: SMatrix, SVector

# Strip the tracking, so that `iszerovalue` and the value-based series thresholds take the
# same branches as they do for the underlying numbers, and so that `dominant_eigenvector`
# treats a tracked element type as a wrapped one and uses its value-refined eigenvector.
# The recursion also strips nested tracking, as in `ReverseDiff.hessian`, and dual numbers
# inside `TrackedReal`s.
Quaternionic.value(x::TrackedReal) = Quaternionic.value(ReverseDiff.value(x))

# ReverseDiff requires the function it differentiates to return a real number (for
# `gradient`) or an array of real numbers (for `jacobian`).  It treats any other `Number`
# that is not a `TrackedReal` as a constant, so a quaternion output, or an array of them,
# would silently give a gradient or Jacobian of zeros.  These methods raise an error
# instead.  Quaternions inside the function are unaffected; only the final output matters.
function ReverseDiff.seeded_reverse_pass!(_, ::AbstractQuaternion, _, _)
    throw(ArgumentError(
        "ReverseDiff.gradient requires a function that returns a real number, but this " *
        "function returns a quaternion.  Return a real number computed from it, or " *
        "differentiate `x -> collect(components(f(x)))` with `ReverseDiff.jacobian`."
    ))
end
function ReverseDiff.seeded_reverse_pass!(
    ::AbstractArray, ::AbstractArray{<:AbstractQuaternion}, ::TrackedArray, _
)
    throw(ArgumentError(
        "ReverseDiff.jacobian requires a function that returns an array of real numbers, " *
        "but this function returns an array of quaternions.  Return their components " *
        "instead, for example with `to_float_array`."
    ))
end

# ReverseDiff specializes `materialize` for `f.(x, A)` and `f.(A, x)`, where `x` is a
# `TrackedReal`, `A` is an array of numbers, and `f` is any binary function with a DiffRules
# rule.  Those methods call `broadcast(f, x, A)`, assuming that `A` holds plain reals, but
# for an array of quaternions that call builds the same `Broadcasted` object and reaches the
# same method again, until the stack overflows.  For the binary functions that quaternions
# support, the methods below broadcast elementwise instead, so that each element records its
# own scalar operations.  The array types are the six that ReverseDiff specializes, so that
# each of its methods is overridden by one of these, which is strictly more specific.

const RDBroadcasted{F,T} = Broadcasted{<:Any,<:Any,F,T}

# Compute the broadcast `bc` elementwise.  This is what Base's `materialize` does, and it
# keeps the style of the broadcast, so that static arrays give static results.
elementwise(bc::Broadcasted) = copy(Base.Broadcast.instantiate(bc))

for f ∈ (:+, :-, :*, :/, :^)
    for A ∈ (:AbstractArray, :AbstractVector, :AbstractMatrix, :Array, :Vector, :Matrix)
        @eval begin
            function Base.Broadcast.materialize(
                bc::RDBroadcasted{typeof($f),<:Tuple{$A{<:AbstractQuaternion},TrackedReal}}
            )
                elementwise(bc)
            end
            function Base.Broadcast.materialize(
                bc::RDBroadcasted{typeof($f),<:Tuple{TrackedReal,$A{<:AbstractQuaternion}}}
            )
                elementwise(bc)
            end
        end
    end
end

# `from_rotation_matrix` and `align` find the dominant eigenvector of a symmetric 4×4
# matrix.  For wrapped element types, `refined_dominant_eigenvector` in src/conversion.jl
# starts from the eigenvector of the matrix of values, which LAPACK computes, and refines it
# by two Newton steps, through which the derivatives are propagated.  ReverseDiff would
# record that starting vector as a constant, so a tape replayed at another input (as by
# `ReverseDiff.gradient!` or `ReverseDiff.hessian!` with a recorded or compiled tape) would
# start from the eigenvector of the recording point, and its results would be wrong unless
# the two points were very close.  The method below is the same computation, except that
# the starting vector is recorded as an instruction of its own, whose forward pass
# recomputes it from the current values and whose reverse pass propagates nothing, because
# the starting vector is a constant at every input.  Each Newton step must stay identical
# to that of the source method.
function Quaternionic.refined_dominant_eigenvector(M::Symmetric{T}) where {T<:TrackedReal}
    A = SMatrix{4,4,T}(M)
    v = SVector{4,T}(starting_vector(A))
    λ = v ⋅ (A * v)
    for _ ∈ 1:2
        r = A * v - λ * v
        δ = Quaternionic.bordered_matrix(A - λ * LinearAlgebra.I, v) \
            SVector{5,T}(-r[1], -r[2], -r[3], -r[4], (1 - v ⋅ v) / 2)
        v += SVector{4,T}(δ[1], δ[2], δ[3], δ[4])
        λ -= δ[5]
    end
    v
end

# The starting vector of `refined_dominant_eigenvector` for the matrix `A` of `TrackedReal`s:
# the dominant unit eigenvector of the matrix of values, as a tuple of `TrackedReal`s on the
# tape of `A`.  The instruction stores the vector of plain numbers computed when it was
# recorded.  On a replay, the recomputed vector is given the sign closer to that one, so
# that the sign chosen downstream at the recording point (in `positive_hemisphere`, whose
# test a tape freezes) stays consistent.  When the values of `A` are themselves tracked, as
# under `ReverseDiff.hessian`, the vector of values is found by the same function one level
# down, so that it is also recomputed when the outer tape (such as a `HessianTape`) is
# replayed.
function starting_vector(A::SMatrix{4,4,<:TrackedReal})
    a = first(A)
    V, D = ReverseDiff.valtype(a), ReverseDiff.derivtype(a)
    v = if V <: TrackedReal
        starting_vector(map(ReverseDiff.value, A))
    else
        Tuple(Quaternionic.dominant_eigenvector(Symmetric(map(Quaternionic.value, A))))
    end
    tp = ReverseDiff.tape(A)
    output = ntuple(k -> TrackedReal{V,D,Nothing}(convert(V, v[k]), zero(D), tp), 4)
    ReverseDiff.record!(
        tp, SpecialInstruction, starting_vector, A, output, map(Quaternionic.value, v)
    )
    output
end

function ReverseDiff.special_forward_exec!(instruction::SpecialInstruction{typeof(starting_vector)})
    A, output, v₀ = instruction.input, instruction.output, instruction.cache
    ReverseDiff.pull_value!(A)
    v = Quaternionic.dominant_eigenvector(Symmetric(map(Quaternionic.value, A)))
    s = sum(v[k] * v₀[k] for k ∈ 1:4) < 0 ? -1 : 1
    for k ∈ 1:4
        ReverseDiff.value!(output[k], convert(ReverseDiff.valtype(output[k]), s * v[k]))
    end
    nothing
end

function ReverseDiff.special_reverse_exec!(instruction::SpecialInstruction{typeof(starting_vector)})
    ReverseDiff.unseed!(instruction.output)
    nothing
end

end # module
