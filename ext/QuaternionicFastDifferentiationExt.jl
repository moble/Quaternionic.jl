module QuaternionicFastDifferentiationExt

using StaticArrays: SVector
import Quaternionic: Quaternionic, absvec,
    AbstractQuaternion, Quaternion, Rotor, QuatVec,
    quaternion, rotor, quatvec,
    QuatVecF64, RotorF64, QuaternionF64,
    wrapper, components, basetype
using PrecompileTools
using FastDifferentiation: FastDifferentiation, Node


Base.abs(q::AbstractQuaternion{Node}) = √sum(x->x^2, components(q))
Base.abs(q::QuatVec{Node}) = √sum(x->x^2, vec(q))
absvec(q::AbstractQuaternion{Node}) = √sum(x->x^2, vec(q))


### Functions that used to appear in quaternion.jl
# The other components are filled with `zero(w)` rather than `false`, because `Node(false)`
# prints as `false`.
quaternion(w::Node) = quaternion(SVector{4}(w, zero(w), zero(w), zero(w)))
function rotor(w::Node)
    # When `w` is a constant, the components are those of src's `rotor(v)` for that
    # constant, so that the result is `±1` for a finite nonzero constant and a NaN rotor for a
    # zero, NaN, or infinite constant.  The sign of a variable or any other expression is
    # unknown, so the identity rotor is returned without symbolic normalization.
    v = FastDifferentiation.is_constant(w) ? FastDifferentiation.value(w) : nothing
    c = if v isa Real
        SVector{4,Node}(components(rotor(v)))
    else
        SVector{4,Node}(one(w), zero(w), zero(w), zero(w))
    end
    Rotor{Node}(c)
end
quatvec(w::Node) = quatvec(SVector{4}(zero(w), zero(w), zero(w), zero(w)))
for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Node}) where {OP} = quaternion
        wrapper(::Type{<:Node}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = quaternion
    end
end
let NT = Node
    for QT ∈ (QuatVec,)
        for OP ∈ (Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = quatvec
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = quatvec
            end
        end
    end
    for QT ∈ (Rotor,)
        for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = quaternion
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = quaternion
            end
        end
    end
end
let T = Node
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval begin
            wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = quaternion
            wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = quaternion
        end
    end
end
Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Node} =
    wrapper(Q){promote_type(basetype(Q), S)}
# As in src/quaternion.jl, a `QuatVec` or a `Rotor` promoted with a scalar becomes a
# `Quaternion`.  These methods also resolve ambiguities with the rules in src.
for QT ∈ (QuatVec, Rotor)
    @eval Base.promote_rule(::Type{$QT{T}}, ::Type{S}) where {T<:Number,S<:Node} =
        Quaternion{promote_type(T, S)}
end


### Functions that used to appear in algebra.jl
for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    let TB = Node
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end
let S = Node
    @eval begin
        Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(*), $S)(s*components(p))
        Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion} = wrapper($S, Val(*), Q)(s*components(p))
        Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(/), $S)(components(p)/s)
        function Base.:/(s::$S, p::Q) where {Q<:AbstractQuaternion}
            f = s / abs2(p)
            wrapper($S, Val(/), Q)(p[1] * f, -p[2] * f, -p[3] * f, -p[4] * f)
        end
    end
end

# Here, we disable FastDifferentiation support for functions that cannot yet be used with FD
# variables, so that they raise an explanatory error.  There are two obstacles.  First, most
# of these functions choose an algorithm with Julia `if` statements that depend on the
# values of the components, which fail for FD variables with errors like
#
#    TypeError: non-boolean (Node) used in boolean context
#
# FastDifferentiation provides `if_else` for such conditionals, but these functions would
# have to be rewritten to use it.  Second, the integer powers do not branch on the
# components, but FastDifferentiation itself returns incorrect Jacobians for some of them
# (such as `q^3`), and throws a `KeyError` for others (such as `rotor(q)^2`).  The same
# upstream problem also affects explicit products like `q*q*q`, and expressions involving
# the normalization of a `Rotor`, which are not disabled here.  Non-integer powers of a
# `Quaternion` are not listed, because they call `log`, which raises its own error.  These
# stubs can be removed once both problems are solved, and these functions are tested.
let conditionals = "it chooses an algorithm with conditionals on the values of the " *
        "components, which FastDifferentiation does not support",
    upstream = "FastDifferentiation computes incorrect derivatives for some integer powers"
    stubs = [
        (:(Quaternionic.to_euler_phases(::AbstractQuaternion{Node})), conditionals),
        (:(Base.log(::Quaternion{Node})), conditionals),
        (:(Base.log(::Rotor{Node})), conditionals),
        (:(Base.exp(::Quaternion{Node})), conditionals),
        (:(Base.exp(::QuatVec{Node})), conditionals),
        (:(Base.sqrt(::Quaternion{Node})), conditionals),
        (:(Base.sqrt(::QuatVec{Node})), conditionals),
        (:(Base.sqrt(::Rotor{Node})), conditionals),
        # The exponent types are listed separately, rather than as `Number`, to avoid
        # ambiguities with Base's `^(::Number, ::Rational)` and with methods in src that
        # take a quaternion exponent.
        (:(Base.:^(::Rotor{Node}, ::Real)), conditionals),
        (:(Base.:^(::Rotor{Node}, ::Rational)), conditionals),
        (:(Base.:^(::Rotor{Node}, ::Complex)), conditionals),
        (:(Base.:^(::Quaternion{Node}, ::Integer)), upstream),
        (:(Base.:^(::QuatVec{Node}, ::Integer)), upstream),
        (:(Base.:^(::Rotor{Node}, ::Integer)), upstream),
    ]
    # `Complex{Node}` is a valid type only when `Node <: Real`, which is true in
    # FastDifferentiation 0.4 but not in 0.3.
    if Node <: Real
        push!(stubs, (
            :(Quaternionic.from_euler_phases(::Complex{Node}, ::Complex{Node}, ::Complex{Node})),
            conditionals
        ))
    end
    for (method, reason) ∈ stubs
        # The message is built here, outside the quoted expression, so that it is complete
        # when the stub is defined.  Interpolating `func` inside the string in the quoted
        # expression would instead look up a global `func` when the stub is called.
        func = first(split(string(method), '('))
        msg = "`$(func)` cannot yet be used with FastDifferentiation variables, because " *
            "$(reason)."
        @eval $(method) = error($msg)
    end
end


# Pre-compilation

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.
    FastDifferentiation.@variables w x y z a b c d e
    s = randn(Float64)
    v = randn(QuatVecF64)
    r = randn(RotorF64)
    q = randn(QuaternionF64)
    𝓈 = w
    𝓋 = quatvec(x, y, z)
    𝓇 = rotor(a, b, c, d)
    𝓆 = quaternion(w, x, y, z)

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether they belong to
        # this package or not
        r(v)
        𝓇(𝓋)
        # Tuples, unlike vectors, keep the type of each element, so that mixed operations
        # are compiled too.
        for a ∈ (s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆)
            conj(a)
            for b ∈ (s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆)
                a * b
                a / b
                a + b
                a - b
            end
        end

    end
end


end # module
