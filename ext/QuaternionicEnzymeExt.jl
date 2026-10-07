# Support for Enzyme.  Enzyme differentiates the package's source natively, in forward and
# reverse mode and at higher order: quaternions are immutable `isbits` structs, so they are
# `Active` arguments in reverse mode, and the shadow of a quaternion has the primal's own
# type.  (The shadow of a `Rotor` is a `Rotor` holding the raw gradient components, so its
# `abs` reports 1; read such shadows with `components` or `vec`.)  So this extension
# imports no ChainRules rules, and it provides only
#
#   * `Enzyme.onehot` for quaternion arguments, so that `Enzyme.gradient(Forward, f, q)`
#     works when `q` is a quaternion;
#   * native rules for `Quaternionic.dominant_eigenvector_lapack`, because Enzyme cannot
#     differentiate the LAPACK call inside it, which `from_rotation_matrix` and `align`
#     reach (these rules are first-order only, so Enzyme cannot compute second
#     derivatives through those two functions);
#   * rules that work around a bug in Enzyme's own reverse rule for `hypot` with three or
#     more arguments, on `Quaternionic.hypotenuse`, through which `abs` and `absvec` of
#     real quaternions call `hypot`.
module QuaternionicEnzymeExt

import Enzyme
import Enzyme: EnzymeRules
using Enzyme: Annotation, Active, Const, Duplicated, BatchDuplicated
using Quaternionic: Quaternionic, AbstractQuaternion, QuatVec
import StaticArrays


############################################################################################
# `Enzyme.onehot`
#
# `Enzyme.gradient(Forward, f, x)` pushes forward one tangent per degree of freedom of `x`,
# and it gets those tangents from `Enzyme.onehot(x)` (or from `Enzyme.chunkedonehot` when a
# `chunk` size is given).  A quaternion is a `Number` with no such methods of its own, so
# they are defined here.  Each tangent is a quaternion of the same type as the argument
# (Enzyme requires a shadow of the primal's own type), with one component equal to 1 and the
# others 0.  A `QuatVec` has only three degrees of freedom, because its scalar part is
# always zero, so it gets three tangents, one for each vector component.  The gradient with
# respect to a quaternion argument is then a tuple of four derivatives (three for a
# `QuatVec`), in the order of `components` (or of `vec`).

# The number of degrees of freedom of `q`, which is the number of tangents `onehot` returns.
ntangents(::AbstractQuaternion) = 4
ntangents(::QuatVec) = 3

# The `i`th tangent that `onehot` returns for `q`: the quaternion of the same type as `q`
# whose `i`th component (`i`th vector component for a `QuatVec`) is 1, and whose others are
# 0.  The type's own constructor stores the components as given, so the tangent of a `Rotor`
# is not normalized.
unittangent(q::AbstractQuaternion{T}, i::Int) where {T} =
    typeof(q)(ntuple(j -> T(j == i), Val(4))...)
unittangent(q::QuatVec{T}, i::Int) where {T} =
    typeof(q)(ntuple(j -> T(j == i + 1), Val(4))...)

Enzyme.onehot(q::AbstractQuaternion) = ntuple(i -> unittangent(q, i), Val(4))
Enzyme.onehot(q::QuatVec) = ntuple(i -> unittangent(q, i), Val(3))
Enzyme.onehot(q::AbstractQuaternion, start::Int, endl::Int) =
    ntuple(i -> unittangent(q, start + i - 1), endl - start + 1)

# Enzyme's generic method counts the degrees of freedom with `length`, which is 1 for any
# `Number`, so the chunks are formed here from `ntangents` instead.
function Enzyme.chunkedonehot(q::AbstractQuaternion, ::Val{chunk}) where {chunk}
    n = ntangents(q)
    ntuple(k -> Enzyme.onehot(q, (k - 1) * chunk + 1, min(k * chunk, n)), cld(n, chunk))
end

# For an array of quaternions, Enzyme's own `onehot` puts `one` of the element type at each
# index, which perturbs only the scalar parts.  So `Enzyme.gradient(Forward, f, x)` and
# `Enzyme.jacobian(Forward, f, x)` (and DifferentiationInterface's forward-mode `gradient`)
# would silently omit the derivatives with respect to the vector parts.  These methods raise
# an error instead.  Reverse mode, and arrays of quaternions inside the function, are
# unaffected.  The `Array` and `SArray` methods resolve ambiguities with Enzyme's methods.
for A ∈ (:AbstractArray, :Array, :(StaticArrays.SArray{S,<:AbstractQuaternion} where {S}))
    E = A isa Symbol ? :($A{<:AbstractQuaternion}) : A
    @eval Enzyme.onehot(::$E) = throw(onehot_array_error())
    @eval Enzyme.onehot(::$E, ::Int, ::Int) = throw(onehot_array_error())
end
function onehot_array_error()
    ArgumentError(
        "Forward-mode gradients and Jacobians with respect to an array of quaternions would " *
        "perturb only their scalar parts.  Pass the components as a real array instead, for " *
        "example with `to_float_array` and `from_float_array`, or use reverse mode."
    )
end


############################################################################################
# Rules for `dominant_eigenvector_lapack`
#
# `Quaternionic.dominant_eigenvector_lapack(A, uplo)` is the only place where
# `from_rotation_matrix` and `align` call LAPACK, and Enzyme cannot differentiate LAPACK.
# The derivative formulas are plain functions in `src/conversion.jl`, shared with the
# ChainRulesCore extension; they read only the triangle of `A` (or of its tangent) named by
# `uplo`, and the pullback returns a dense matrix that is zero in the other triangle.  The
# rules handle any batch width.
#
# A matrix that is inactive gets a zero tangent and no cotangent.  It may arrive as a
# `Const`, but with runtime activity enabled, Enzyme passes a matrix that turns out to be
# inactive at run time as a `Duplicated` (or `BatchDuplicated`) whose shadow is the primal
# matrix itself.  Accumulating a cotangent into that shadow would overwrite the caller's
# matrix, so such a shadow is recognized by identity, as in Enzyme's own rules for
# `LinearAlgebra`, and skipped.
#
# These rules are first-order only.  Their bodies call LAPACK (through `eigen` and the
# dense solve in `bordered_solve`), which Enzyme cannot differentiate, so Enzyme cannot
# compute second derivatives through `from_rotation_matrix` or `align`; it throws an error
# there.  ForwardDiff can, through the value-refined eigenvector in `src`.

const LAPACKFloat = Union{Float32,Float64}

# The `k`th of the `N` shadows of the argument `A`.
shadowof(A::Annotation, k, N) = N == 1 ? A.dval : A.dval[k]

# Whether the `k`th of the `N` shadows of the matrix argument `A` stands for an inactive
# matrix: `A` is a `Const`, or runtime activity is enabled and the shadow is `A` itself.
isinactive(config, ::Const, k, N) = true
isinactive(config, A::Annotation, k, N) =
    EnzymeRules.runtime_activity(config) && shadowof(A, k, N) === A.val

# The tangent of `v` along the `k`th of the `N` shadows of the matrix argument `A`.
function eigenvectortangent(config, A::Annotation, uplo, v, k, N)
    isinactive(config, A, k, N) && return zero(v)
    Quaternionic.dominant_eigenvector_pushforward(A.val, uplo, v, shadowof(A, k, N))
end

function EnzymeRules.forward(
    config::EnzymeRules.FwdConfig,
    func::Const{typeof(Quaternionic.dominant_eigenvector_lapack)},
    ::Type{RT},
    A::Annotation{Matrix{T}},
    uplo::Const{Char}
) where {RT,T<:LAPACKFloat}
    v = func.val(A.val, uplo.val)
    N = EnzymeRules.width(config)
    if EnzymeRules.needs_shadow(config)
        dv = if N == 1
            eigenvectortangent(config, A, uplo.val, v, 1, N)
        else
            ntuple(k -> eigenvectortangent(config, A, uplo.val, v, k, N), Val(N))
        end
        if EnzymeRules.needs_primal(config)
            return N == 1 ? Duplicated(v, dv) : BatchDuplicated(v, dv)
        else
            return dv
        end
    elseif EnzymeRules.needs_primal(config)
        return v
    else
        return nothing
    end
end

# The forward pass stores copies of `A` and of `v` on the tape, in case the caller
# overwrites either one before the reverse pass, along with the shadow of `v`, into which
# the caller accumulates the cotangent of `v`.
function EnzymeRules.augmented_primal(
    config::EnzymeRules.RevConfig,
    func::Const{typeof(Quaternionic.dominant_eigenvector_lapack)},
    ::Type{RT},
    A::Annotation{Matrix{T}},
    uplo::Const{Char}
) where {RT,T<:LAPACKFloat}
    v = func.val(A.val, uplo.val)
    N = EnzymeRules.width(config)
    dv = N == 1 ? zero(v) : ntuple(_ -> zero(v), Val(N))
    primal = EnzymeRules.needs_primal(config) ? v : nothing
    shadow = EnzymeRules.needs_shadow(config) ? dv : nothing
    return EnzymeRules.AugmentedReturn(primal, shadow, (copy(A.val), copy(v), dv))
end

# Accumulate the cotangent of `A` into its `k`th shadow, given the `k`th cotangent `v̄` of
# `v`, unless that shadow stands for an inactive matrix, and zero `v̄`, which the forward
# pass allocated.
function accumulate_eigenvector_cotangent!(config, A, k, N, A0, uplo, v, v̄)
    if !isinactive(config, A, k, N)
        shadowof(A, k, N) .+= Quaternionic.dominant_eigenvector_pullback(A0, uplo, v, v̄)
    end
    fill!(v̄, 0)
    nothing
end

function EnzymeRules.reverse(
    config::EnzymeRules.RevConfig,
    func::Const{typeof(Quaternionic.dominant_eigenvector_lapack)},
    ::Type{RT},
    tape,
    A::Annotation{Matrix{T}},
    uplo::Const{Char}
) where {RT,T<:LAPACKFloat}
    A0, v, dv = tape
    N = EnzymeRules.width(config)
    for k in 1:N
        accumulate_eigenvector_cotangent!(config, A, k, N, A0, uplo.val, v, N == 1 ? dv : dv[k])
    end
    return (nothing, nothing)
end


############################################################################################
# Rules for `hypotenuse`
#
# WORKAROUND for two Enzyme bugs, one of which crashes the process
# (EnzymeAD/Enzyme.jl#3775).  These rules, and the function `hypotenuse` in
# `src/math.jl`, exist only to avoid those bugs, and should be removed once both are fixed.
#
# For real components, `abs(q)` is `hypotenuse(components(q)...)` and `absvec(q)` is
# `hypotenuse(vec(q)...)`, where `Quaternionic.hypotenuse` is `hypot` of three or four
# numbers, and `rotor` normalizes with `abs`.  Enzyme's own reverse rule for `hypot` with
# three or more arguments (`_hypotreverse` in Enzyme's `src/internal_rules/math.jl`, as of
# Enzyme 0.13.211) reads the cotangent as `dret.val[i]` when the batch width is greater than
# 1, but Enzyme passes a tuple of `Active` values there, so it throws `FieldError: type
# Tuple has no field val`.  Batched reverse mode, which DifferentiationInterface uses for
# Jacobians, therefore failed for nearly every function that builds a `Rotor`.  These rules
# compute the same derivatives as Enzyme's rules for `hypot`: the derivative with respect to
# each argument `x` is `x / h`, with `h` the result, or zero where `h` is zero, and they
# handle any batch width.  Forward mode has no such bug, but a function with reverse rules
# needs a forward rule as well, because Enzyme calls it when it differentiates the reverse
# pass in forward mode (forward-over-reverse Hessians).  These rules can be removed once
# Enzyme's reverse rule is fixed.
#
# The rules are defined on `hypotenuse`, whose arguments are numbers, rather than on `abs`
# and `absvec`, whose argument is a quaternion.  Enzyme (as of 0.13.211) adds the cotangent
# that a reverse rule returns for an `Active` struct argument into the shadow of that
# argument with a vector load and store that claim the alignment of a vector of all the
# components (32 bytes for four `Float64`s), but the shadow is aligned only to 8 bytes.  On
# x86_64 processors with AVX, that store crashes the process (with a segmentation fault or
# an access violation) whenever the shadow happens not to be 32-byte aligned.  A number is
# passed by value, so its cotangent needs no such store.
#
# The workaround helps only with recent releases of Enzyme (it was verified with 0.13.205
# and 0.13.209).  Older releases of Enzyme 0.13 (at least up to 0.13.85) refuse every custom
# reverse rule with an active result at a batch width greater than 1 ("Not yet supported:
# Enzyme custom rule of batch size ..."), which includes these rules as well as Enzyme's own
# rule for `hypot`, so batched reverse mode through `abs` fails there either way.

const HypotFloat = Base.IEEEFloat

# The divisor in the derivatives of `h = hypot(xs...)`, which is `h` itself, or 1 where `h`
# is zero, as in Enzyme's rules for `hypot`.
hypotdivisor(h) = iszero(h) ? one(h) : h

# The `k`th of the `N` tangents of the argument `x`, which is zero for a constant.
argumenttangent(x::Const, k, N) = zero(x.val)
argumenttangent(x::Annotation, k, N) = N == 1 ? x.dval : x.dval[k]

# The tangent of `h = hypotenuse(xs...)` along the `k`th of the `N` tangents of the
# arguments `xs`.
hypottangent(xs, h, k, N) =
    +(map(x -> x.val * argumenttangent(x, k, N), xs)...) / hypotdivisor(h)

# The cotangent of the result in reverse mode: a single `Active` value, or a tuple of them
# for a batch width greater than 1, or the type of an inactive result, whose cotangent is
# zero.
resultcotangent(dret::Active, h, k) = dret.val
resultcotangent(dret::Tuple, h, k) = dret[k].val
resultcotangent(::Type, h, k) = zero(h)

# The cotangent of the argument `x` of `h = hypotenuse(xs...)`, given the cotangent `d` of
# `h` (a tuple of them for a batch width greater than 1): nothing for a constant argument.
argumentcotangent(::Const, h, d) = nothing
argumentcotangent(x::Active, h, d::Number) = x.val * d / hypotdivisor(h)
argumentcotangent(x::Active, h, d::Tuple) = map(dₖ -> argumentcotangent(x, h, dₖ), d)

# The first argument `x` of each rule is separate from the others, `xs`, only so that the
# type parameter `T` is bound.
function EnzymeRules.forward(
    config::EnzymeRules.FwdConfig, func::Const{typeof(Quaternionic.hypotenuse)},
    ::Type{RT}, x::Annotation{T}, xs::Annotation{T}...
) where {RT,T<:HypotFloat}
    args = (x, xs...)
    h = func.val(map(a -> a.val, args)...)
    N = EnzymeRules.width(config)
    if EnzymeRules.needs_shadow(config)
        dh = if N == 1
            hypottangent(args, h, 1, N)
        else
            ntuple(k -> hypottangent(args, h, k, N), Val(N))
        end
        if EnzymeRules.needs_primal(config)
            return N == 1 ? Duplicated(h, dh) : BatchDuplicated(h, dh)
        else
            return dh
        end
    elseif EnzymeRules.needs_primal(config)
        return h
    else
        return nothing
    end
end

function EnzymeRules.augmented_primal(
    config::EnzymeRules.RevConfig, func::Const{typeof(Quaternionic.hypotenuse)},
    ::Type{RT}, x::Union{Const{T},Active{T}}, xs::Union{Const{T},Active{T}}...
) where {RT,T<:HypotFloat}
    h = func.val(x.val, map(a -> a.val, xs)...)
    primal = EnzymeRules.needs_primal(config) ? h : nothing
    return EnzymeRules.AugmentedReturn(primal, nothing, h)
end

function EnzymeRules.reverse(
    config::EnzymeRules.RevConfig, func::Const{typeof(Quaternionic.hypotenuse)}, dret, h,
    x::Union{Const{T},Active{T}}, xs::Union{Const{T},Active{T}}...
) where {T<:HypotFloat}
    N = EnzymeRules.width(config)
    d = if N == 1
        resultcotangent(dret, h, 1)
    else
        ntuple(k -> resultcotangent(dret, h, k), Val(N))
    end
    return map(a -> argumentcotangent(a, h, d), (x, xs...))
end

end # module
