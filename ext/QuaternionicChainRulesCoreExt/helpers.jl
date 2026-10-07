# Helpers shared by every rule file of `QuaternionicChainRulesCoreExt`.
#
# This is the fixed helper API of the AD specification (section 8).  The rule files
# (`arithmetic.jl`, `optouts.jl`, `elementary.jl`, and `geometry.jl`) may rely on exactly
# these names:
#
#   AQ, RealQ                 aliases for `AbstractQuaternion` and `AbstractQuaternion{<:Real}`
#   cot(Δ)                    normalize a cotangent (or tangent) of a quaternion-valued
#                             output to a `Quaternion`, or to `ZeroTangent()`
#   scalarcot(Δ)              the same for a scalar-valued output
#   q_(x)                     `quaternion(x)` for a quaternion, the identity for a number
#   conjcomponents(q)         the `Quaternion` with every component complex-conjugated
#   qadjoint(q)               `conj(conjcomponents(q))`, the adjoint of multiplication by `q`
#   inner(a, b)               the Hermitian inner product `Σᵢ conj(aᵢ) bᵢ` of the components
#   bilinear(a, b)            the bilinear product `Σᵢ aᵢ bᵢ` of the components (an addition
#                             to the fixed API, for pushforwards with complex components)
#   zerofill(x)               `false` for `nothing` or an `AbstractZero`, otherwise `x`
#   componentcot(i, Δ)        the quaternion whose only nonzero component is `Δ`, at index `i`
#   normalize_pullback(a, Δ)  the pullback of `a ↦ a / abs(a)` at `a`
#   normalize_pushforward(a, ȧ)  the pushforward of `a ↦ a / abs(a)` at `a`
#   projectargs(args, Δs)     project each cotangent in `Δs` onto its argument in `args`
#
# The name `q_` is kept as the specification gives it; it has no leading underscore.  Other
# internal functions of Quaternionic should be referred to by qualified name, as in
# `Quaternionic.csqrt`, rather than added to the imports of the module file.  The module file
# also makes available `LinearAlgebra`, `dot`, `norm`, `normalize`, `SVector`, `basetype`,
# `value`, and `iszerovalue`, together with everything that ChainRulesCore exports.

const AQ = AbstractQuaternion
const RealQ = AbstractQuaternion{<:Real}

"""
    cot(Δ)

Normalize a cotangent (or a tangent) that a backend may pass for a quaternion-valued output
into a `Quaternion`, or into `ZeroTangent()` when it is zero.

Thunks are unthunked; `nothing` and every `AbstractZero` become `ZeroTangent()`; a quaternion
of any type (including a stray `Rotor` or `QuatVec`) is converted to a `Quaternion` without
being renormalized; a real or complex number becomes the quaternion with only that scalar
part; a `Tangent`, a `NamedTuple` with the single field `components` or `data`, a `Tuple`, or
an `AbstractVector` of length 4 is unpacked recursively, with `nothing` and `AbstractZero`
entries treated as zeros.
"""
cot(Δ::AbstractThunk) = cot(unthunk(Δ))
cot(::Union{Nothing,AbstractZero}) = ZeroTangent()
cot(Δ::AbstractQuaternion) = quaternion(Δ)
cot(Δ::Union{Real,Complex}) = quaternion(Δ)
cot(Δ::Tangent) = cot(ChainRulesCore.backing(Δ))
cot(Δ::NamedTuple{(:components,)}) = cot(Δ.components)
cot(Δ::NamedTuple{(:data,)}) = cot(Δ.data)
function cot(Δ::Tuple)
    length(Δ) == 4 ||
        throw(DimensionMismatch("A quaternion cotangent must have 4 components, not $(length(Δ))"))
    return cotcomponents(Δ[1], Δ[2], Δ[3], Δ[4])
end
function cot(Δ::AbstractVector)
    length(Δ) == 4 ||
        throw(DimensionMismatch("A quaternion cotangent must have 4 components, not $(length(Δ))"))
    return cotcomponents(Δ[begin], Δ[begin+1], Δ[begin+2], Δ[begin+3])
end

"""
    cotcomponents(a, b, c, d)

Return the `Quaternion` with components `a`, `b`, `c`, and `d`, each of which may be a
thunk, `nothing`, or an `AbstractZero`, or `ZeroTangent()` when all four are zero.  This is
the worker of `cot` for tuples and vectors.
"""
function cotcomponents(a, b, c, d)
    a, b, c, d = unthunk(a), unthunk(b), unthunk(c), unthunk(d)
    if all(x -> x isa Union{Nothing,AbstractZero}, (a, b, c, d))
        return ZeroTangent()
    end
    return quaternion(zerofill(a), zerofill(b), zerofill(c), zerofill(d))
end

"""
    scalarcot(Δ)

Normalize a cotangent that a backend may pass for a scalar-valued output: thunks are
unthunked, `nothing` becomes `ZeroTangent()`, and a `Tangent` of a `Complex` becomes the
`Complex` number.
"""
scalarcot(Δ) = Δ
scalarcot(Δ::AbstractThunk) = scalarcot(unthunk(Δ))
scalarcot(::Union{Nothing,AbstractZero}) = ZeroTangent()
scalarcot(Δ::Tangent{<:Complex}) = Complex(zerofill(Δ.re), zerofill(Δ.im))

"""
    q_(x)

Return `quaternion(x)` for a quaternion `x` of any type, and `x` itself for any other number.
"""
q_(x::AbstractQuaternion) = quaternion(x)
q_(x::Number) = x

"""
    conjcomponents(q)

Return `q` as a `Quaternion` with each component complex-conjugated.  This is merely
`quaternion(q)` for real components.  For a number `x`, return `conj(x)`.
"""
conjcomponents(q::AbstractQuaternion{<:Real}) = quaternion(q)
conjcomponents(q::AbstractQuaternion) = quaternion(conj(q[1]), conj(q[2]), conj(q[3]), conj(q[4]))
conjcomponents(x::Number) = conj(x)

"""
    qadjoint(q)

Return the quaternion `a` such that multiplication on the left by `a` is the adjoint
(conjugate transpose, as a map on ℝ⁴ or ℂ⁴) of multiplication on the left by `q`; the same
holds for multiplication on the right.  This is `conj(conjcomponents(q))`, which is just
`conj(q)` for real components.  For a number `x`, return `conj(x)`.
"""
qadjoint(q::AbstractQuaternion) = conj(conjcomponents(q))
qadjoint(x::Number) = conj(x)

"""
    inner(a, b)

Return the Hermitian inner product `Σᵢ conj(aᵢ) bᵢ` of the components of `a` and `b`.  A
number that is not a quaternion is treated as a quaternion with only a scalar part.  For
real components, this is the same as `real(a ⋅ b)`.
"""
inner(a::AbstractQuaternion, b::AbstractQuaternion) =
    conj(a[1]) * b[1] + conj(a[2]) * b[2] + conj(a[3]) * b[3] + conj(a[4]) * b[4]
inner(a::AbstractQuaternion, b::Number) = conj(a[1]) * b
inner(a::Number, b::AbstractQuaternion) = conj(a) * b[1]
inner(a::Number, b::Number) = conj(a) * b

"""
    bilinear(a, b)

Return the bilinear product `Σᵢ aᵢ bᵢ` of the components of `a` and `b`, without complex
conjugation.  A number that is not a quaternion is treated as a quaternion with only a
scalar part.  For real components, this is the same as `inner(a, b)` and `real(a ⋅ b)`.
"""
bilinear(a::AbstractQuaternion, b::AbstractQuaternion) =
    a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + a[4] * b[4]
bilinear(a::AbstractQuaternion, b::Number) = a[1] * b
bilinear(a::Number, b::AbstractQuaternion) = a * b[1]
bilinear(a::Number, b::Number) = a * b

"""
    zerofill(x)

Return `false` (a strong zero) when `x` is `nothing` or an `AbstractZero`, the unthunked
value when `x` is a thunk, and `x` itself otherwise.
"""
zerofill(x) = x
zerofill(::Union{Nothing,AbstractZero}) = false
zerofill(x::AbstractThunk) = zerofill(unthunk(x))

"""
    componentcot(i, Δ)

Return the `Quaternion` whose `i`th component is `Δ` and whose other components are zero.
An `AbstractZero` is returned unchanged.
"""
function componentcot(i::Integer, Δ)
    z = zero(Δ)
    return quaternion(
        ifelse(i == 1, Δ, z), ifelse(i == 2, Δ, z), ifelse(i == 3, Δ, z), ifelse(i == 4, Δ, z)
    )
end
componentcot(::Integer, Δ::AbstractZero) = Δ

"""
    normalize_pullback(a, Δ)

Return the pullback of the normalization `a ↦ a / abs(a)` at the quaternion `a`, applied to
the cotangent `Δ`, as a `Quaternion`.  With `n = abs(quaternion(a))` (computed by `hypot`, so
that it neither overflows nor underflows where the primal does not) and `â = a / n`, this is
`(Δ - conjcomponents(â) * inner(â, Δ)) / conj(n)`.  For complex components, `n` is the
holomorphic square root of the spinor norm, as in the source, so this is the conjugate
transpose of the complex Jacobian.  An `AbstractZero` cotangent is returned unchanged.
"""
function normalize_pullback(a::AbstractQuaternion, Δ)
    Δq = cot(Δ)
    Δq isa AbstractZero && return Δq
    aq = quaternion(a)
    n = abs(aq)
    â = aq / n
    return (Δq - conjcomponents(â) * inner(â, Δq)) / conj(n)
end

"""
    normalize_pushforward(a, ȧ)

Return the pushforward of the normalization `a ↦ a / abs(a)` at the quaternion `a`, applied
to the tangent `ȧ`, as a `Quaternion`.  With `n` and `â` as in `normalize_pullback`,
this is `(ȧ - â * bilinear(â, ȧ)) / n`.  An `AbstractZero` tangent is returned unchanged.
"""
function normalize_pushforward(a::AbstractQuaternion, ȧ)
    ȧq = cot(ȧ)
    ȧq isa AbstractZero && return ȧq
    aq = quaternion(a)
    n = abs(aq)
    â = aq / n
    return (ȧq - â * bilinear(â, ȧq)) / n
end

"""
    projectargs(args, Δs)

Project each cotangent in the tuple `Δs` with `ProjectTo` of the corresponding argument in
the tuple `args`, leaving every `AbstractZero` unchanged.
"""
projectargs(args::Tuple, Δs::Tuple) =
    map((a, d) -> d isa AbstractZero ? d : ProjectTo(a)(d), args, Δs)
