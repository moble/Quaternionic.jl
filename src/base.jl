# Useful functions from Base

function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:AbstractQuaternion,S<:AbstractQuaternion}
    Base.rtoldefault(basetype(x), basetype(y), atol)
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:AbstractQuaternion,S<:Number}
    Base.rtoldefault(basetype(x), y, atol)
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:Number,S<:AbstractQuaternion}
    Base.rtoldefault(x, basetype(y), atol)
end

# For complex components, `abs` returns the complex spinor norm, which cannot be compared
# with a tolerance, so Base's `isapprox(::Number, ::Number)` would throw.  Instead, these
# methods compare the components as vectors, using the real Euclidean norm of the
# difference, just as `isapprox` does for vectors of complex numbers.
for (T1, T2) ∈ (
    (AbstractQuaternion{<:Complex}, AbstractQuaternion),
    (AbstractQuaternion, AbstractQuaternion{<:Complex}),
    (AbstractQuaternion{<:Complex}, AbstractQuaternion{<:Complex}),
    (AbstractQuaternion{<:Complex}, Number),
    (Number, AbstractQuaternion{<:Complex}),
)
    @eval function Base.isapprox(q1::$T1, q2::$T2; kwargs...)
        isapprox(components(quaternion(q1)), components(quaternion(q2)); kwargs...)
    end
end

Base.:(==)(q1::AbstractQuaternion{<:Number}, q2::AbstractQuaternion{<:Number}) = (q1[1]==q2[1]) && (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q::AbstractQuaternion{<:Number}, x::Number) = isreal(q) && real(q) == x
Base.:(==)(x::Number, q::AbstractQuaternion) = isreal(q) && real(q) == x
Base.isequal(q1::AbstractQuaternion, q2::AbstractQuaternion) = isequal(q1[1],q2[1]) && isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])

# A `QuatVec` has zero scalar part, so it equals another quaternion only if that quaternion
# also has zero scalar part, and it equals a number only if both are zero.  The scalar part
# that a `QuatVec` stores is always zero, but these methods treat it as zero regardless.
Base.:(==)(q1::QuatVec{<:Number}, q2::AbstractQuaternion{<:Number}) = iszero(q2[1]) && (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q1::AbstractQuaternion{<:Number}, q2::QuatVec{<:Number}) = iszero(q1[1]) && (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q1::QuatVec{<:Number}, q2::QuatVec{<:Number}) = (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q::QuatVec{<:Number}, x::Number) = iszero(q[2]) && iszero(q[3]) && iszero(q[4]) && iszero(x)
Base.:(==)(x::Number, q::QuatVec) = iszero(q[2]) && iszero(q[3]) && iszero(q[4]) && iszero(x)

# `isequal` must agree with `hash`, which hashes all four components and hashes a quaternion
# with zero vector part like its scalar part.  So, as for `Complex`, signed zeros are
# distinguished, NaN equals NaN, and a quaternion is `isequal` to a number only if its
# vector components are positive zeros and its scalar part is `isequal` to the number.  The
# scalar part of a `QuatVec` is a positive zero.  The methods comparing with a number are
# restricted to the standard number types, because packages such as Symbolics define
# `isequal` between their own types and any `Number`, which would otherwise be ambiguous.
const StandardNumber = Union{AbstractFloat,Integer,Rational,AbstractIrrational,Complex}
Base.isequal(q1::AbstractQuaternion, q2::QuatVec) = isequal(q1[1],zero(q1[1])) && isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])
Base.isequal(q1::QuatVec, q2::AbstractQuaternion) = isequal(zero(q2[1]),q2[1]) && isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])
Base.isequal(q1::QuatVec, q2::QuatVec) = isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])
Base.isequal(q::AbstractQuaternion, x::StandardNumber) = isequal(q[1],x) && isequal(q[2],zero(q[2])) && isequal(q[3],zero(q[3])) && isequal(q[4],zero(q[4]))
Base.isequal(x::StandardNumber, q::AbstractQuaternion) = isequal(q, x)
Base.isequal(q::QuatVec, x::StandardNumber) = isequal(zero(x),x) && isequal(q[2],zero(q[2])) && isequal(q[3],zero(q[3])) && isequal(q[4],zero(q[4]))

# Note that, for a quaternion with complex components, `isreal` checks only that the vector
# part is zero; the scalar part may be complex.  This matches `real(q)`, which returns the
# scalar part `q[1]`, and `q == x` for a number `x`, which compares `x` with that scalar
# part.  In the spacetime-algebra reading of complex quaternions, such an element is the sum
# of a scalar and a pseudoscalar.

Base.isreal(q::AbstractQuaternion{T}) where {T<:Number} = iszero(q[2]) && iszero(q[3]) && iszero(q[4])
Base.isinteger(q::AbstractQuaternion{T}) where {T<:Number} = isreal(q) && isinteger(real(q))
Base.isfinite(q::AbstractQuaternion{T}) where {T<:Number} = isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3]) && isfinite(q[4])
Base.isnan(q::AbstractQuaternion{T}) where {T<:Number} = isnan(q[1]) || isnan(q[2]) || isnan(q[3]) || isnan(q[4])
Base.isinf(q::AbstractQuaternion{T}) where {T<:Number} = isinf(q[1]) || isinf(q[2]) || isinf(q[3]) || isinf(q[4])
Base.iszero(q::AbstractQuaternion{T}) where {T<:Number} = iszero(q[1]) && iszero(q[2]) && iszero(q[3]) && iszero(q[4])
Base.isone(q::AbstractQuaternion{T}) where {T<:Number} = isone(q[1]) && iszero(q[2]) && iszero(q[3]) && iszero(q[4])

"""
    value(x)

This is essentially the identity function, but intended to be overridden for types like
`ForwardDiff.Dual`, where the value is the part of the dual number that corresponds to the
original function value, and the tangent part is the part that corresponds to the
derivative.  For nested types, such as the nested dual numbers used to compute higher-order
derivatives, this should strip every level of nesting.
"""
value(x) = x
value(z::Complex) = complex(value(real(z)), value(imag(z)))

"""
    iszerovalue(x)

This is essentially `Base.iszero`, but intended to be overridden for types like
`ForwardDiff.Dual`, where the value may be zero, but if the tangent is not zero, then
`Base.iszero` will return `false`.  This function, on the other hand, ignores the tangent
part, and will return `true` if and only if the value part is zero.

This is needed internally in math functions like `exp`, `log`, and `sqrt`, where we
frequently need to switch the algorithm based on whether some components are zero.  In those
isolated cases, a Taylor series should be provided that will be exactly zero for, e.g.,
`Float64`, but will also work correctly for `ForwardDiff.Dual` and other ADs.

Like `Base.iszero`, this function is defined recursively for arrays, quaternions, and
complex numbers.

"""
iszerovalue(x) = iszero(value(x))
iszerovalue(x::AbstractArray) = all(iszerovalue, x)
iszerovalue(z::Complex) = iszerovalue(real(z)) && iszerovalue(imag(z))
iszerovalue(q::AbstractQuaternion) = iszerovalue(components(q))
iszerovalue(q::QuatVec) = iszerovalue(vec(q))

Base.round(q::QT, r::RoundingMode=RoundNearest; kwargs...) where {QT<:AbstractQuaternion} = QT(round.(components(q), r; kwargs...))
# Rounding the components of a `Rotor` generally gives a quaternion whose norm is not 1, so
# the result is a `Quaternion`.
Base.round(q::Rotor{T}, r::RoundingMode=RoundNearest; kwargs...) where {T} = Quaternion(round.(components(q), r; kwargs...))

Base.in(q::AbstractQuaternion, r::AbstractRange{<:Number}) = isreal(q) && real(q) in r

Base.flipsign(q::AbstractQuaternion, x::Number) = ifelse(signbit(x), -q, q)

Base.bswap(q::Q) where {T<:Number, Q<:AbstractQuaternion{T}} = Q(bswap(q[1]), bswap(q[2]), bswap(q[3]), bswap(q[4]))

if UInt === UInt64
    const h_imagx = 0xdf13da9384000582
    const h_imagy = 0x437d0726f1028bcd
    const h_imagz = 0xcf13f7ab1f367e01
else
    const h_imagx = 0x27a4bf84
    const h_imagy = 0xccefdeeb
    const h_imagz = 0x1683854f
end
const hash_0_imagx = hash(0, h_imagx)
const hash_0_imagy = hash(0, h_imagy)
const hash_0_imagz = hash(0, h_imagz)

# A quaternion with zero vector part hashes like its scalar part, as `isequal` requires.
function Base.hash(q::AbstractQuaternion, h::UInt)
    # TODO: with default argument specialization, this would be better:
    # hash(q[1], h ⊻ hash(q[2], h ⊻ h_imagx) ⊻ hash(0, h ⊻ h_imagx) ⊻ hash(q[3], h ⊻ h_imagy) ⊻ hash(0, h ⊻ h_imagy) ⊻ hash(q[4], h ⊻ h_imagz) ⊻ hash(0, h ⊻ h_imagz))
    hash(q[1], h ⊻ hash(q[2], h_imagx) ⊻ hash_0_imagx ⊻ hash(q[3], h_imagy) ⊻ hash_0_imagy ⊻ hash(q[4], h_imagz) ⊻ hash_0_imagz)
end
# The scalar part of a `QuatVec` is hashed as a positive zero, matching `isequal`.
function Base.hash(q::QuatVec, h::UInt)
    hash(zero(q[2]), h ⊻ hash(q[2], h_imagx) ⊻ hash_0_imagx ⊻ hash(q[3], h_imagy) ⊻ hash_0_imagy ⊻ hash(q[4], h_imagz) ⊻ hash_0_imagz)
end

# These utility functions print a component `x` of a quaternion as a signed term, such as `"
# + 2.0"` or `" - 3.0"`, to which `show` appends the basis element.  A component whose
# printed form contains an operator is wrapped in parentheses.  In general, the sign must
# then stay inside the parentheses: moving the minus sign of the `Complex` number `-1 + 2im`
# outside would print `- (1 + 2im)`, which is `-1 - 2im`.  Only for a single signed real
# number, as printed by the standard real types, may the sign be moved outside, as in `-
# (3.0e-9)`.  Non-finite and `Bool` components are followed by `*`, as in `show` for
# `Complex`, so that `NaN*𝐢` and `true*𝐣` can be parsed back.  The methods taking `io`
# print the component with the properties of `io`, such as `:compact`.  Other types, such as
# symbolic types, are printed through the method taking only `x`, which extensions may
# specialize.
function _pm_ascii(s::AbstractString)
    if s[1] ∉ "+-"
        s = "+" * s
    end
    if occursin(r"[+^/-]", s[2:end])
        s = s[1] == '+' ? " + (" * s[2:end] * ")" : " + (" * s * ")"
    else
        s = " " * s[1] * " " * s[2:end]
    end
    s
end
_pm_ascii(x) = _pm_ascii(string(x))
_pm_ascii(io::IO, x) = _pm_ascii(x)
_pm_ascii(io::IO, x::Complex) = _pm_ascii(sprint(print, x; context=io))
function _pm_ascii(io::IO, x::Union{AbstractFloat,Integer,Rational})
    s = sprint(print, x; context=io)
    sign, s = s[1] == '-' ? ('-', s[2:end]) : ('+', s)
    if occursin(r"[+^/-]", s)
        s = "(" * s * ")"
    elseif x isa Bool || !isfinite(x)
        s = s * "*"
    end
    " " * sign * " " * s
end

function Base.show(io::IO, q::AbstractQuaternion)
    print(
        io,
        q isa QuatVec ? "" : q[1],
        _pm_ascii(io, q[2]), "𝐢",
        _pm_ascii(io, q[3]), "𝐣",
        _pm_ascii(io, q[4]), "𝐤"
    )
end

function Base.show(io::IO, q::Rotor)
    print(io, "rotor(")
    invoke(Base.show, Tuple{IO, AbstractQuaternion}, io, q)
    print(io, ")")
end


function Base.read(s::IO, QT::Type{Q}) where {T<:Number, Q<:AbstractQuaternion{T}}
    w = read(s,T)
    x = read(s,T)
    y = read(s,T)
    z = read(s,T)
    QT(w,x,y,z)
end

function Base.write(s::IO, q::AbstractQuaternion)
    write(s,q[1],q[2],q[3],q[4])
end

#Broadcast.broadcasted(f, q::QT, args...) where {QT<:AbstractQuaternion{<:Number}} = wrapper(QT)(f.(components(q), args...))
