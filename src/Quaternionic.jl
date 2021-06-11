module Quaternionic

export Quaternion, QuaternionF64, QuaternionF32, QuaternionF16, imx, imy, imz
export abs2vec, absvec
export randn_rotor
export as_quat_array, as_float_array, to_euler_phases!, to_euler_phases

using StaticArrays, Latexify, LaTeXStrings
import Random: AbstractRNG, default_rng, randn!
import Symbolics


abstract type AbstractQuaternion{T<:Real} <: Number end

"""
    Quaternion{T<:Real} <: Number

Quaternionic number type with elements of type `T`.

`QuaternionF16`, `QuaternionF32` and `QuaternionF64` are aliases for `Quaternion{Float16}`,
`Quaternion{Float32}` and `Quaternion{Float64}` respectively.

See also: [`Quaternion`](@ref)
"""
struct Quaternion{T<:Real} <: AbstractQuaternion{T}
    components::SVector{4, T}
end
function Quaternion{T}(w, x, y, z) where {T<:Real}
    Quaternion(T.([w, x, y, z])...)
end

"""
    Quaternion(w, x, y, z)
    Quaternion(x, y, z)
    Quaternion(w)
    Quaternion(:z)
    Quaternion{T}(w, x, y, z)

Creates a new quaternion with the given components.  The first argument `w` is
the scalar component, and `x`, `y`, and `z` are the corresponding "vector"
components.  The type of the returned quaternion will be inferred from the
input arguments.  If numeric arguments are missing, they will be set to zero.
It is also possible to pass one of the symbols `:w`, `:x`, `:y`, or `:z` to
obtain a unit vector (by default, with eltype `Float64`) along the corresponding
direction.  It is also possible to specify the element type `T`, by passing the
type parameter as usual.

# Examples

```jldoctest
julia> Quaternion(1, 2, 3, 4)
1 + 2𝐢 + 3𝐣 + 4𝐤
julia> Quaternion{Float64}(1, 2, 3, 4)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> Quaternion(1.0, 2.0, 3.0, 4.0)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> Quaternion(2, 3, 4)
0 + 2𝐢 + 3𝐣 + 4𝐤
julia> Quaternion(1)
1 + 0𝐢 + 0𝐣 + 0𝐤
julia> Quaternion(:z)
0.0 + 0.0𝐢 + 0.0𝐣 + 1.0𝐤

```
"""
function Quaternion(w::Real, x::Real, y::Real, z::Real)
    Quaternion(SVector{4}(w, x, y, z))
end
function Quaternion(w::Real)
    Quaternion(w, zero(w), zero(w), zero(w))
end
function Quaternion(x::Real, y::Real, z::Real)
    Quaternion(zero(x), x, y, z)
end
function Quaternion{T}(sym::Symbol) where {T<:Real}
    if sym === :w
        return Quaternion{T}(one(T), zero(T), zero(T), zero(T))
    elseif sym === :x
        return Quaternion{T}(zero(T), one(T), zero(T), zero(T))
    elseif sym === :y
        return Quaternion{T}(zero(T), zero(T), one(T), zero(T))
    elseif sym === :z
        return Quaternion{T}(zero(T), zero(T), zero(T), one(T))
    else
        throw(ArgumentError(
            "Only :w, :x, :y, or :z are accepted, not `$sym`"
        ))
    end
end
Quaternion(sym::Symbol) = Quaternion{Float64}(sym)
Quaternion(q::Quaternion) = q
function Quaternion(q)
    v = SVector{4}(q)
    Quaternion{eltype(v)}(v)
end

Quaternion(::Type{T}) where {T<:Real} = Quaternion{T}
Quaternion(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}

const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}

Base.zero(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(false, false, false, false)
Base.zero(q::Quaternion{T}) where {T<:Real} = Base.zero(Quaternion{T})
Base.one(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(true, false, false, false)
Base.one(q::Quaternion{T}) where {T<:Real} = Base.one(Quaternion{T})
"""
    imx

The quaternionic unit associated with rotation about the `x` axis.

# Examples
```jldoctest
julia> imx * imx
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imx
0.0 + 1.2𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
const imx = Quaternion(false, true, false, false)
"""
    imy

The quaternionic unit associated with rotation about the `y` axis.

# Examples
```jldoctest
julia> imy * imy
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imy
0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤
```
"""
const imy = Quaternion(false, false, true, false)
"""
    imz

The quaternionic unit associated with rotation about the `z` axis.

# Examples
```jldoctest
julia> imz * imz
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imz
0.0 + 0.0𝐢 + 0.0𝐣 + 1.2𝐤
```
"""
const imz = Quaternion(false, false, false, true)

function Base.getproperty(q::Quaternion, sym::Symbol)
    @inbounds begin
        if sym === :w
            return q.components[1]
        elseif sym === :x
            return q.components[2]
        elseif sym === :y
            return q.components[3]
        elseif sym === :z
            return q.components[4]
        elseif sym === :re
            return q.components[1]
        elseif sym === :im
            return q.components[2:4]
        elseif sym === :vec
            return q.components[2:4]
        else # fallback to getfield
            return getfield(q, sym)
        end
    end
end

Base.getindex(q::Quaternion, i::Int) = q.components[i]
# Base.getindex(q::Quaternion, i::Number) = q[convert(Int, i)]
# Base.getindex(q::Quaternion, I) = [q[i] for i in I]
Base.eltype(::Type{Quaternion{T}}) where {T} = T

# COV_EXCL_START
Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}
Base.bswap(q::Quaternion) = Quaternion(bswap(q.w), bswap(q.x), bswap(q.y), bswap(q.z))
# COV_EXCL_STOP

Base.float(::Type{Quaternion{T}}) where {T<:AbstractFloat} = Quaternion{T}
Base.float(::Type{Quaternion{T}}) where {T} = Quaternion{float(T)}
Base.float(q::Quaternion{T}) where T<:Real = Quaternion(float(q.components))

Base.real(::Type{Quaternion{T}}) where {T<:Real} = real(T)
Base.real(q::Quaternion) = q.re
Base.imag(q::Quaternion) = q.im

Base.isreal(q::Quaternion) = iszero(q.x) && iszero(q.y) && iszero(q.z)
Base.isinteger(q::Quaternion) = isreal(q) && isinteger(real(q))
Base.isfinite(q::Quaternion) = isfinite(q.w) && isfinite(q.x) && isfinite(q.y) && isfinite(q.z)
Base.isnan(q::Quaternion) = isnan(q.w) || isnan(q.x) || isnan(q.y) || isnan(q.z)
Base.isinf(q::Quaternion) = isinf(q.w) || isinf(q.x) || isinf(q.y) || isinf(q.z)
Base.iszero(q::Quaternion) = iszero(q.w) && iszero(q.x) && iszero(q.y) && iszero(q.z)
Base.isone(q::Quaternion) = isone(q.w) && iszero(q.x) && iszero(q.y) && iszero(q.z)

Base.:-(q::Quaternion) = Quaternion(-q.components)

Base.:+(q::Quaternion, p::Quaternion) = Quaternion(q.components+p.components)
Base.:-(q::Quaternion, p::Quaternion) = Quaternion(q.components-p.components)

Base.:+(q::Quaternion, p::Number) = Quaternion(q.w+p, q.x, q.y, q.z)
Base.:-(q::Quaternion, p::Number) = Quaternion(q.w-p, q.x, q.y, q.z)

Base.:+(q::Number, p::Quaternion) = Quaternion(q+p.w, p.x, p.y, p.z)
Base.:-(q::Number, p::Quaternion) = Quaternion(q-p.w, -p.x, -p.y, -p.z)

Base.flipsign(q::Quaternion, x::Real) = ifelse(signbit(x), -q, q)

function Base.:*(q::Quaternion, p::Quaternion)
    Quaternion(
        q.w*p.w - q.x*p.x - q.y*p.y - q.z*p.z,
        q.w*p.x + q.x*p.w + q.y*p.z - q.z*p.y,
        q.w*p.y - q.x*p.z + q.y*p.w + q.z*p.x,
        q.w*p.z + q.x*p.y - q.y*p.x + q.z*p.w
    )
end

function Base.:/(q::Quaternion, p::Quaternion)
    pnorm = p.w^2 + p.x^2 + p.y^2 + p.z^2
    Quaternion(
        (+q.w*p.w + q.x*p.x + q.y*p.y + q.z*p.z) / pnorm,
        (-q.w*p.x + q.x*p.w - q.y*p.z + q.z*p.y) / pnorm,
        (-q.w*p.y + q.x*p.z + q.y*p.w - q.z*p.x) / pnorm,
        (-q.w*p.z - q.x*p.y + q.y*p.x + q.z*p.w) / pnorm
    )
end

Base.:*(s::Number, p::Quaternion) = Quaternion(s*p.components)

function Base.:/(s::Number, p::Quaternion)
    pnorm = s / (p.w^2 + p.x^2 + p.y^2 + p.z^2)
    Quaternion(
        p.w * pnorm,
        -p.x * pnorm,
        -p.y * pnorm,
        -p.z * pnorm
    )
end

function Base.:*(q::Quaternion, s::Number)
    Quaternion(s*q.components)
end

function Base.:/(q::Quaternion, s::Number)
    Quaternion(q.components / s)
end

Base.:*(s::Symbolics.Num, p::Quaternion) = Quaternion(s*p.components)

function Base.:/(s::Symbolics.Num, p::Quaternion)
    pnorm = s / (p.w^2 + p.x^2 + p.y^2 + p.z^2)
    Quaternion(
        p.w * pnorm,
        -p.x * pnorm,
        -p.y * pnorm,
        -p.z * pnorm
    )
end

function Base.:*(q::Quaternion, s::Symbolics.Num)
    Quaternion(s*q.components)
end

function Base.:/(q::Quaternion, s::Symbolics.Num)
    Quaternion(q.components / s)
end

function Base.:(==)(q1::Quaternion{Symbolics.Num}, q2::Quaternion{Symbolics.Num})
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
Base.:(==)(q1::Quaternion, q2::Quaternion) = (q1.w==q2.w) && (q1.x==q2.x) && (q1.y==q2.y) && (q1.z==q2.z)
Base.:(==)(q::Quaternion, x::Real) = isreal(q) && real(q) == x
Base.:(==)(x::Real, q::Quaternion) = isreal(q) && real(q) == x
Base.:(==)(q::Quaternion, x::Symbolics.Num) = isreal(q) && real(q) == x
Base.:(==)(x::Symbolics.Num, q::Quaternion) = isreal(q) && real(q) == x

Base.isequal(q1::Quaternion, q2::Quaternion) = isequal(q1.w,q2.w) && isequal(q1.x,q2.x) && isequal(q1.y,q2.y) && isequal(q1.z,q2.z)
Base.in(q::Quaternion, r::AbstractRange{<:Real}) = isreal(q) && real(q) in r

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

function Base.hash(q::Quaternion, h::UInt)
    # TODO: with default argument specialization, this would be better:
    # hash(q.w, h ⊻ hash(q.x, h ⊻ h_imagx) ⊻ hash(0, h ⊻ h_imagx) ⊻ hash(q.y, h ⊻ h_imagy) ⊻ hash(0, h ⊻ h_imagy) ⊻ hash(q.z, h ⊻ h_imagz) ⊻ hash(0, h ⊻ h_imagz))
    hash(q.w, h ⊻ hash(q.x, h_imagx) ⊻ hash_0_imagx ⊻ hash(q.y, h_imagy) ⊻ hash_0_imagy ⊻ hash(q.z, h_imagz) ⊻ hash_0_imagz)
end

"""
    conj(q)

Return the quaternion conjugate, which flips the sign of each "vector"
component.

# Examples
```jldoctest
julia> conj(Quaternion(1,2,3,4))
1 - 2𝐢 - 3𝐣 - 4𝐤
```
"""
Base.conj(q::Quaternion) = Quaternion(q.w, -q.x, -q.y, -q.z)

"""
    abs2(q)

Sum the squares of the components of the quaternion

# Examples
```jldoctest
julia> abs2(Quaternion(1,2,4,10))
121
```
"""
Base.abs2(q::Quaternion) = sum(q.components.^2)

"""
    abs(q)

Square-root of the sum the squares of the components of the quaternion

# Examples
```jldoctest
julia> abs(Quaternion(1,2,4,10))
11.0
```
"""
Base.abs(q::Quaternion) = sqrt(abs2(q))

"""
    abs2vec(q)

Sum the squares of the "vector" components of the quaternion

# Examples
```jldoctest
julia> abs2vec(Quaternion(1,2,3,6))
49
```
"""
abs2vec(q::Quaternion) = @inbounds q.components[2]^2 + q.components[3]^2 + q.components[4]^2

"""
    absvec(q)

Square-root of the sum of the squares of the "vector" components of the quaternion

# Examples
```jldoctest
julia> absvec(Quaternion(1,2,3,6))
7.0
```
"""
absvec(q::Quaternion) = sqrt(abs2vec(q))
# norm(q::Quaternion) = Base.abs2(q)  ## This might just be confusing
Base.inv(q::Quaternion) = conj(q) / abs2(q)

"""
    log(q)

Logarithm of a quaternion.

As with the usual complex logarithm, the quaternion logarithm has multiple
branches, though the quaternion branches are three-dimensional: for any unit
"vector" quaternion q̂, you could add any integer multiple of 2πq̂ to the result
of this function and still get the same result after exponentiating (within
numerical accuracy).  This function is the principal logarithm.

This function has discontinuous (and fairly arbitrary) behavior along the
negative real axis: if the "vector" components of the quaternion are precisely
zero *and* the scalar component is negative, the returned quaternion will have
scalar component `log(-q.w)`, but will also have a `z` component of π.  The
choice of the `z` direction is arbitrary; the "vector" component of the
returned quaternion could be π times any unit vector.

Note that this function is not specialized to unit-quaternion inputs, so the
scalar component of the returned value will be nonzero unless the input has
*precisely* unit magnitude.

# Examples
```jldoctest
julia> log(exp(1.2Quaternion(:y)))
0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤

julia> log(Quaternion(exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> log(Quaternion(-exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 3.141592653589793𝐤
```

"""
function Base.log(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        if q.w < 0
            return Quaternion(log(-q.w), 0, 0, π)
        end
        return Quaternion(log(q.w), 0, 0, 0)
    end
    absolutevec = sqrt(absolute2vec)
    f = atan(absolutevec, q.w) / absolutevec  # acos((w^2-absolutevec^2) / (w^2+absolutevec^2)) / 2absolutevec
    Quaternion(log(abs2(q))/2, f*q.x, f*q.y, f*q.z)
end
# function Base.log(q::UnitQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if iszero(absolute2vec)
#         if q.w < 0
#             return Quaternion{T}(0, 0, 0, π)
#         end
#         return Quaternion{T}(0, 0, 0, 0)
#     end
#     absolutevec = sqrt(absolute2vec)
#     f = atan(absolutevec, q.w) / absolutevec  # acos((w^2-absolutevec^2) / (w^2+absolutevec^2)) / 2absolutevec
#     Quaternion(0, f*q.x, f*q.y, f*q.z)
# end

"""
    exp(q)

Exponential of a quaternion

# Examples
```jldoctest
julia> exp(π/4*Quaternion(:x))  # Rotation by π/2 about the x axis
0.7071067811865476 + 0.7071067811865475𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
function Base.exp(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        return Quaternion(exp(q.w), 0, 0, 0)
    end
    absolutevec = sqrt(absolute2vec)
    s = sin(absolutevec) / absolutevec
    e = exp(q.w)
    Quaternion(e*cos(absolutevec), e*s*q.x, e*s*q.y, e*s*q.z)
end
# function Base.exp(q::VectorQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if iszero(absolute2vec)
#         return Quaternion{T}(0, 0, 0, 0)
#     end
#     absolutevec = sqrt(absolute2vec)
#     s = sin(absolutevec) / absolutevec
#     Quaternion(cos(absolutevec), s*q.x, s*q.y, s*q.z)
# end

@doc raw"""
    sqrt(q)

Square-root of a quaternion.

The general formula whenever the denominator is nonzero is

``
\sqrt{q} = \frac{|q| + q} {\sqrt{2|q| + 2q.w}}
``

This can be proven by expanding `q` as `q.w + q.vec` and multiplying the
expression above by itself.

When the denominator is zero, this function has discontinuous (and fairly
arbitrary) behavior, just as with the quaternion [`log`](@ref) function.  In
this case, either all components are zero — in which case the result is simply
the zero quaternion — or the "vector" components of the quaternion are
precisely zero and the scalar component is negative.  If the latter is true,
the denominator above will be a pure-imaginary number.  Because the quaternions
come with infinitely many elements that square to -1, it is not clear *which*
imaginary should be used, so we arbitrarily choose to set the result
proportional to the `z` quaternion.  The choice of the `z` direction is
arbitrary; the "vector" component of the returned quaternion could be in any
direction.

# Examples
```jldoctest
julia> q = Quaternion(1.2, 3.4, 5.6, 7.8);

julia> sqrtq = √q;

julia> sqrtq^2 ≈ q
true

julia> √Quaternion(4)
2.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> √Quaternion(-4)
0.0 + 0.0𝐢 + 0.0𝐣 + 2.0𝐤
```
"""
function Base.sqrt(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        if q.w < 0
            return Quaternion(0, 0, 0, sqrt(-q.w))
        end
        return Quaternion(sqrt(q.w), 0, 0, 0)
    end
    absolute2 = absolute2vec + q.w^2
    c1 = sqrt(absolute2) + q.w
    c2 = sqrt(inv(2*c1))
    Quaternion(c1*c2, q.x*c2, q.y*c2, q.z*c2)
end
# function Base.sqrt(q::UnitQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if absolute2vec == zero(T)
#         if q.w < 0
#             return Quaternion{T}(0, 0, 0, 1)
#         end
#         return Quaternion{T}(1, 0, 0, 0)
#     end
#     c1 = 1 + q.w
#     c2 = sqrt(inv(2*c1))
#     Quaternion(c1*c2, q.x*c2, q.y*c2, q.z*c2)
# end

"""
    angle(q)

Phase angle in radians of the rotation represented by this quaternion.

Note that this may be different from your interpretation of the angle of a
complex number in an important way.  Because quaternions act on vectors by
conjugation — as in `q*v*conj(q)` — there are *two* copies of `q` involved in
that expression; in some sense, a quaternion acts "twice".  Therefore, this
angle may be twice what you expect from an analogy with complex numbers —
dpending on how you interpret the correspondence between complex numbers and
quaternions.  Also, while rotations in the complex plane have a natural choice
of axis (the positive `z` direction), that is not the case for quaternions,
which means that the sign of this angle is arbitrary, and we always choose it
to be positive.

# Examples
```jldoctest
julia> θ=1.2;

julia> R=exp(θ * Quaternion(:z) / 2);

julia> angle(R)
1.2

```
"""
Base.angle(q::Quaternion) = 2 * absvec(log(q))
# Base.angle(q::UnitQuaternion) = 2 * absvec(log(q))


function Base.show(io::IO, q::Quaternion)
    function pm(x)
        s = "$x"
        if s[1] ∉ "+-"
            s = "+" * s
        end
        if occursin(r"[+-]", s[2:end])
            s = " " * s[1] * " " * "{" * s[2:end] * "}"
        else
            s = " " * s[1] * " " * s[2:end]
        end
        s
    end
    print(
        io,
        q.w,
        pm(q.x), "𝐢",
        pm(q.y), "𝐣",
        pm(q.z), "𝐤"
    )
end

function Base.show(io::IO, ::MIME"text/latex", q::Quaternion)
    function pm(x)
        s = latexify(x, env=:raw, bracket=true)
        if s[1] ∉ "+-"
            s = "+" * s
        end
        if occursin(r"[+-]", s[2:end])
            s = " " * s[1] * " " * "\\left\\{" * s[2:end] * "\\right\\}"
        else
            s = " " * s[1] * " " * s[2:end]
        end
        s
    end
    s = latexstring(
        latexify(q.w, env=:raw, bracket=true),
        pm(q.x), "\\,\\mathbf{i}",
        pm(q.y), "\\,\\mathbf{j}",
        pm(q.z), "\\,\\mathbf{k}"
    )
    print(io, s)
end

function Base.read(s::IO, ::Type{Quaternion{T}}) where T<:Real
    w = read(s,T)
    x = read(s,T)
    y = read(s,T)
    z = read(s,T)
    Quaternion{T}(w,x,y,z)
end
function Base.write(s::IO, q::Quaternion)
    write(s,q.w,q.x,q.y,q.z)
end


Broadcast.broadcasted(f, q::Quaternion, args...; kwargs...) = Quaternion(f.(q.components, args...; kwargs...))


include("conversion.jl")
include("random.jl")


end  # module
