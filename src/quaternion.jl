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
Quaternion{T}(q::Quaternion) where {T<:Real} = Quaternion{T}(q.components)
Quaternion{T}(w, x, y, z) where {T<:Real} = Quaternion(T.([w, x, y, z])...)


"""
    Quaternion(w, x, y, z)
    Quaternion(x, y, z)
    Quaternion(w)
    Quaternion{T}(w, x, y, z)

Creates a new quaternion with the given components.  The first argument `w` is
the scalar component, and `x`, `y`, and `z` are the corresponding "vector"
components.  The type of the returned quaternion will be inferred from the
input arguments.  If numeric arguments are missing, they will be set to zero.
The element type `T` can also be specified, by passing the type parameter as
usual.

Note that the constants [`imx`](@ref), [`imy`](@ref), and [`imz`](@ref) can
also be used like the complex `im` to create new `Quaternion` object.

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

"""
    imx

The quaternionic unit associated with rotation about the `x` axis.  Can also be entered as Unicode
bold: `𝐢`.

# Examples
```jldoctest
julia> imx * imx
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imx
0.0 + 1.2𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
const imx = Quaternion(false, true, false, false)
const 𝐢 = imx

"""
    imy

The quaternionic unit associated with rotation about the `y` axis.  Can also be entered as Unicode
bold: `𝐣`.

# Examples
```jldoctest
julia> imy * imy
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imy
0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤
```
"""
const imy = Quaternion(false, false, true, false)
const 𝐣 = imy

"""
    imz

The quaternionic unit associated with rotation about the `z` axis.  Can also be entered as Unicode
bold: `𝐤`.

# Examples
```jldoctest
julia> imz * imz
-1 + 0𝐢 + 0𝐣 + 0𝐤
julia> 1.2imz
0.0 + 0.0𝐢 + 0.0𝐣 + 1.2𝐤
```
"""
const imz = Quaternion(false, false, false, true)
const 𝐤 = imz

Base.zero(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(false, false, false, false)
Base.zero(q::Quaternion{T}) where {T<:Real} = Base.zero(Quaternion{T})
Base.one(::Type{Quaternion{T}}) where {T<:Real} = Quaternion{T}(true, false, false, false)
Base.one(q::Quaternion{T}) where {T<:Real} = Base.one(Quaternion{T})

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

@inline Base.getindex(q::Quaternion, i::Int) = (@boundscheck checkbounds(q.components,i); q.components[i])
# Base.getindex(q::Quaternion, i::Number) = q[convert(Int, i)]
Base.@propagate_inbounds Base.getindex(q::Quaternion, I) = [q[i] for i in I]
Base.eltype(::Type{Quaternion{T}}) where {T} = T

Base.promote_rule(::Type{Quaternion{T}}, ::Type{S}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.promote_rule(::Type{Quaternion{T}}, ::Type{Quaternion{S}}) where {T<:Real,S<:Real} =
    Quaternion{promote_type(T,S)}
Base.widen(::Type{Quaternion{T}}) where {T} = Quaternion{widen(T)}

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

Base.bswap(q::Quaternion) = Quaternion(bswap(q.w), bswap(q.x), bswap(q.y), bswap(q.z))

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
