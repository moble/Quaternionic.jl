"""
    Quaternion{T<:Real} <: Number

Quaternionic number type with elements of type `T`.

`QuaternionF16`, `QuaternionF32` and `QuaternionF64` are aliases for `Quaternion{Float16}`,
`Quaternion{Float32}` and `Quaternion{Float64}` respectively.  See also [`Rotor`](@ref) and
[`QuatVec`](@ref).

The functions

    Quaternion(w, x, y, z)
    Quaternion(x, y, z)
    Quaternion(w)
    Quaternion{T}(w, x, y, z)

create a new quaternion with the given components.  The argument `w` is the
scalar component, and `x`, `y`, and `z` are the corresponding "vector"
components.  If any of these arguments is missing, it will be set to zero.  The
type of the returned quaternion will be inferred from the input arguments, or
can be specified, by passing the type parameter `T` as above.

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
struct Quaternion{T<:Real} <: AbstractQuaternion{T}
    components::SVector{4, T}
end

struct Rotor{T<:Real} <: AbstractQuaternion{T}
    components::SVector{4, T}
end

struct QuatVec{T<:Real} <: AbstractQuaternion{T}
    components::SVector{4, T}
end


# Constructor from all 4 components
(::Type{QT})(w, x, y, z) where {QT<:AbstractQuaternion} = QT(SVector{4}(w, x, y, z))
(::Type{QT})(w, x, y, z) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(w, x, y, z))
Rotor(w, x, y, z) = (n=√(w^2+x^2+y^2+z^2); Rotor(SVector{4}(w/n, x/n, y/n, z/n)))
Rotor{T}(w, x, y, z) where {T<:Real} = (n=√(w^2+x^2+y^2+z^2); Rotor{T}(SVector{4, T}(w/n, x/n, y/n, z/n)))
QuatVec(w, x, y, z) = QuatVec(SVector{4}(false, x, y, z))
QuatVec{T}(w, x, y, z) where {T<:Real} = QuatVec{T}(SVector{4, T}(false, x, y, z))

# Constructor from vector components
(::Type{QT})(x, y, z) where {QT<:AbstractQuaternion} = QT(SVector{4}(false, x, y, z))
(::Type{QT})(x, y, z) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(false, x, y, z))
Rotor(x, y, z) = (n=√(x^2+y^2+z^2); Rotor(SVector{4}(false, x/n, y/n, z/n)))
Rotor{T}(x, y, z) where {T<:Real} = (n=√(x^2+y^2+z^2); Rotor{T}(SVector{4, T}(false, x/n, y/n, z/n)))
QuatVec(x, y, z) = QuatVec(SVector{4}(false, x, y, z))
QuatVec{T}(x, y, z) where {T<:Real} = QuatVec{T}(SVector{4, T}(false, x, y, z))

# Constructor from scalar component
(::Type{QT})(w::Real) where {QT<:AbstractQuaternion} = QT(SVector{4}(w, false, false, false))
(::Type{QT})(w::Real) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(w, false, false, false))
Quaternion(w::Real) = Quaternion(SVector{4}(w, false, false, false))
Quaternion{T}(w::Real) where {T<:Real} = Quaternion{T}(SVector{4, T}(w, false, false, false))
Rotor(w::Real) = Rotor(SVector{4}(one(w), false, false, false))
Rotor{T}(w::Real) where {T<:Real} = Rotor{T}(SVector{4, T}(one(T), false, false, false))
QuatVec(w::Real) = QuatVec(SVector{4, typeof(w)}(false, false, false, false))
QuatVec{T}(w::Real) where {T<:Real} = QuatVec{T}(SVector{4, T}(false, false, false, false))

# Copy constructor
(::Type{QT})(q::AbstractQuaternion) where {QT<:AbstractQuaternion} = QT(q.components)
(::Type{QT})(q::AbstractQuaternion{S}) where {T<:Real, S<:Real, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(q.components))
(::Type{QT})(q::AbstractQuaternion{T}) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(SVector{4, T}(q.components))
Quaternion(q::AbstractQuaternion{T}) where {T<:Real} = Quaternion(q.components...)
Quaternion{T}(q::AbstractQuaternion{S}) where {T<:Real, S<:Real} = Quaternion{T}(q.components...)
Rotor(q::AbstractQuaternion{T}) where {T<:Real} = Rotor(q.components...)
Rotor{T}(q::AbstractQuaternion{S}) where {T<:Real, S<:Real} = Rotor{T}(q.components...)
QuatVec(q::AbstractQuaternion{T}) where {T<:Real} = QuatVec(q.components...)
QuatVec{T}(q::AbstractQuaternion{S}) where {T<:Real, S<:Real} = QuatVec{T}(q.components...)

# Type constructors
(::Type{QT})(::Type{T}) where {T<:Real, QT<:AbstractQuaternion} = QT{T}
(::Type{QT})(::Type{<:AbstractQuaternion{T}}) where {T<:Real, QT<:AbstractQuaternion} = QT{T}

# Handy aliases like `ComplexF64`, etc.
const QuaternionF64 = Quaternion{Float64}
const QuaternionF32 = Quaternion{Float32}
const QuaternionF16 = Quaternion{Float16}
const RotorF64 = Rotor{Float64}
const RotorF32 = Rotor{Float32}
const RotorF16 = Rotor{Float16}
const QuatVecF64 = QuatVec{Float64}
const QuatVecF32 = QuatVec{Float32}
const QuatVecF16 = QuatVec{Float16}

# Handy constants like `im`
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

# Essential constructors
Base.zero(::Type{QT}) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(false, false, false, false)
Base.zero(q::QT) where {T<:Real, QT<:AbstractQuaternion{T}} = Base.zero(QT)
Base.zero(::Type{Rotor}) = throw(DomainError("Rotor", "Zero is not a possible rotor."))
Base.zero(::Type{Rotor{T}}) where T = throw(DomainError("Rotor", "Zero is not a possible rotor."))

Base.one(::Type{QT}) where {T<:Real, QT<:AbstractQuaternion{T}} = QT(true, false, false, false)
Base.one(q::QT) where {T<:Real, QT<:AbstractQuaternion{T}} = Base.one(QT)
Base.one(::Type{QuatVec}) = throw(DomainError("QuatVec", "One is not a possible 3-vector."))
Base.one(::Type{QuatVec{T}}) where T = throw(DomainError("QuatVec", "One is not a possible 3-vector."))

# Getting pieces of quaternions
function Base.getproperty(q::AbstractQuaternion, sym::Symbol)
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
@inline Base.getindex(q::AbstractQuaternion, i::Int) = (@boundscheck checkbounds(q.components,i); q.components[i])
# Base.getindex(q::AbstractQuaternion, i::Number) = q[convert(Int, i)]
Base.@propagate_inbounds Base.getindex(q::AbstractQuaternion, I) = [q[i] for i in I]
Base.real(::Type{AbstractQuaternion{T}}) where {T<:Real} = real(T)
Base.real(q::AbstractQuaternion{T}) where {T<:Real} = q.re
Base.imag(q::AbstractQuaternion{T}) where {T<:Real} = q.im

# Type games
Base.eltype(::Type{<:AbstractQuaternion{T}}) where {T} = T
Base.widen(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){widen(eltype(Q))}
Base.float(::Type{Q}) where {Q<:AbstractQuaternion{<:AbstractFloat}} = Q
Base.float(::Type{Q}) where {Q<:AbstractQuaternion} = wrapper(Q){float(eltype(Q))}
Base.float(q::AbstractQuaternion{T}) where {T<:AbstractFloat} = q
Base.float(q::AbstractQuaternion{T}) where {T} = wrapper(typeof(q)){float(T)}(float(q.components))

Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Real} =
    wrapper(Q){promote_type(eltype(Q),S)}
Base.promote_rule(::Type{Q1}, ::Type{Q2}) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion} =
    wrapper(Q1, Q2){promote_type(eltype(Q1),eltype(Q2))}

wrapper(::T) where {T} = wrapper(T)

wrapper(::Type{<:AbstractQuaternion}, ::Val{OP}, ::Type{<:AbstractQuaternion}) where {OP} = Quaternion
wrapper(::Type{<:QuatVec}, ::Val{+}, ::Type{<:QuatVec}) = QuatVec
wrapper(::Type{<:QuatVec}, ::Val{-}, ::Type{<:QuatVec}) = QuatVec
wrapper(::Type{<:Rotor}, ::Val{*}, ::Type{<:Rotor}) = Rotor
wrapper(::Type{<:Rotor}, ::Val{/}, ::Type{<:Rotor}) = Rotor
