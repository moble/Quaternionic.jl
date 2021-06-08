module Quaternionic
export Quaternion, randn, random_rotors
export as_quat_array, as_float_array, to_euler_phases!, to_euler_phases

using StaticArrays, Latexify, LaTeXStrings
import Random: AbstractRNG, default_rng
import Symbolics


abstract type AbstractQuaternion{T<:Real} <: Number end

struct Quaternion{T<:Real} <: AbstractQuaternion{T}
    components::SVector{4, T}
end

function Quaternion(w, x, y, z)
    Quaternion(SVector{4}(w, x, y, z))
end

function Quaternion(q)
    v = SVector{4}(q)
    Quaternion{eltype(v)}(v)
end

function Base.getproperty(q::Quaternion, sym::Symbol)
    if sym === :w
        return q.components[1]
    elseif sym === :x
        return q.components[2]
    elseif sym === :y
        return q.components[3]
    elseif sym === :z
        return q.components[4]
    else # fallback to getfield
        return getfield(q, sym)
    end
end

Base.copy(q::Quaternion) = Quaternion(Base.copy(q.components))

Base.getindex(q::Quaternion, i::Int) = q.components[i]
Base.getindex(q::Quaternion, i::Number) = q[convert(Int, i)]
Base.getindex(q::Quaternion, I) = [q[i] for i in I]

Base.iterate(q::Quaternion, state=1) = state > 4 ? nothing : (q[state], state+1)
Base.IteratorSize(::Type{Quaternion{T}}) where {T} = Base.HasLength()
Base.eltype(::Type{Quaternion{T}}) where {T} = T
Base.length(::Quaternion{T}) where {T} = 4

# function Base.setindex!(q::Quaternion, v, i::Int)
#     q.components[i] = v
# end

Base.firstindex(q::Quaternion) = 1
Base.lastindex(q::Quaternion) = 4

Base.:-(q::Quaternion) = Quaternion(-q.components)
Base.:+(q::Quaternion, p::Quaternion) = Quaternion(q.components+p.components)
Base.:-(q::Quaternion, p::Quaternion) = Quaternion(q.components-p.components)

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
Base.float(q::Quaternion{T}) where T<:Real = convert(Quaternion{float(T)}, q)
Base.real(::Type{Quaternion{T}}) where {T<:Real} = real(T)
Base.real(q::Quaternion) = q.w
Base.zero(::Type{Quaternion{T}}) where {T<:Real} = Quaternion(zero(T), zero(T), zero(T), zero(T))
Base.zero(q::Quaternion{T}) where {T<:Real} = Base.zero(Quaternion{T})
Base.one(::Type{Quaternion{T}}) where {T<:Real} = Quaternion(one(T), zero(T), zero(T), zero(T))
Base.one(q::Quaternion{T}) where {T<:Real} = Base.one(Quaternion{T})
Base.isfinite(q::Quaternion) = (
    isfinite(q.w)
    && isfinite(q.x)
    && isfinite(q.y)
    && isfinite(q.z)
)

Base.conj(q::Quaternion) = Quaternion(q.w, -q.x, -q.y, -q.z)
Base.abs2(q::Quaternion) = sum(q.components.^2)
Base.abs(q::Quaternion) = sqrt(abs2(q))
abs2vec(q::Quaternion) = sum(q.components[2:4].^2)
absvec(q::Quaternion) = sqrt(abs2vec(q))
Base.inv(q::Quaternion) = conj(q) / abs2(q)
# norm(q::Quaternion) = Base.abs2(q)  ## This might just be too confusing

function Base.log(q::Quaternion{T}) where {T}
    absolute2vec = abs2vec(q)
    if absolute2vec == zero(T)
        if q.w < 0
            return Quaternion(log(-q.w), zero(T), zero(T), T(π))
        end
        return Quaternion(log(q.w), zero(T), zero(T), zero(T))
    end
    absolutevec = sqrt(absolute2vec)
    f = atan(absolutevec, q.w) / absolutevec  # acos((w^2-absolutevec^2) / (w^2+absolutevec^2)) / 2absolutevec
    Quaternion(log(abs2(q))/2, f*q.x, f*q.y, f*q.z)
end

function Base.exp(q::Quaternion{T}) where {T}
    absolute2vec = abs2vec(q)
    if absolute2vec == zero(T)
        return Quaternion(exp(q.w), zero(T), zero(T), zero(T))
    end
    absolutevec = sqrt(absolute2vec)
    s = sin(absolutevec) / absolutevec
    e = exp(q.w)
    Quaternion(e*cos(absolutevec), e*s*q.x, e*s*q.y, e*s*q.z)
end

function Base.sqrt(q::Quaternion{T}) where {T}
    absolute2vec = abs2vec(q)
    if absolute2vec == zero(T)
        if q.w < 0
            return Quaternion(zero(T), zero(T), zero(T), sqrt(-q.w))
        end
        return Quaternion(sqrt(q.w), zero(T), zero(T), zero(T))
    end
    absolute2 = absolute2vec + q.w^2
    if absolute2 == zero(T)
        return zero(q)
    end
    c1 = sqrt(absolute2) + q.w
    c2 = sqrt(inv(2*c1))
    Quaternion(c1*c2, q.x*c2, q.y*c2, q.z*c2)
end

Base.angle(q::Quaternion) = 2 * absvec(log(q))


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


Broadcast.broadcasted(f, q::Quaternion, args...) = Quaternion(f.(q.components, args...))


include("conversion.jl")
include("random.jl")


end  # module
