# Useful functions from Base

function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:AbstractQuaternion,S<:AbstractQuaternion}
    Base.rtoldefault(eltype(x), eltype(y), atol)
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:AbstractQuaternion,S<:Number}
    Base.rtoldefault(eltype(x), y, atol)
end
function Base.rtoldefault(x::Union{T,Type{T}}, y::Union{S,Type{S}}, atol::Real) where {T<:Number,S<:AbstractQuaternion}
    Base.rtoldefault(x, eltype(y), atol)
end

function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1.w-q2.w; expand=true)) &&
        iszero(Symbolics.simplify(q1.x-q2.x; expand=true)) &&
        iszero(Symbolics.simplify(q1.y-q2.y; expand=true)) &&
        iszero(Symbolics.simplify(q1.z-q2.z; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Real)
    (
        iszero(Symbolics.simplify(q1.w-q2; expand=true)) &&
        iszero(Symbolics.simplify(q1.x; expand=true)) &&
        iszero(Symbolics.simplify(q1.y; expand=true)) &&
        iszero(Symbolics.simplify(q1.z; expand=true))
    )
end
function Base.:(==)(q1::Real, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1-q2.w; expand=true)) &&
        iszero(Symbolics.simplify(q2.x; expand=true)) &&
        iszero(Symbolics.simplify(q2.y; expand=true)) &&
        iszero(Symbolics.simplify(q2.z; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Symbolics.Num)
    (
        iszero(Symbolics.simplify(q1.w-q2; expand=true)) &&
        iszero(Symbolics.simplify(q1.x; expand=true)) &&
        iszero(Symbolics.simplify(q1.y; expand=true)) &&
        iszero(Symbolics.simplify(q1.z; expand=true))
    )
end
function Base.:(==)(q1::Symbolics.Num, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1-q2.w; expand=true)) &&
        iszero(Symbolics.simplify(q2.x; expand=true)) &&
        iszero(Symbolics.simplify(q2.y; expand=true)) &&
        iszero(Symbolics.simplify(q2.z; expand=true))
    )
end
Base.:(==)(q1::AbstractQuaternion{<:Real}, q2::AbstractQuaternion{<:Real}) = (q1.w==q2.w) && (q1.x==q2.x) && (q1.y==q2.y) && (q1.z==q2.z)
Base.:(==)(q::AbstractQuaternion{<:Real}, x::Real) = isreal(q) && real(q) == x
Base.:(==)(x::Real, q::AbstractQuaternion) = isreal(q) && real(q) == x

Base.isequal(q1::AbstractQuaternion, q2::AbstractQuaternion) = isequal(q1.w,q2.w) && isequal(q1.x,q2.x) && isequal(q1.y,q2.y) && isequal(q1.z,q2.z)

Base.isreal(q::AbstractQuaternion{T}) where {T<:Real} = iszero(q.x) && iszero(q.y) && iszero(q.z)
Base.isinteger(q::AbstractQuaternion{T}) where {T<:Real} = isreal(q) && isinteger(real(q))
Base.isfinite(q::AbstractQuaternion{T}) where {T<:Real} = isfinite(q.w) && isfinite(q.x) && isfinite(q.y) && isfinite(q.z)
Base.isnan(q::AbstractQuaternion{T}) where {T<:Real} = isnan(q.w) || isnan(q.x) || isnan(q.y) || isnan(q.z)
Base.isinf(q::AbstractQuaternion{T}) where {T<:Real} = isinf(q.w) || isinf(q.x) || isinf(q.y) || isinf(q.z)
Base.iszero(q::AbstractQuaternion{T}) where {T<:Real} = iszero(q.w) && iszero(q.x) && iszero(q.y) && iszero(q.z)
Base.isone(q::AbstractQuaternion{T}) where {T<:Real} = isone(q.w) && iszero(q.x) && iszero(q.y) && iszero(q.z)

Base.in(q::AbstractQuaternion, r::AbstractRange{<:Real}) = isreal(q) && real(q) in r

Base.flipsign(q::AbstractQuaternion, x::Real) = ifelse(signbit(x), -q, q)

Base.bswap(q::Q) where {T<:Real, Q<:AbstractQuaternion{T}} = Q(bswap(q.w), bswap(q.x), bswap(q.y), bswap(q.z))

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

function Base.hash(q::AbstractQuaternion, h::UInt)
    # TODO: with default argument specialization, this would be better:
    # hash(q.w, h ⊻ hash(q.x, h ⊻ h_imagx) ⊻ hash(0, h ⊻ h_imagx) ⊻ hash(q.y, h ⊻ h_imagy) ⊻ hash(0, h ⊻ h_imagy) ⊻ hash(q.z, h ⊻ h_imagz) ⊻ hash(0, h ⊻ h_imagz))
    hash(q.w, h ⊻ hash(q.x, h_imagx) ⊻ hash_0_imagx ⊻ hash(q.y, h_imagy) ⊻ hash_0_imagy ⊻ hash(q.z, h_imagz) ⊻ hash_0_imagz)
end

function Base.show(io::IO, q::AbstractQuaternion)
    function pm(x)
        s = "$x"
        if s[1] ∉ "+-"
            s = "+" * s
        end
        if occursin(r"[+^/-]", s[2:end])
            s = " " * s[1] * " " * "(" * s[2:end] * ")"
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

function Base.show(io::IO, ::MIME"text/latex", q::AbstractQuaternion)
    function pm(x)
        s = latexify(x, env=:raw, bracket=true)
        if s[1] ∉ "+-"
            s = "+" * s
        end
        if occursin(r"[+^/-]", s[2:end])
            s = " " * s[1] * " " * "\\left(" * s[2:end] * "\\right)"
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

function Base.read(s::IO, QT::Type{Q}) where {T<:Real, Q<:AbstractQuaternion{T}}
    w = read(s,T)
    x = read(s,T)
    y = read(s,T)
    z = read(s,T)
    QT(w,x,y,z)
end

function Base.write(s::IO, q::AbstractQuaternion)
    write(s,q.w,q.x,q.y,q.z)
end

#Broadcast.broadcasted(f, q::QT, args...) where {QT<:AbstractQuaternion{<:Real}} = wrapper(QT)(f.(q.components, args...))

# Broadcast-like operations from Symbolics
# (d::Symbolics.Operator)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q.w), d(q.x), d(q.y), d(q.z))
# (d::Symbolics.Operator)(q::QuatVec) = QuatVec(d(q.x), d(q.y), d(q.z))
(d::Symbolics.Differential)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q.w), d(q.x), d(q.y), d(q.z))
(d::Symbolics.Differential)(q::Rotor) = Quaternion(d(q.w), d(q.x), d(q.y), d(q.z))
(d::Symbolics.Differential)(q::QuatVec) = QuatVec(d(q.x), d(q.y), d(q.z))
