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

Base.:(==)(q1::AbstractQuaternion{<:Number}, q2::AbstractQuaternion{<:Number}) = (q1[1]==q2[1]) && (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q::AbstractQuaternion{<:Number}, x::Number) = isreal(q) && real(q) == x
Base.:(==)(x::Number, q::AbstractQuaternion) = isreal(q) && real(q) == x
Base.isequal(q1::AbstractQuaternion, q2::AbstractQuaternion) = isequal(q1[1],q2[1]) && isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])

Base.:(==)(q1::QuatVec{<:Number}, q2::AbstractQuaternion{<:Number}) = (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q1::AbstractQuaternion{<:Number}, q2::QuatVec{<:Number}) = (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q1::QuatVec{<:Number}, q2::QuatVec{<:Number}) = (q1[2]==q2[2]) && (q1[3]==q2[3]) && (q1[4]==q2[4])
Base.:(==)(q::QuatVec{<:Number}, x::Number) = false
Base.:(==)(x::Number, q::QuatVec) = false
Base.isequal(q1::AbstractQuaternion, q2::QuatVec) = isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])
Base.isequal(q1::QuatVec, q2::AbstractQuaternion) = isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])
Base.isequal(q1::QuatVec, q2::QuatVec) = isequal(q1[2],q2[2]) && isequal(q1[3],q2[3]) && isequal(q1[4],q2[4])

Base.isreal(q::AbstractQuaternion{T}) where {T<:Number} = iszero(q[2]) && iszero(q[3]) && iszero(q[4])
Base.isinteger(q::AbstractQuaternion{T}) where {T<:Number} = isreal(q) && isinteger(real(q))
Base.isfinite(q::AbstractQuaternion{T}) where {T<:Number} = isfinite(q[1]) && isfinite(q[2]) && isfinite(q[3]) && isfinite(q[4])
Base.isnan(q::AbstractQuaternion{T}) where {T<:Number} = isnan(q[1]) || isnan(q[2]) || isnan(q[3]) || isnan(q[4])
Base.isinf(q::AbstractQuaternion{T}) where {T<:Number} = isinf(q[1]) || isinf(q[2]) || isinf(q[3]) || isinf(q[4])
Base.iszero(q::AbstractQuaternion{T}) where {T<:Number} = iszero(q[1]) && iszero(q[2]) && iszero(q[3]) && iszero(q[4])
Base.isone(q::AbstractQuaternion{T}) where {T<:Number} = isone(q[1]) && iszero(q[2]) && iszero(q[3]) && iszero(q[4])

Base.round(q::QT, r::RoundingMode=RoundNearest; kwargs...) where {QT<:AbstractQuaternion} = QT(round.(components(q), r; kwargs...))

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

function Base.hash(q::AbstractQuaternion, h::UInt)
    # TODO: with default argument specialization, this would be better:
    # hash(q[1], h ⊻ hash(q[2], h ⊻ h_imagx) ⊻ hash(0, h ⊻ h_imagx) ⊻ hash(q[3], h ⊻ h_imagy) ⊻ hash(0, h ⊻ h_imagy) ⊻ hash(q[4], h ⊻ h_imagz) ⊻ hash(0, h ⊻ h_imagz))
    hash(q[1], h ⊻ hash(q[2], h_imagx) ⊻ hash_0_imagx ⊻ hash(q[3], h_imagy) ⊻ hash_0_imagy ⊻ hash(q[4], h_imagz) ⊻ hash_0_imagz)
end

function _pm_ascii(x)
    # Utility function to print a component of a quaternion
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

function Base.show(io::IO, q::AbstractQuaternion)
    print(
        io,
        q isa QuatVec ? "" : q[1],
        _pm_ascii(q[2]), "𝐢",
        _pm_ascii(q[3]), "𝐣",
        _pm_ascii(q[4]), "𝐤"
    )
end

function Base.show(io::IO, q::Rotor)
    print(io, "rotor(")
    invoke(Base.show, Tuple{IO, AbstractQuaternion}, io, q)
    print(io, ")")
end

_pm_latex(x) = _pm_ascii(x)


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
