Base.:-(q::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(-q.components)

Base.:+(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion} = wrapper(Q1, Q2)(q.components+p.components)
Base.:-(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion} = wrapper(Q1, Q2)(q.components-p.components)

Base.:+(q::Q, p::Number) where {Q<:AbstractQuaternion} = wrapper(Q)(q.w+p, q.x, q.y, q.z)
Base.:-(q::Q, p::Number) where {Q<:AbstractQuaternion} = wrapper(Q)(q.w-p, q.x, q.y, q.z)

Base.:+(q::Number, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(q+p.w, p.x, p.y, p.z)
Base.:-(q::Number, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(q-p.w, -p.x, -p.y, -p.z)

Base.:+(q::Q, p::Symbolics.Num) where {Q<:AbstractQuaternion} = wrapper(Q)(q.w+p, q.x, q.y, q.z)
Base.:-(q::Q, p::Symbolics.Num) where {Q<:AbstractQuaternion} = wrapper(Q)(q.w-p, q.x, q.y, q.z)

Base.:+(q::Symbolics.Num, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(q+p.w, p.x, p.y, p.z)
Base.:-(q::Symbolics.Num, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(q-p.w, -p.x, -p.y, -p.z)

Base.flipsign(q::AbstractQuaternion, x::Real) = ifelse(signbit(x), -q, q)

function Base.:*(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    wrapper(Q1, Q2)(
        q.w*p.w - q.x*p.x - q.y*p.y - q.z*p.z,
        q.w*p.x + q.x*p.w + q.y*p.z - q.z*p.y,
        q.w*p.y - q.x*p.z + q.y*p.w + q.z*p.x,
        q.w*p.z + q.x*p.y - q.y*p.x + q.z*p.w
    )
end

function Base.:/(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    if p == q
        return one(promote_type(typeof(q), float(typeof(p))))
    end
    den = (p.w^2 + p.x^2 + p.y^2 + p.z^2)
    wrapper(Q1, Q2)(
        (+q.w*p.w + q.x*p.x + q.y*p.y + q.z*p.z) / den,
        (-q.w*p.x + q.x*p.w - q.y*p.z + q.z*p.y) / den,
        (-q.w*p.y + q.x*p.z + q.y*p.w - q.z*p.x) / den,
        (-q.w*p.z - q.x*p.y + q.y*p.x + q.z*p.w) / den
    )
end

Base.:*(s::Number, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(s*p.components)

function Base.:/(s::Number, p::Q) where {Q<:AbstractQuaternion}
    f = s / (p.w^2 + p.x^2 + p.y^2 + p.z^2)
    wrapper(Q)(
        p.w * f,
        -p.x * f,
        -p.y * f,
        -p.z * f
    )
end

function Base.:*(q::Q, s::Number) where {Q<:AbstractQuaternion}
    wrapper(Q)(s*q.components)
end

function Base.:/(q::Q, s::Number) where {Q<:AbstractQuaternion}
    wrapper(Q)(q.components / s)
end

Base.:*(s::Symbolics.Num, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(s*p.components)

function Base.:/(s::Symbolics.Num, p::Q) where {Q<:AbstractQuaternion}
    f = s / (p.w^2 + p.x^2 + p.y^2 + p.z^2)
    wrapper(Q)(
        p.w * f,
        -p.x * f,
        -p.y * f,
        -p.z * f
    )
end

function Base.:*(q::Q, s::Symbolics.Num) where {Q<:AbstractQuaternion}
    wrapper(Q)(s*q.components)
end

function Base.:/(q::Q, s::Symbolics.Num) where {Q<:AbstractQuaternion}
    wrapper(Q)(q.components / s)
end

function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::AbstractQuaternion{Symbolics.Num})
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Real)
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
function Base.:(==)(q1::Real, q2::AbstractQuaternion{Symbolics.Num})
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Symbolics.Num)
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
function Base.:(==)(q1::Symbolics.Num, q2::AbstractQuaternion{Symbolics.Num})
    qdiff = Symbolics.simplify.(Symbolics.simplify(q1-q2; expand=true); expand=true)
    iszero(qdiff.w) && iszero(qdiff.x) && iszero(qdiff.y) && iszero(qdiff.z)
end
Base.:(==)(q1::AbstractQuaternion{<:Real}, q2::AbstractQuaternion{<:Real}) = (q1.w==q2.w) && (q1.x==q2.x) && (q1.y==q2.y) && (q1.z==q2.z)
Base.:(==)(q::AbstractQuaternion{<:Real}, x::Real) = isreal(q) && real(q) == x
Base.:(==)(x::Real, q::AbstractQuaternion) = isreal(q) && real(q) == x
# Base.:(==)(q::AbstractQuaternion, x::Symbolics.Num) = isreal(q) && real(q) == x
# Base.:(==)(x::Symbolics.Num, q::AbstractQuaternion) = isreal(q) && real(q) == x

Base.isequal(q1::AbstractQuaternion, q2::AbstractQuaternion) = isequal(q1.w,q2.w) && isequal(q1.x,q2.x) && isequal(q1.y,q2.y) && isequal(q1.z,q2.z)
Base.in(q::AbstractQuaternion, r::AbstractRange{<:Real}) = isreal(q) && real(q) in r

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
Base.conj(q::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(q.w, -q.x, -q.y, -q.z)
