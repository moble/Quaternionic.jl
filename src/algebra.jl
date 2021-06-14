Base.:-(q::Quaternion) = Quaternion(-q.components)

Base.:+(q::Quaternion, p::Quaternion) = Quaternion(q.components+p.components)
Base.:-(q::Quaternion, p::Quaternion) = Quaternion(q.components-p.components)

Base.:+(q::Quaternion, p::Number) = Quaternion(q.w+p, q.x, q.y, q.z)
Base.:-(q::Quaternion, p::Number) = Quaternion(q.w-p, q.x, q.y, q.z)

Base.:+(q::Number, p::Quaternion) = Quaternion(q+p.w, p.x, p.y, p.z)
Base.:-(q::Number, p::Quaternion) = Quaternion(q-p.w, -p.x, -p.y, -p.z)

Base.:+(q::Quaternion, p::Symbolics.Num) = Quaternion(q.w+p, q.x, q.y, q.z)
Base.:-(q::Quaternion, p::Symbolics.Num) = Quaternion(q.w-p, q.x, q.y, q.z)

Base.:+(q::Symbolics.Num, p::Quaternion) = Quaternion(q+p.w, p.x, p.y, p.z)
Base.:-(q::Symbolics.Num, p::Quaternion) = Quaternion(q-p.w, -p.x, -p.y, -p.z)

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
