# Essential elements of making quaternions into an a algebra

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

Base.:+(q::Rotor, p::Rotor) = Quaternion(q.components+p.components)
Base.:-(q::Rotor, p::Rotor) = Quaternion(q.components-p.components)
Base.:+(q::Rotor, p::Q) where {Q<:AbstractQuaternion} = Quaternion(q.components+p.components)
Base.:-(q::Rotor, p::Q) where {Q<:AbstractQuaternion} = Quaternion(q.components-p.components)
Base.:+(q::Q, p::Rotor) where {Q<:AbstractQuaternion} = Quaternion(q.components+p.components)
Base.:-(q::Q, p::Rotor) where {Q<:AbstractQuaternion} = Quaternion(q.components-p.components)

Base.:+(q::Rotor, p::Number) = Quaternion(q.w+p, q.x, q.y, q.z)
Base.:-(q::Rotor, p::Number) = Quaternion(q.w-p, q.x, q.y, q.z)

Base.:+(q::Number, p::Rotor) = Quaternion(q+p.w, p.x, p.y, p.z)
Base.:-(q::Number, p::Rotor) = Quaternion(q-p.w, -p.x, -p.y, -p.z)

Base.:+(q::Rotor, p::Symbolics.Num) = Quaternion(q.w+p, q.x, q.y, q.z)
Base.:-(q::Rotor, p::Symbolics.Num) = Quaternion(q.w-p, q.x, q.y, q.z)

Base.:+(q::Symbolics.Num, p::Rotor) = Quaternion(q+p.w, p.x, p.y, p.z)
Base.:-(q::Symbolics.Num, p::Rotor) = Quaternion(q-p.w, -p.x, -p.y, -p.z)

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
