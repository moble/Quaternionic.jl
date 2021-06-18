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


for TA ∈ [AbstractQuaternion, Rotor, QuatVec]
    for TB ∈ [AbstractQuaternion, Rotor, QuatVec]
        @eval begin
            Base.:+(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(+), $TB)(q.components+p.components)
            Base.:-(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(-), $TB)(q.components-p.components)
        end
    end
    for TB ∈ [Number, Symbolics.Num]
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q.w+p, q.x, q.y, q.z)
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q.w-p, q.x, q.y, q.z)

            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q.w, q.x, q.y, q.z)
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p-q.w, -q.x, -q.y, -q.z)
        end
    end
end


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
    den = abs2(p)
    wrapper(Q1, Q2)(
        (+q.w*p.w + q.x*p.x + q.y*p.y + q.z*p.z) / den,
        (-q.w*p.x + q.x*p.w - q.y*p.z + q.z*p.y) / den,
        (-q.w*p.y + q.x*p.z + q.y*p.w - q.z*p.x) / den,
        (-q.w*p.z - q.x*p.y + q.y*p.x + q.z*p.w) / den
    )
end

Base.:*(s::Number, p::Q) where {Q<:AbstractQuaternion} = wrapper(Q)(s*p.components)
Base.:*(s::Number, p::Rotor) = wrapper(Q)(signbit(s) ? -p.components : p.components)

function Base.:/(s::Number, p::Q) where {Q<:AbstractQuaternion}
    f = s / abs2(p)
    wrapper(Q)(
        p.w * f,
        -p.x * f,
        -p.y * f,
        -p.z * f
    )
end

Base.:*(q::Q, s::Number) where {Q<:AbstractQuaternion} = wrapper(Q)(s*q.components)

Base.:/(q::Q, s::Number) where {Q<:AbstractQuaternion} = wrapper(Q)(q.components / s)

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

Base.:*(q::Q, s::Symbolics.Num) where {Q<:AbstractQuaternion} = wrapper(Q)(s*q.components)

Base.:/(q::Q, s::Symbolics.Num) where {Q<:AbstractQuaternion} = wrapper(Q)(q.components / s)
