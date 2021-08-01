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
    for TB ∈ [Real, Symbolics.Num]
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q.w+p, q.x, q.y, q.z)
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q.w-p, q.x, q.y, q.z)

            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q.w, q.x, q.y, q.z)
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q.w, -q.x, -q.y, -q.z)
        end
    end
end


function Base.:*(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    wrapper(Q1, Val(*), Q2)(
        q.w*p.w - q.x*p.x - q.y*p.y - q.z*p.z,
        q.w*p.x + q.x*p.w + q.y*p.z - q.z*p.y,
        q.w*p.y - q.x*p.z + q.y*p.w + q.z*p.x,
        q.w*p.z + q.x*p.y - q.y*p.x + q.z*p.w
    )
end


function Base.:/(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    if p == q
        return one(promote_type(Q1, float(Q2)))
    end
    den = abs2(p)
    wrapper(Q1, Val(/), Q2)(
        (+q.w*p.w + q.x*p.x + q.y*p.y + q.z*p.z) / den,
        (-q.w*p.x + q.x*p.w - q.y*p.z + q.z*p.y) / den,
        (-q.w*p.y + q.x*p.z + q.y*p.w - q.z*p.x) / den,
        (-q.w*p.z - q.x*p.y + q.y*p.x + q.z*p.w) / den
    )
end


for S ∈ [Real, Symbolics.Num]
    @eval begin
        Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(*), $S)(s*p.components)
        Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion} = wrapper($S, Val(*), Q)(s*p.components)
        Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(/), $S)(p.components/s)
        function Base.:/(s::$S, p::Q) where {Q<:AbstractQuaternion}
            f = s / abs2(p)
            wrapper($S, Val(/), Q)(p.w * f, -p.x * f, -p.y * f, -p.z * f)
        end
    end
end


"""
    p ⋅ q

Evaluate the inner ("dot") product between two quaternions.  Equal to the
scalar part of `p * conj(q)`.

Note that this function is not very commonly used, except as a quick way to
determine whether the two quaternions are more anti-parallel than parallel, for
functions like [`unflip`](@ref).
"""
@inline function ⋅(p::AbstractQuaternion, q::AbstractQuaternion)
    p.w*q.w + p.x*q.x + p.y*q.y + p.z*q.z
end

"""
    a × b

Return the cross product of two pure-vector quaternions.  Equal to ½ of the
commutator product `a*b-b*a`.
"""
@inline function ×(a::QuatVec, b::QuatVec)
    QuatVec(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    )
end

"""
    a ×̂ b

Return the *direction* of the cross product between `a` and `b`; the normalized
vector along `a×b` — unless the magnitude is zero, in which case the zero
vector is returned.
"""
@inline function ×̂(a::QuatVec, b::QuatVec)
    axb = a × b
    av = absvec(axb)
    if av == zero(av)
        return axb
    else
        return axb / av
    end
end

"""
    normalize(q)

Return a copy of this quaternion, normalized.

Note that this returns the same type as the input quaternion.  If you want to
convert to a `Rotor`, just call `Rotor(q)`, which includes a normalization
step.
"""
@inline function normalize(q::AbstractQuaternion)
    return q / abs(q)
end
@inline function normalize(q::Rotor)
    return Rotor(q)  # already normalizes
end
