# Essential elements of making quaternions into an algebra

"""
    conj(q)

Return the quaternion conjugate, which flips the sign of each "vector"
component.

# Examples
```jldoctest
julia> conj(quaternion(1,2,3,4))
1 - 2𝐢 - 3𝐣 - 4𝐤
```
"""
Base.conj(q::Q) where {Q<:AbstractQuaternion} = Q(q[1], -q[2], -q[3], -q[4])
Base.conj(q::Q) where {Q<:AbstractQuaternion{Bool}} = wrapper(Q)(q[1], -q[2], -q[3], -q[4])

Base.:-(q::Q) where {Q<:AbstractQuaternion} = Q(-components(q))
Base.:-(q::Q) where {Q<:AbstractQuaternion{Bool}} = wrapper(Q)(-Int.(components(q)))

# Note that, in the two definitions above, the `wrapper` function is used to
# ensure that the result is of the same type as the input quaternion, but with a different
# basetype.  This is necessary because `-true` is not a `Bool`, it's the `Int` -1.  This is
# the only type I can think of where `-` changes the type of the input, so that's the only
# special case we need.  (I'm assuming `Unsigned` types are not used here.)



for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    for TB ∈ (AbstractQuaternion, Rotor, QuatVec)
        @eval begin
            Base.:+(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(+), $TB)(components(q)+components(p))
            Base.:-(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(-), $TB)(components(q)-components(p))
        end
    end
    let TB = Number
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end


function Base.:*(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    wrapper(Q1, Val(*), Q2)(
        q[1]*p[1] - q[2]*p[2] - q[3]*p[3] - q[4]*p[4],
        q[1]*p[2] + q[2]*p[1] + q[3]*p[4] - q[4]*p[3],
        q[1]*p[3] - q[2]*p[4] + q[3]*p[1] + q[4]*p[2],
        q[1]*p[4] + q[2]*p[3] - q[3]*p[2] + q[4]*p[1]
    )
end


function Base.:/(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    if p == q
        return one(promote_type(Q1, float(Q2)))
    end
    den = abs2(p)
    wrapper(Q1, Val(/), Q2)(
        (+q[1]*p[1] + q[2]*p[2] + q[3]*p[3] + q[4]*p[4]) / den,
        (-q[1]*p[2] + q[2]*p[1] - q[3]*p[4] + q[4]*p[3]) / den,
        (-q[1]*p[3] + q[2]*p[4] + q[3]*p[1] - q[4]*p[2]) / den,
        (-q[1]*p[4] - q[2]*p[3] + q[3]*p[2] + q[4]*p[1]) / den
    )
end


let S = Number
    @eval begin
        Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(*), $S)(s*components(p))
        Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion} = wrapper($S, Val(*), Q)(s*components(p))
        Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(/), $S)(components(p)/s)
        function Base.:/(s::$S, p::Q) where {Q<:AbstractQuaternion}
            f = s / abs2(p)
            wrapper($S, Val(/), Q)(p[1] * f, -p[2] * f, -p[3] * f, -p[4] * f)
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
@inline function LinearAlgebra.:⋅(p::AbstractQuaternion, q::AbstractQuaternion)
    p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4]
end

"""
    a × b

Return the cross product of two pure-vector quaternions.  Equal to ½ of the
commutator product `a*b-b*a`.
"""
@inline function ×(a::QuatVec, b::QuatVec)
    quatvec(
        a[3] * b[4] - a[4] * b[3],
        a[4] * b[2] - a[2] * b[4],
        a[2] * b[3] - a[3] * b[2]
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
convert to a `Rotor`, just call `rotor(q)`, which includes a normalization
step.
"""
@inline function normalize(q::AbstractQuaternion)
    return q / abs(q)
end
@inline function normalize(q::Rotor)
    return rotor(q)  # already normalizes
end

function (R::Rotor)(v::QuatVec)
    quatvec(SA[
        false,
        ((R[1]^2 + R[2]^2 - R[3]^2 - R[4]^2)*v[2]
            + (R[1]*R[3] + R[2]*R[4])*2v[4] + (R[2]*R[3] - R[1]*R[4])*2v[3]),
        ((R[1]^2 - R[2]^2 + R[3]^2 - R[4]^2)*v[3]
            + (R[2]*R[3] + R[1]*R[4])*2v[2] + (R[3]*R[4] - R[1]*R[2])*2v[4]),
        ((R[1]^2 + R[4]^2 - R[2]^2 - R[3]^2)*v[4]
            + (R[1]*R[2] + R[3]*R[4])*2v[3] + (R[2]*R[4] - R[1]*R[3])*2v[2])
    ])
end
function (R::Quaternion)(v::QT) where {QT<:AbstractQuaternion}
    wrapper(QT)(
        abs2(R) * v[1],
        ((R[1]^2 + R[2]^2 - R[3]^2 - R[4]^2)*v[2]
            + (R[1]*R[3] + R[2]*R[4])*2v[4] + (R[2]*R[3] - R[1]*R[4])*2v[3]),
        ((R[1]^2 - R[2]^2 + R[3]^2 - R[4]^2)*v[3]
            + (R[2]*R[3] + R[1]*R[4])*2v[2] + (R[3]*R[4] - R[1]*R[2])*2v[4]),
        ((R[1]^2 + R[4]^2 - R[2]^2 - R[3]^2)*v[4]
            + (R[1]*R[2] + R[3]*R[4])*2v[3] + (R[2]*R[4] - R[1]*R[3])*2v[2])
    )
end
