
"""
    abs2(q)

Sum the squares of the components of the quaternion

# Examples
```jldoctest
julia> abs2(Quaternion(1,2,4,10))
121
```
"""
Base.abs2(q::Quaternion) = sum(q.components.^2)

"""
    abs(q)

Square-root of the sum the squares of the components of the quaternion

# Examples
```jldoctest
julia> abs(Quaternion(1,2,4,10))
11.0
```
"""
Base.abs(q::Quaternion) = sqrt(abs2(q))

"""
    abs2vec(q)

Sum the squares of the "vector" components of the quaternion

# Examples
```jldoctest
julia> abs2vec(Quaternion(1,2,3,6))
49
```
"""
abs2vec(q::Quaternion) = @inbounds q.components[2]^2 + q.components[3]^2 + q.components[4]^2

"""
    absvec(q)

Square-root of the sum of the squares of the "vector" components of the quaternion

# Examples
```jldoctest
julia> absvec(Quaternion(1,2,3,6))
7.0
```
"""
absvec(q::Quaternion) = sqrt(abs2vec(q))
# norm(q::Quaternion) = Base.abs2(q)  ## This might just be confusing
Base.inv(q::Quaternion) = conj(q) / abs2(q)

"""
    log(q)

Logarithm of a quaternion.

As with the usual complex logarithm, the quaternion logarithm has multiple
branches, though the quaternion branches are three-dimensional: for any unit
"vector" quaternion q̂, you could add any integer multiple of 2πq̂ to the result
of this function and still get the same result after exponentiating (within
numerical accuracy).  This function is the principal logarithm.

This function has discontinuous (and fairly arbitrary) behavior along the
negative real axis: if the "vector" components of the quaternion are precisely
zero *and* the scalar component is negative, the returned quaternion will have
scalar component `log(-q.w)`, but will also have a `z` component of π.  The
choice of the `z` direction is arbitrary; the "vector" component of the
returned quaternion could be π times any unit vector.

Note that this function is not specialized to unit-quaternion inputs, so the
scalar component of the returned value will be nonzero unless the input has
*precisely* unit magnitude.

# Examples
```jldoctest
julia> log(exp(1.2imy))
0.0 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤

julia> log(Quaternion(exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> log(Quaternion(-exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 3.141592653589793𝐤
```

"""
function Base.log(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        if q.w < 0
            return Quaternion(log(-q.w), 0, 0, π)
        end
        return Quaternion(log(q.w), 0, 0, 0)
    end
    absolutevec = sqrt(absolute2vec)
    f = atan(absolutevec, q.w) / absolutevec  # acos((w^2-absolutevec^2) / (w^2+absolutevec^2)) / 2absolutevec
    Quaternion(log(abs2(q))/2, f*q.x, f*q.y, f*q.z)
end
# function Base.log(q::UnitQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if iszero(absolute2vec)
#         if q.w < 0
#             return Quaternion{T}(0, 0, 0, π)
#         end
#         return Quaternion{T}(0, 0, 0, 0)
#     end
#     absolutevec = sqrt(absolute2vec)
#     f = atan(absolutevec, q.w) / absolutevec  # acos((w^2-absolutevec^2) / (w^2+absolutevec^2)) / 2absolutevec
#     Quaternion(0, f*q.x, f*q.y, f*q.z)
# end

"""
    exp(q)

Exponential of a quaternion

# Examples
```jldoctest
julia> exp(imx*π/4)  # Rotation through π/2 (note the extra 1/2) about the x axis
0.7071067811865476 + 0.7071067811865475𝐢 + 0.0𝐣 + 0.0𝐤
```
"""
function Base.exp(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        return Quaternion(exp(q.w), 0, 0, 0)
    end
    absolutevec = sqrt(absolute2vec)
    s = sin(absolutevec) / absolutevec
    e = exp(q.w)
    Quaternion(e*cos(absolutevec), e*s*q.x, e*s*q.y, e*s*q.z)
end
# function Base.exp(q::VectorQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if iszero(absolute2vec)
#         return Quaternion{T}(1, 0, 0, 0)
#     end
#     absolutevec = sqrt(absolute2vec)
#     s = sin(absolutevec) / absolutevec
#     Quaternion(cos(absolutevec), s*q.x, s*q.y, s*q.z)
# end

@doc raw"""
    sqrt(q)

Square-root of a quaternion.

The general formula whenever the denominator is nonzero is

``
\sqrt{q} = \frac{|q| + q} {\sqrt{2|q| + 2q.w}}
``

This can be proven by expanding `q` as `q.w + q.vec` and multiplying the
expression above by itself.

When the denominator is zero, this function has discontinuous (and fairly
arbitrary) behavior, just as with the quaternion [`log`](@ref) function.  In
this case, either all components are zero — in which case the result is simply
the zero quaternion — or the "vector" components of the quaternion are
precisely zero and the scalar component is negative.  If the latter is true,
the denominator above will be a pure-imaginary number.  Because the quaternions
come with infinitely many elements that square to -1, it is not clear *which*
imaginary should be used, so we arbitrarily choose to set the result
proportional to the `z` quaternion.  The choice of the `z` direction is
arbitrary; the "vector" component of the returned quaternion could be in any
direction.

# Examples
```jldoctest
julia> q = Quaternion(1.2, 3.4, 5.6, 7.8);

julia> sqrtq = √q;

julia> sqrtq^2 ≈ q
true

julia> √Quaternion(4)
2.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> √Quaternion(-4)
0.0 + 0.0𝐢 + 0.0𝐣 + 2.0𝐤
```
"""
function Base.sqrt(q::Quaternion{T}) where {T}
    q = float(q)
    absolute2vec = abs2vec(q)
    if iszero(absolute2vec)
        if q.w < 0
            return Quaternion(0, 0, 0, sqrt(-q.w))
        end
        return Quaternion(sqrt(q.w), 0, 0, 0)
    end
    absolute2 = absolute2vec + q.w^2
    c1 = sqrt(absolute2) + q.w
    c2 = sqrt(inv(2*c1))
    Quaternion(c1*c2, q.x*c2, q.y*c2, q.z*c2)
end
# function Base.sqrt(q::UnitQuaternion{T}) where {T}
#     q = float(q)
#     absolute2vec = abs2vec(q)
#     if absolute2vec == zero(T)
#         if q.w < 0
#             return Quaternion{T}(0, 0, 0, 1)
#         end
#         return Quaternion{T}(1, 0, 0, 0)
#     end
#     c1 = 1 + q.w
#     c2 = sqrt(inv(2*c1))
#     Quaternion(c1*c2, q.x*c2, q.y*c2, q.z*c2)
# end

"""
    angle(q)

Phase angle in radians of the rotation represented by this quaternion.

Note that this may be different from your interpretation of the angle of a
complex number in an important way.  Because quaternions act on vectors by
conjugation — as in `q*v*conj(q)` — there are *two* copies of `q` involved in
that expression; in some sense, a quaternion acts "twice".  Therefore, this
angle may be twice what you expect from an analogy with complex numbers —
dpending on how you interpret the correspondence between complex numbers and
quaternions.  Also, while rotations in the complex plane have a natural choice
of axis (the positive `z` direction), that is not the case for quaternions,
which means that the sign of this angle is arbitrary, and we always choose it
to be positive.

# Examples
```jldoctest
julia> θ=1.2;

julia> R=exp(θ * imz / 2);

julia> angle(R)
1.2

```
"""
Base.angle(q::Quaternion) = 2 * absvec(log(q))
# Base.angle(q::UnitQuaternion) = 2 * absvec(log(q))
