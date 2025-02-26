# General math functions of quaternions

"""
    abs2(q)

Sum the squares of the components of the quaternion

# Examples
```jldoctest
julia> abs2(quaternion(1,2,4,10))
121
```
"""
Base.abs2(q::AbstractQuaternion) = sum(abs2, components(q))
Base.abs2(q::AbstractQuaternion{T}) where {T<:Real} = sum(x->x^2, components(q))
Base.abs2(q::QuatVec) = sum(abs2, vec(q))
Base.abs2(q::QuatVec{T}) where {T<:Real} = sum(x->x^2, vec(q))
Base.abs2(::Rotor{T}) where {T<:Number} = one(T)

"""
    abs(q)

Square-root of the sum the squares of the components of the quaternion

This function uses Julia's built-in `hypot` function to avoid overflow and underflow.

# Examples
```jldoctest
julia> abs(quaternion(1,2,4,10))
11.0
```
"""
Base.abs(q::AbstractQuaternion) = hypot(components(q)...)
Base.abs(q::QuatVec) = hypot(vec(q)...)
Base.abs(::Rotor{T}) where {T<:Number} = one(real(T))

"""
    abs2vec(q)

Sum the squares of the "vector" components of the quaternion

# Examples
```jldoctest
julia> abs2vec(quaternion(1,2,3,6))
49
```
"""
abs2vec(q::AbstractQuaternion) = sum(abs2, vec(q))
abs2vec(q::AbstractQuaternion{T}) where {T<:Real} = sum(x->x^2, vec(q))

"""
    absvec(q)

Square-root of the sum of the squares of the "vector" components of the quaternion.

This function uses Julia's built-in `hypot` function to avoid overflow and underflow.

# Examples
```jldoctest
julia> absvec(quaternion(1,2,3,6))
7.0
```
"""
absvec(q::AbstractQuaternion) = hypot(vec(q)...)

# norm(q::Quaternion) = Base.abs2(q)  ## This might just be confusing

Base.inv(q::AbstractQuaternion) = conj(q) / abs2(q)
Base.inv(q::Rotor) = conj(q)  # Specialize to ensure output is also a Rotor


@doc raw"""
    log(q)

Logarithm of a quaternion.

As with the usual complex logarithm, the quaternion logarithm has multiple branches, though
the quaternion branches are three-dimensional: for any unit "vector" quaternion q̂, you
could add any integer multiple of 2πq̂ to the result of this function and still get the same
result after exponentiating (within numerical accuracy).  This function is the principal
logarithm.

This function has discontinuous (and fairly arbitrary) behavior along the negative real
axis: if the "vector" components of the quaternion are precisely zero *and* the scalar
component is negative, the returned quaternion will have scalar component `log(-q[1])`, but
will also have a `z` component of π.  The choice of the `z` direction is arbitrary; the
"vector" component of the returned quaternion could be π times any unit vector.

Note that `q` may be either a `Quaternion` or a `Rotor`.  If it is a `Quaternion`, this
function does not assume that it has unit norm, so the scalar component of the returned
value will generally be nonzero unless the input has *precisely* unit magnitude (which is
impossible with Float64 about 52.07% of the time due to finite machine precision), and the
return type is a `Quaternion`.  If the input is a `Rotor`, a `QuatVec` is returned, which
has scalar part exactly 0.

# Examples
```jldoctest
julia> log(exp(1.2imy))
 + 0.0𝐢 + 1.2𝐣 + 0.0𝐤

julia> log(quaternion(exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> log(quaternion(-exp(7)))
7.0 + 0.0𝐢 + 0.0𝐣 + 3.141592653589793𝐤
```

# Notes

This function uses an accurate algorithm for finding the logarithm.  For simplicity, we will
assume that the input quaternion is a unit quaternion, so that the logarithm will be a pure
vector ``\vec{v}``, which we write as ``v\hat{v}``, where ``v`` is just the scalar norm.
Note that, because of the periodicity of the `exp` function, we can assume that ``v \in [0,
\pi]`` and, in particular, ``\sin(v) \geq 0``.  Now, expand the exponential as
```math
\exp\left(\vec{v}\right) = \exp\left(v \hat{v}\right) = \cos(v) + \hat{v} \sin(v).
```
The input to this function is the right-hand side, but we do not yet know its decomposition
into ``v`` and ``\hat{v}``.  But we can find ``\cos(v)`` as the scalar part of the input,
and ``\sin(v)`` as the `absvec` (since we know that ``\sin(v) \geq 0``).  Then, we can
compute ``v = \mathrm{atan}(\sin(v), \cos(v))``.  And finally, we simply multiply the vector
part of the input by ``v / \sin(v)`` to obtain the logarithm.  This factor is given
accurately by `invsinc(v)` whenever ``|v| \leq \pi/2``.

When that condition is not satisfied (which also implies ``\cos(v)<0``, we can rewrite the
problem as
```math
\exp\left(v \hat{v}\right) = \cos(v) + \hat{v} \sin(v) = -\cos(v-\pi) - \hat{v} \sin(v-\pi)
= -\cos(v') - \hat{v} \sin(v').
```
Here, we want to multiply the vector component by ``-v / \sin(v') = -(v'+\pi) / sin(v')``.
Note that we can easily compute ``v' = \mathrm{atan}(\sin(v), -\cos(v))``.  This algorithm
is surprisingly accurate, even when ``v`` is extremely close to ``\pi``, which implies that
the vector part of the input is extremely small.

The only special case remaining to handle is when ``\cos(v) < 0`` but ``\sin(v)`` is
*identically* zero.  In this case, we could throw an error, but this is not usually helpful.
Instead, we arbitrarily choose to return ``\pi z``.

Finally, we can return to our assumption that the input has unit magnitude.  In the
preceding, this didn't matter because we only computed ``v`` as the *ratio* of the scalar
part and the magnitude of the vector part, so the overall magnitude cancelled out.  So that
part of the computation remains unchanged.  Instead, we note that for scalar ``s``, we have
```math
\exp\left(s + \vec{v}\right) = \exp\left(s\right) \exp\left(\vec{v}\right),
```
so the logarithm is just the sum of the logarithms of the scalar and vector parts.

"""
function Base.log(q::Quaternion{T}) where {T}
    cosv = q[1]
    sinv = absvec(q)
    a² = abs2(q)
    if iszero(a²)
        return Quaternion{T}(-Inf, false, false, false)
    elseif cosv ≥ 0
        v = atan(sinv, cosv)
        f = invsinc(v) / √a²
        return quaternion(log(a²)/2, f * vec(q)...)
    elseif iszero(sinv)
        return Quaternion{T}(log(a²)/2, false, false, π)
    else
        v′ = atan(sinv, -cosv)
        f = -invsinc(v′) * (v′-π) / v′ / √a²
        return quaternion(log(a²)/2, f * vec(q)...)
    end
end
function Base.log(q::Rotor{T}) where {T}
    cosv = q[1]
    sinv = absvec(q)
    if cosv ≥ 0
        v = atan(sinv, cosv)
        f = invsinc(v)
        return quatvec(f * vec(q))
    elseif iszero(sinv)
        return QuatVec{T}(false, false, false, π)
    else
        v′ = atan(sinv, -cosv)
        f = -invsinc(v′) * (v′-π) / v′
        return quatvec(f * vec(q))
    end
end

"""
    exp(q)

Exponential of a quaternion

# Examples
```jldoctest
julia> exp(imx*π/4)  # Rotation through π/2 (note the extra 1/2) about the x axis
rotor(0.7071067811865476 + 0.7071067811865475𝐢 + 0.0𝐣 + 0.0𝐤)
```
"""
function Base.exp(q::Quaternion{T}) where {T}
    a = abs2vec(q)
    e = exp(q[1])
    # let ecos, esinc = ifelse(iszero(a), e*(1 - a*(1 - a/12)/2, e), (e*cos(√a), e*_sincu(√a)))
    #     Rotor{typeof(ecos)}(ecos, esinc*v⃗[2], esinc*v⃗[3], esinc*v⃗[4])
    # end
    if iszero(a)
        # Take this a little seriously, to obtain accurate ForwardDiff derivative
        Quaternion{typeof(e)}(e*(1 - a*(1 - a/12)/2), e*q[2], e*q[3], e*q[4])
    else
        esinc = e * _sincu(√a)
        Quaternion{typeof(esinc)}(e*cos(√a), esinc*q[2], esinc*q[3], esinc*q[4])
    end
end
function Base.exp(v⃗::QuatVec{T}) where {T}
    a = abs2vec(v⃗)
    c, s = if iszero(a)
        # Take this a little seriously, to obtain accurate ForwardDiff derivative
        1 - a*(1 - a/12)/2, one(a)
    else
        cos(√a), _sincu(√a)
    end
    Rotor{typeof(c)}(c, s*v⃗[2], s*v⃗[3], s*v⃗[4])
end

@doc raw"""
    sqrt(q)

Square-root of a quaternion.

The general formula whenever the denominator is nonzero is

```math
\sqrt{q} = \frac{|q| + q} {\sqrt{2|q| + 2q[1]}}
```

This can be proven by expanding `q` as `q[1] + vec(q)` and multiplying the expression above
by itself.

Note that whenever the vector part is zero and the scalar part is negative, the solution is
not unique (and the denominator above is zero), because it necessarily involves the
square-root of -1, of which there are infinitely many in the space of quaternions.  In this
case, we arbitrarily choose the vector part of the result to be in the `z` direction.  A
reasonable alternative would be to throw an error; instead it is left to the user to check
for that condition.

# Examples
```jldoctest
julia> q = quaternion(1.2, 3.4, 5.6, 7.8);

julia> sqrtq = √q;

julia> sqrtq^2 ≈ q
true

julia> √quaternion(4.0)
2.0 + 0.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> √quaternion(-4.0)
0.0 + 0.0𝐢 + 0.0𝐣 + 2.0𝐤
```

# Notes

This function uses an algorithm for finding the square root that is very accurate (typically
achieving the correct result to within machine precision) for *most* values.  However, naive
application of the formula above can lead to catastrophic cancellation when the scalar part
is negative and significantly larger in magnitude to the vector part.  Therefore, when `q[1]
 < 0`, we transform the problem into the case where `q[1] > 0` as
```math
\sqrt{q} = \bar{\sqrt{-\bar{q}}}\, \sqrt{-1},
```
where the bar denotes quaternionic conjugation, and we interpret ``\sqrt{-1}`` to be a unit
imaginary that commutes with ``q``.  The obvious candidate is the normalized vector part of
``q`` if the vector part is nonzero; otherwise, as noted above, we arbitrarily choose it to
be the unit vector in `z` direction.  The calculation of ``\sqrt{-\bar{q}}`` can use the
naive approach, and the result is still accurate to within machine precision.

"""
function Base.sqrt(q::T) where {T<:AbstractQuaternion}
    if q[1] ≥ 0
        if all(iszero, vec(q))
            T(√(q[1]), false, false, false)
        else
            c₁ = abs(q) + q[1]
            c₂ = √inv(2c₁)
            T(c₁*c₂, q[2]*c₂, q[3]*c₂, q[4]*c₂)
        end
    else # q[1] < 0
        if all(iszero, vec(q))
            T(false, false, false, √(-q[1]))
        else
            c₁ = abs(q) - q[1]
            c₂ = √inv(2c₁)
            T(c₁*c₂, -q[2]*c₂, -q[3]*c₂, -q[4]*c₂) * T(normalize(vec(q))...)
        end
    end
end


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
Base.angle(q::Quaternion{T}) where T = 2 * absvec(log(q))
Base.angle(q::Rotor{T}) where T = 2 * absvec(log(q))


function Base.:^(q::Quaternion, s::Number)
    exp(s * log(q))
end
function Base.:^(q::Rotor, s::Number)
    q = float(q)
    absolutevec = absvec(q)
    if absolutevec ≤ eps(typeof(absolutevec))
        if q[1] < 0
            # log(q) ≈ π𝐤
            sin_πs, cos_πs = sincospi(oftype(absolutevec, s))
            return Rotor{eltype(q)}([cos_πs, 0, 0, sin_πs])
        end
        # log(q) ≈ 0
        return one(q)
    end
    f1 = s * atan(absolutevec, q[1])
    sin_f1, cos_f1 = sincos(f1)
    f2 = sin_f1 / absolutevec
    Rotor{typeof(cos_f1)}([cos_f1, f2*q[2], f2*q[3], f2*q[4]])
end

# We need to be more specific about the quaternion types here because we had
# to be specific about the quaternion type above, and Integer<:Number
for QT ∈ [AbstractQuaternion, Quaternion, QuatVec, Rotor]
    @eval function Base.:^(q::$QT, s::Integer)
        if s ≥ 0
            Base.power_by_squaring(q, s)
        else
            inv(Base.power_by_squaring(q, -s))
        end
    end
end
