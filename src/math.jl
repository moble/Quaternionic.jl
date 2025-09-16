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

Square root of the sum the squares of the components of the quaternion

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

Square root of the sum of the squares of the "vector" components of the quaternion.

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
    a = abs(q)
    if iszero(a)
        return Quaternion{T}(-Inf, false, false, false)
    elseif cosv ≥ 0
        v = atan(sinv, cosv)
        f = invsinc(v) / a
        return log(a) + f * quatvec(q)
    elseif iszero(sinv)  # i.e., q is a negative real number
        # Note that we check this branch only after ruling out cosv≥0 because this could
        # otherwise correspond to *positive* real numbers, which are treated correctly by
        # the preceding branch, but only the preceding branch will behave correctly for AD.
        return Quaternion{T}(log(a), false, false, π)
    else
        v′ = atan(sinv, -cosv)
        f = -invsinc(v′) * (v′-π) / v′ / a
        return log(a) + f * quatvec(q)
    end
end
function Base.log(q::Rotor{T}) where {T}
    cosv = q[1]
    sinv = absvec(q)
    if cosv ≥ 0
        v = atan(sinv, cosv)
        f = invsinc(v)
        return f * quatvec(q)
    elseif iszero(sinv)  # i.e., q is a negative real number
        # Note that we check this branch only after ruling out cosv≥0 because this could
        # otherwise correspond to *positive* real numbers, which are treated correctly by
        # the preceding branch, but only the preceding branch will behave correctly for AD.
        return QuatVec{T}(false, false, false, π)
    else
        v′ = atan(sinv, -cosv)
        f = -invsinc(v′) * (v′-π) / v′
        return f * quatvec(q)
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
    a² = abs2vec(q)
    e = exp(q[1])
    if iszero(a²)
        # Take this a little seriously, to obtain accurate AD
        ec = e*(1 - a²*(1 - a²/12)/2)
        es = e*(1 - a²*(1 - a²/20)/6)
        return Quaternion{typeof(ec)}(ec, es*q[2], es*q[3], es*q[4])
    else
        a = absvec(q)
        esinc = e * _sincu(a)
        Quaternion{typeof(esinc)}(e*cos(a), esinc*q[2], esinc*q[3], esinc*q[4])
    end
end
function Base.exp(v⃗::QuatVec{T}) where {T}
    a² = abs2vec(v⃗)
    c, s = if iszero(a²)
        # Take this a little seriously, to obtain accurate AD
        1 - a²*(1 - a²/12)/2, 1 - a²*(1 - a²/20)/6
    else
        a = absvec(v⃗)
        cos(a), _sincu(a)
    end
    Rotor{typeof(c)}(c, s*v⃗[2], s*v⃗[3], s*v⃗[4])
end

@doc raw"""
    sqrt(q)

Square root of a quaternion.  For `q` equal to a negative real number, we arbitrarily choose
the result to be proportional to `𝐤` (though any unit vector would be correct).

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

The general formula whenever the denominator is nonzero is

```math
\sqrt{q} = \frac{|q| + q} {\sqrt{2|q| + 2q[1]}}
```

This can be proven by squaring the numerator, using `q = q[1] + q⃗` and `|q|^2 = q[1]^2 -
q⃗^2`.  When the denominator is zero — or quite simply whenever `q[1] < 0` so that the
denominator may be subject to cancellation — we can use the fact that

```math
|q| + q[1] = \frac{|q⃗|^2} {|q| − q[1]}
```

to evaluate the expression in a more stable form.

Note that whenever the vector part is zero and the scalar part is negative, the solution is
not unique (and the denominator above is zero), because it necessarily involves the square
root of -1, of which there are infinitely many in the space of quaternions.  In this case,
we arbitrarily choose the vector part of the result to be proportional to `𝐤`, as mentioned
above.  A reasonable alternative would be to throw an error; instead it is left to the user
to check for that condition if it would be a problem.

Analytically, any derivative of this function will blow up as you approach the negative real
axis, because the function is discontinuous there.  Therefore, you should not expect to be
able to accurately compute derivatives of this function at points near the negative real
axis (including `q=0`) using automatic differentiation.  Specifically *on* the negative real
axis the derivative is not defined at all, but because of our fixed choice of result here,
automatic differentiation will typically (and incorrectly) return the derivative as zero.
Ultimately, the reason for this is geometrical, so it should be avoided by the user in any
case.  Therefore, for the sake of efficiency and accuracy when the problem is geometrically
well conditioned, we do not attempt to special-case the derivative at these points.

"""
function Base.sqrt(q::T) where {T<:AbstractQuaternion}
    if q[1] <= 0 && iszero(abs2vec(q))
        return T(false, false, false, √(-q[1]))
    end
    c₁ = ifelse(q[1] ≥ 0, (abs(q) + q[1]), (abs2vec(q) / (abs(q) - q[1])))
    c₂ = √inv(2c₁)
    return T(c₁*c₂, q[2]*c₂, q[3]*c₂, q[4]*c₂)
end


"""
    angle(q)

Phase angle in radians of the rotation represented by this quaternion.

Note that this may be different from your interpretation of the angle of a
complex number in an important way.  Because quaternions act on vectors by
conjugation — as in `q*v*conj(q)` — there are *two* copies of `q` involved in
that expression; in some sense, a quaternion acts "twice".  Therefore, this
angle may be twice what you expect from an analogy with complex numbers —
depending on how you interpret the correspondence between complex numbers and
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
            return Rotor{basetype(q)}([cos_πs, 0, 0, sin_πs])
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
