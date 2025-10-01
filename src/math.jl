# General math functions of quaternions

"""
    abs2(q)

Sum of the squares of all four components of the quaternion.

Note that the result for a `Rotor` is identically 1, even if that is not true numerically.

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
Base.abs2(::Rotor{T}) where {T<:Number} = one(real(T))

"""
    abs(q)

Square root of the sum of the squares of all four components of the quaternion.

This function uses Julia's built-in `hypot` function to avoid overflow and underflow.

Note that the result for a `Rotor` is identically 1, even if that is not true numerically.

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

Sum of the squares of the three "vector" components of the quaternion.

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

Square root of the sum of the squares of the three "vector" components of the quaternion.

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

!!! note "Branch-cut behavior"

    As with the complex logarithm, the quaternion logarithm is multi-valued: you
    could add any integer multiple of ``2πq̂`` (for some unit vector ``q̂``) to the result of
    this function and get the same result after exponentiating.  This function is the
    principal logarithm: the choice with the smallest-magnitude vector part.

    Similarly, a branch cut is imposed along the negative real axis: if the vector
    components of `q` are precisely zero *and* the scalar component is negative, the
    returned quaternion will be `log(-q[1]) + π𝐤`.  Unlike in the complex case, the choice
    of `𝐤` is arbitrary; the vector component of the returned quaternion could be π
    times any unit vector.

!!! warning "Automatic differentiation caveats"

    Values of `q` that are very close to (but not on) the negative real axis will produce
    accurate results, but automatic differentiation is likely to become numerically
    unstable.  *Precisely on* the non-positive real axis, derivatives should not be defined
    at all, but will typically return incorrect (finite) values.  These regions should be
    avoided in any case, because of the analytic discontinuities and geometric ambiguities.

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

This function uses an accurate algorithm for finding the logarithm, although automatic
differentiation may become inaccurate near the negative real axis because the analytic
derivatives blow up there.

The `log` function is very analogous to the [`sqrt`](@ref) function [which is essentially
`exp(log(q)/2)`], in that both functions have discontinuous behavior along the negative real
axis, and derivatives that blow up as you approach that axis.  Therefore, the same caveats
about computing derivatives using automatic differentiation apply here as well.
Specifically, those derivatives will be numerically unstable near the negative real axis,
and will be incorrect *on* the negative real axis itself; rather than being undefined, the
derivatives will typically be returned as zero because of our fixed choice of result there.

If we decompose the result of this function as ``\log(q) = s + \vec{v}``, where ``s`` is the
scalar part and ``\vec{v}`` is the pure-vector part, then clearly ``s`` and ``\vec{v}``
commute, so their exponentials also commute, and we have
```math
q = \exp\left(\log(q)\right) = \exp\left(s\right) \exp\left(\vec{v}\right).
```
Note that the exponential of a pure-vector quaternion is a unit quaternion, so we have
decomposed ``q`` into a product of a positive real number and a unit quaternion.  We already
know how to take the logarithm of a positive real number.  Therefore, in the following we
assume that ``q`` is a unit quaternion, so that ``s=0`` and we only need to find
``\vec{v}``.

We now write ``\log(q)`` as ``v\hat{v}``, where ``v`` is just the scalar norm.  Note that,
because of the periodicity of the `exp` function, we can assume that ``v \in [0, \pi]`` and,
in particular, ``\sin(v) \geq 0``.  Now, expand the exponential as
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
Here, we want to multiply the vector component by ``-v / \sin(v') = -(v'+\pi) / \sin(v')``.
Note that we can easily compute ``v' = \mathrm{atan}(\sin(v), -\cos(v))``.  This algorithm
is surprisingly accurate, even when ``v`` is extremely close to ``\pi``, which implies that
the vector part of the input is extremely small.

The only special case remaining to handle is when ``\cos(v) < 0`` but ``\sin(v)`` is
*identically* zero.  In this case, we could throw an error, but this is not usually helpful.
Instead, we arbitrarily choose to return ``\pi 𝐤``.

If `q` is a `Rotor`, we return a `QuatVec`; if `q` is a general `Quaternion`, we return a
general `Quaternion` — though if `q` happens to have norm exactly 1, the result will have a
scalar part of exactly 0.  Note that, because there is no geometric reason to take the
logarithm of a `QuatVec`, that case is not implemented; if you really need to compute it,
you can convert the `QuatVec` to a `Quaternion` first.

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

@doc raw"""
    exp(q)

Exponential of a quaternion.

The exponential of a quaternion is defined as usual by its power series, which converges for
all finite quaternions:
```math
\exp(q) = \sum_{k=0}^\infty \frac{q^k}{k!}.
```

!!! note "Derivatives at 0"

    Note that automatic differentiation of this function *at exactly* `q=0` will only be
    accurate through fifth order.  Using sixth-order derivatives or higher is so unusual
    that this will likely not be a problem in practice.

# Examples
```jldoctest
julia> R = exp(imx*π/4)  # Rotation through π/2 (note the extra 1/2) about the x axis
rotor(0.7071067811865476 + 0.7071067811865475𝐢 + 0.0𝐣 + 0.0𝐤)

julia> R * imx * conj(R)
0.0 + 1.0𝐢 + 0.0𝐣 + 0.0𝐤

julia> R * imy * conj(R)
0.0 + 0.0𝐢 + (2.220446049250313e-16)𝐣 + 1.0𝐤

julia> R * imz * conj(R)
0.0 + 0.0𝐢 - 1.0𝐣 + (2.220446049250313e-16)𝐤
```

# Notes

The quaternionic exponential is very easy to calculate by analogy with the complex
exponential.  We can write a quaternion as ``q = s + v\hat{v}``, where ``s`` is the scalar
part, ``v`` is the norm of the pure-vector part, and ``\hat{v}`` is a unit vector.  Then,
``\hat{v}^2 = -1``, so it acts exactly like the imaginary unit ``i`` in complex numbers.
Obviously, ``s``, ``v``, and ``\hat{v}`` all commute with each other, so the math is simple
and we can immediately calculate in analogy with Euler's formula for complex numbers:
```math
\exp(q) = \exp(s) \left(\cos(v) + \hat{v}\sin(v)\right).
```

The only special case is when ``v=0``, in which case there is no unique choice of
``\hat{v}``, but then `exp(q) = exp(s)`.  Unfortunately, this means that we need to use a
separate branch for this isolated case, which means that automatic differentiation with
respect to the vector components will incorrectly produce a derivative of zero at that point
if we do nothing but return the value.  Therefore, we actually use a Taylor expansion to
compute the result through fifth order in the magnitude of the vector part.  This will
produce correct derivatives up to fifth order.
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

Square root of a quaternion.

!!! note "Branch-cut behavior"

    As with the logarithm, the quaternionic square-root has a branch cut along the
    non-positive real axis: if the vector components of `q` are precisely zero *and* the
    scalar component is negative, the returned quaternion will be `√(-q[1]) * 𝐤`.  Unlike in
    the complex case, the choice of `𝐤` is arbitrary; the vector component of the returned
    quaternion could be any unit vector.

!!! warning "Automatic differentiation caveats"

    Values of `q` that are very close to (but not on) the negative real axis will produce
    accurate results, but automatic differentiation is likely to become numerically
    unstable.  *Precisely on* the non-positive real axis, derivatives should not be defined
    at all, but will typically return incorrect (finite) values.  These regions should be
    avoided in any case, because of the analytic discontinuities and geometric ambiguities.


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
\sqrt{q} = \frac{|q| + q} {\sqrt{2(|q| + q[1])}}
```

This can be proven by squaring the numerator, using `q = q[1] + q⃗` and `|q|^2 = q[1]^2 -
q⃗^2 = q[1]^2 + |q⃗|^2`.  When the denominator is zero — or quite simply whenever `q[1] < 0`
so that the denominator may be subject to cancellation — we can use the fact that

```math
|q| + q[1] = \frac{|q⃗|^2} {|q| − q[1]}
```

to evaluate the expression in a more stable form:

```math
\sqrt{q} = \left(|q| + q\right)\sqrt{\frac{|q| − q[1]} {2|q⃗|^2}}.
```

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
    ## Work around https://github.com/chalk-lab/Mooncake.jl/issues/794
    # c₁ = ifelse(q[1] ≥ 0, (abs(q) + q[1]), (abs2vec(q) / (abs(q) - q[1])))
    c₁ = if q[1] ≥ 0
        (abs(q) + q[1])
    else
        (abs2vec(q) / (abs(q) - q[1]))
    end
    c₂ = √inv(2c₁)
    return T(c₁*c₂, q[2]*c₂, q[3]*c₂, q[4]*c₂)
end


"""
    angle(q)

Phase angle in radians of the (spinorial) rotation represented by this quaternion.

Note that this may be different from your interpretation of the angle of a complex number in
an important way.  Because quaternions act on vectors by conjugation — as in `q*v*conj(q)` —
there are *two* copies of `q` involved in that expression; in some sense, a quaternion acts
"twice".  Therefore, this angle may be twice what you expect from an analogy with complex
numbers — depending on how you interpret the correspondence between complex numbers and
quaternions.

Also, because quaternions are spinors, rotation by 2π is not the identity operation; only a
full rotation through 4π is.  Rotation by 2π changes the sign of the quaternion, which has
no effect on the rotation of *vectors*, but does have an effect on rotations of more general
objects (spinors).  This issue of sign is important when considering the continuity of
quaternionic functions (as in interpolations, differentiation, and integration, for
example).  Therefore, the angle returned by this function will be in the range `[0, 2π]`,
rather than the range `[0, π]` used for complex numbers.  If you only care about the effects
on vectors, and want to map the result of this function to the range `[0, π]`, you can use
`θ -> min(θ, 2π-θ)`.

# Examples
```jldoctest
julia> θ=1.2;

julia> R=exp(θ * imz / 2);  # R*v*conj(R) rotates v by θ about the z axis

julia> angle(R)
1.2

julia> angle(exp(6.2 * imz / 2))  # Note that this is greater than π
6.2

```
"""
Base.angle(q::Quaternion{T}) where T = 2 * absvec(log(q))
Base.angle(q::Rotor{T}) where T = 2 * absvec(log(q))


@doc raw"""
    q ^ s
    ^(q, s)

Exponentiation operator, equivalent to ``\exp(s \log(q))``.

When `s` is a real number, this is useful for natural "linear" interpolation/extrapolation
through quaternion space, starting from `q^0 = 1` and going directly to `q^1 = q`.
Specifically, when `q` is a `Rotor`, this provides a geodesic on the unit 3-sphere going
between those two points.  For more general interpolations, see the [`slerp`](@ref)
function.

"""
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
