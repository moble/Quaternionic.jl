# General math functions of quaternions

"""
    abs2(q)

Sum of the squares of all four components of the quaternion.

Note that the result for a `Rotor` is identically 1, even if that is not true numerically.

For complex components, as used to represent elements of the spacetime algebra (see
[`Lorentz`](@ref)), this is the complex "spinor norm" ``\\sum_i z_i^2``, which is a complex
number that may even vanish for a nonzero quaternion; it is *not* the Euclidean sum
``\\sum_i |z_i|^2``.  Use `norm` to obtain the Euclidean norm.

# Examples
```jldoctest
julia> abs2(quaternion(1,2,4,10))
121
```
"""
Base.abs2(q::AbstractQuaternion) = sum(abs2, components(q))
Base.abs2(q::AbstractQuaternion{T}) where {T<:Real} = sum(x->x^2, components(q))
Base.abs2(q::AbstractQuaternion{Complex{T}}) where {T<:Real} = sum(z->z*z, components(q))
Base.abs2(q::QuatVec) = sum(abs2, vec(q))
Base.abs2(q::QuatVec{T}) where {T<:Real} = sum(x->x^2, vec(q))
# For complex components, the sum is written out rather than taken over the view that `vec`
# returns; see the comment on `absvec` below.
Base.abs2(q::QuatVec{Complex{T}}) where {T<:Real} = q[2]*q[2] + q[3]*q[3] + q[4]*q[4]
Base.abs2(::Rotor{T}) where {T<:Number} = one(real(T))

# WORKAROUND for Enzyme bugs (EnzymeAD/Enzyme.jl#ISSUE_AVX).  `hypotenuse` is just `hypot`
# of three or four numbers.  Except for complex components, `abs` and `absvec` call this
# function rather than `hypot` itself, only so that `QuaternionicEnzymeExt` can give it
# rules whose arguments are the components, which are numbers.  Rules on `abs` and `absvec`
# themselves would take a quaternion argument, and Enzyme's handling of such rules crashes
# on x86_64 processors with AVX; see the extension.  Once that bug and the bug in Enzyme's
# own rule for `hypot` (also described in the extension) are fixed, `abs` and `absvec`
# should call `hypot` directly again, and this function and its rules should be removed.
hypotenuse(x, y, z) = hypot(x, y, z)
hypotenuse(w, x, y, z) = hypot(w, x, y, z)

"""
    abs(q)

Square root of the sum of the squares of all four components of the quaternion.

For real components, this function uses Julia's built-in `hypot` function to avoid overflow
and underflow.

Note that the result for a `Rotor` is identically 1, even if that is not true numerically.

For complex components, this is the principal square root of the complex spinor norm
[`abs2`](@ref), computed with a rescaling that likewise avoids overflow and underflow.  It
is a complex number in general, and may vanish for a nonzero quaternion.  Use `norm` to
obtain the real Euclidean norm.

# Examples
```jldoctest
julia> abs(quaternion(1,2,4,10))
11.0
```
"""
Base.abs(q::AbstractQuaternion) = hypotenuse(components(q)...)
Base.abs(q::AbstractQuaternion{Complex{T}}) where {T<:Real} = _hypot(components(q))
Base.abs(q::QuatVec) = hypotenuse(vec(q)...)
# See the comment on `absvec` below.
Base.abs(q::QuatVec{Complex{T}}) where {T<:Real} = _hypot(SVector{3}(q[2], q[3], q[4]))
Base.abs(::Rotor{T}) where {T<:Number} = one(real(T))

"""
    abs2vec(q)

Sum of the squares of the three "vector" components of the quaternion.

For complex components, this is the complex sum ``\\sum_i z_i^2`` over the vector
components, as with [`abs2`](@ref).

# Examples
```jldoctest
julia> abs2vec(quaternion(1,2,3,6))
49
```
"""
abs2vec(q::AbstractQuaternion) = sum(abs2, vec(q))
abs2vec(q::AbstractQuaternion{T}) where {T<:Real} = sum(x->x^2, vec(q))
# For complex components, the sum is written out rather than taken over the view that `vec`
# returns; see the comment on `absvec` below.
abs2vec(q::AbstractQuaternion{Complex{T}}) where {T<:Real} = q[2]*q[2] + q[3]*q[3] + q[4]*q[4]

"""
    absvec(q)

Square root of the sum of the squares of the three "vector" components of the quaternion.

For real components, this function uses Julia's built-in `hypot` function to avoid
overflow and underflow.  For complex components, this is the principal square root of the
complex sum [`abs2vec`](@ref), as with [`abs`](@ref).

# Examples
```jldoctest
julia> absvec(quaternion(1,2,3,6))
7.0
```
"""
absvec(q::AbstractQuaternion) = hypotenuse(vec(q)...)
# For complex components, the vector part is passed to `_hypot` as an `SVector` rather than
# as the view that `vec` returns, because Enzyme 0.13.209 aborts in reverse mode on
# reductions over that view of complex numbers (an upstream bug).  The value is the same.
absvec(q::AbstractQuaternion{Complex{T}}) where {T<:Real} = _hypot(SVector{3}(q[2], q[3], q[4]))

"""
    norm(q)

Euclidean norm of the quaternion `q`.

For real components, this equals [`abs(q)`](@ref abs); in particular, it is identically 1
for a `Rotor`.  For complex components, it is instead the real, nonnegative Euclidean norm
``\\sqrt{\\sum_i |z_i|^2}`` of the eight real numbers making up the four components, whereas
`abs` and `abs2` return the complex spinor norm.  This is the norm with which `≈`
(`isapprox`) compares quaternions with complex components.

# Examples
```jldoctest
julia> norm(quaternion(1,2,4,10))
11.0

julia> norm(Quaternion{ComplexF64}(1, 1im, 0, 0))
1.4142135623730951
```
"""
LinearAlgebra.norm(q::AbstractQuaternion) = abs(q)
LinearAlgebra.norm(q::AbstractQuaternion{<:Complex}) = norm(components(q))

Base.inv(q::AbstractQuaternion) = conj(q) / abs2(q)
function Base.inv(q::Union{Quaternion{T},QuatVec{T}}) where {T<:AbstractFloat}
    a² = abs2(q)
    if floatmin(T) ≤ a² < T(Inf)
        return conj(q) / a²
    end
    m = maximum(abs, components(q))
    if iszero(m) || !isfinite(m)
        return conj(q) / a²  # Inf or NaN, as appropriate
    end
    # `abs2` overflowed or underflowed, so rescale `q` by 2⁻ᵉ, which is exact, to bring its
    # largest component near 1.  Then inv(q) = 2⁻ᵉ inv(2⁻ᵉ q).  The factor 2⁻ᵉ itself may
    # not be representable, so `ldexp` is applied to each component instead.
    e = exponent(m)
    p = typeof(q)(ldexp.(components(q), -e))
    r = conj(p) / abs2(p)
    return typeof(r)(ldexp.(components(r), -e))
end
Base.inv(q::Rotor) = conj(q)  # Specialize to ensure output is also a Rotor


# The following functions implement the Taylor series used by `exp`, `log`, and `^` near the
# identity, along with the test for whether the series applies.  The coefficients are
# written as ratios of small integers, applied by nested multiplication and division, so
# that they are exact and cheap for every element type, including `Float16`, `BigFloat`, and
# `Double64`, and so that they do not lose accuracy when converted.

# The threshold ε^(1/8) below which the series are used, where `ε` is the machine epsilon.
# The truncation error of `atanseries` and `logseries` is then about x⁸/17 ≤ ε, while the
# general branches are kept away from the region where their derivatives lose accuracy.
seriestolerance(::Type{T}) where {T<:AbstractFloat} = sqrt(sqrt(sqrt(eps(T))))

# The threshold for `cosseries` and `sincseries`, whose truncation error is about y⁶/12!.
# For `Float64`, `Float32`, and `Float16`, ε^(1/8) is small enough that this error is below
# ε, but for types with more than about 114 bits of precision, such as `BigFloat`, it is
# not, so the threshold is reduced to 27ε^(1/6) ≤ (12! ε)^(1/6) for them.  The sixth root is
# computed with `exp` and `log`, because `cbrt(::Double64)` returns a tuple.
function trigseriestolerance(::Type{T}) where {T<:AbstractFloat}
    min(seriestolerance(T), 27 * exp(log(eps(T)) / 6))
end

# Return `true` if `|x|` is at most `tolerance(T)`, where `T` is the floating-point type of
# the value of `x`.  The test is applied to `value(x)`, so that dual numbers take the same
# branch as their values, whatever their partials are.  For types whose values are not
# floating-point numbers, such as symbolic types, the size of `x` cannot be tested, so this
# returns `iszerovalue(x)`.
issmallvalue(x, tolerance=seriestolerance) = issmallvalue(value(x), x, tolerance)
function issmallvalue(v::Union{T,Complex{T}}, _, tolerance) where {T<:AbstractFloat}
    abs(v) ≤ tolerance(T)
end
issmallvalue(_, x, _) = iszerovalue(x)

# Return `true` if the vector part of `q` is small compared to its scalar part, in the sense
# that `|abs2vec(q)| ≤ ε^(1/8) |q[1]|²`.  As with `issmallvalue`, the test is applied to the
# values of the components, and reduces to `iszerovalue(vec(q))` for other types.  The
# components are rescaled first, so that the squares cannot underflow.
isnearscalar(q::AbstractQuaternion) = isnearscalar(value.(components(q)), q)
function isnearscalar(c::SVector{4,<:Union{T,Complex{T}}}, _) where {T<:AbstractFloat}
    m = maximum(abs, c)
    iszero(m) && return true
    c = c / m
    abs(c[2]*c[2] + c[3]*c[3] + c[4]*c[4]) ≤ seriestolerance(T) * abs(c[1]*c[1])
end
isnearscalar(_, q) = iszerovalue(vec(q))

# cos(√y), accurate when `issmallvalue(y, trigseriestolerance)`
cosseries(y) = 1 - y/2 * (1 - y/12 * (1 - y/30 * (1 - y/56 * (1 - y/90))))

# sin(√y)/√y, accurate when `issmallvalue(y, trigseriestolerance)`
sincseries(y) = 1 - y/6 * (1 - y/20 * (1 - y/42 * (1 - y/72 * (1 - y/110))))

# atan(√x)/√x, accurate when `issmallvalue(x)`
function atanseries(x)
    1 - x/3 * (1 - 3x/5 * (1 - 5x/7 * (1 - 7x/9 * (1 - 9x/11 * (1 - 11x/13 * (1 - 13x/15))))))
end

# √(1+x) atan(√x)/√x, accurate when `issmallvalue(x)`.  This is v/sin(v) for the angle v
# with tan(v) = √x.  For a complex `Rotor`, `x` is complex, so its square root goes
# through `csqrt`.
logseries(x) = atanseries(x) * csqrt(1 + x)

# log(abs(q)) for a `Quaternion` with `a = abs(q)`.  When `a` is close to 1, `log(a)` loses
# relative accuracy because `a` itself is rounded, so we evaluate log1p(|q|² - 1)/2, with
# |q|² - 1 computed in a form that is accurate near |q| = 1.
function logabs(q, a)
    if abs2(a - 1) < 0.25
        log1p((q[1] - 1) * (q[1] + 1) + abs2vec(q)) / 2
    else
        log(a)
    end
end


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

The `log` function is very analogous to the [`sqrt`](@ref) function [which is essentially
`exp(log(q)/2)`], in that both functions have discontinuous behavior along the negative real
axis, and derivatives that blow up as you approach that axis.  Therefore, the same caveats
about computing derivatives using automatic differentiation apply to both functions.

Integer and rational components are converted to floating-point numbers first.  Because the
algorithm chooses its branch by comparing component values, symbolic component types such as
`Symbolics.Num` instead get the general formula ``\log|q| + \arctan(|\vec{v}|, w)\,
\vec{v}/|\vec{v}|`` (from the Symbolics extension), which should not be evaluated on the
real axis by direct substitution.

If we decompose the result of this function as ``\log(q) = s + \vec{v}``, where ``s`` is the
scalar part and ``\vec{v}`` is the pure-vector part, then clearly ``s`` and ``\vec{v}``
commute, so their exponentials also commute, and we have
```math
q = \exp\left(\log(q)\right) = \exp\left(s\right) \exp\left(\vec{v}\right).
```
Note that the exponential of a pure-vector quaternion is a unit quaternion, so we have
decomposed ``q`` into a product of a positive real number and a unit quaternion.  We already
know how to take the logarithm of a positive real number; when ``|q|`` is close to 1, we
compute it as `log1p(abs2(q) - 1)/2`, with `abs2(q) - 1` evaluated without cancellation, for
accuracy.  Therefore, in the following we assume that ``q`` is a unit quaternion, so that
``s=0`` and we only need to find ``\vec{v}``.

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
accurately by `invsinc(v)` whenever ``|v| \leq \pi/2``.  Near the identity, where the vector
part is small compared to the scalar part, ``\sin(v)`` is the square root of a small number,
whose derivatives lose accuracy; there, we instead evaluate ``v/\sin(v)`` as a Taylor series
in ``x = |\vec{q}|^2 / q[1]^2 = \tan^2(v)``, so that automatic differentiation remains
accurate.

When ``|v| > \pi/2`` (which also implies ``\cos(v)<0``), we instead define ``v' = \pi -
v``, so that ``\cos(v') = -\cos(v)`` and ``\sin(v') = \sin(v)``, and we can rewrite the
problem as
```math
\exp\left(v \hat{v}\right) = \cos(v) + \hat{v} \sin(v) = -\cos(v') + \hat{v} \sin(v').
```
We can compute ``v' = \mathrm{atan}(\sin(v), -\cos(v))`` accurately, and then we need to
multiply the unit vector ``\hat{v}`` by ``v = \pi - v'``.  This algorithm is surprisingly
accurate, even when ``v`` is extremely close to ``\pi``, which implies that the vector part
of the input is extremely small.

The only special case remaining to handle is when ``\cos(v) < 0`` but ``\sin(v)`` is
*identically* zero.  In this case, we could throw an error, but this is not usually helpful.
Instead, we arbitrarily choose to return ``\pi 𝐤``.

For complex components, the same formulas are used with the complex spinor norm
[`abs`](@ref) in place of the Euclidean norm.  The branch is then chosen by the sign of the
real part of ``\cos(v) = q[1]/|q|``.  When it is nonnegative, the angle is computed as ``v =
-i \log\left(\cos(v) + i \sin(v)\right)``, where ``i`` is the imaginary unit of the complex
components and the logarithm is the principal complex logarithm; this is smooth even where
``q[1] = 0``.  Otherwise, ``v'`` is computed with the principal complex branch of `atan`.

If `q` is a `Rotor`, we return a `QuatVec`; if `q` is a general `Quaternion`, we return a
general `Quaternion` — though if `q` has unit norm, the scalar part of the result will be
zero up to roundoff.  Note that, because there is no geometric reason to take the logarithm
of a `QuatVec`, that case is not implemented; if you really need to compute it, you can
convert the `QuatVec` to a `Quaternion` first.

"""
function Base.log(q::Quaternion)
    q = float(q)
    T = basetype(q)
    if iszerovalue(q)  # q == 0
        return Quaternion{T}(-Inf, false, false, false)
    end
    # When the largest component is beyond about 1e154 or below about 1e-146 (for Float64),
    # the values are still right, but the derivatives of `x / y` and `atan(y, x)` that AD
    # tools use contain `y²` and `x² + y²`, which overflow or underflow.  In that case `q`
    # is multiplied by a power of two σ = 2ᵏ that brings its largest component near 1, and
    # log(q) = log(σq) - k log(2) is returned.  This scaling is exact, and the vector part
    # of the logarithm does not depend on |q|.  As in division, the test uses only the
    # values of the components, so that dual numbers are rescaled in the same way as the
    # floats they contain, and the power is clamped to the range of normal numbers, so that
    # the recursive call scales again if one step was not enough.
    if value(q[1]) isa AbstractFloat
        R = typeof(value(q[1]))
        m = maximum(x -> abs(value(x)), components(q))
        if isfinite(m) && !(sqrt(floatmin(R) / eps(R)) ≤ m ≤ sqrt(floatmax(R)) / 4)
            k = clamp(-exponent(m), exponent(floatmin(R)), exponent(floatmax(R)))
            return log(ldexp(one(R), k) * q) - k * log(R(2))
        end
    end
    a = abs(q)
    cosv = q[1]
    if cosv ≥ 0  # q[1] ≥ 0
        f = if isnearscalar(q)
            # Near the identity, `absvec` is the square root of a tiny number, whose
            # derivatives are huge and lose all accuracy when they cancel, or are infinite
            # when it is exactly 0; use a series in tan²(v) instead.
            logseries(abs2vec(quatvec(q) / cosv))
        else
            invsinc(atan(absvec(q), cosv))
        end
        return logabs(q, a) + f * (quatvec(q) / a)
    elseif iszerovalue(vec(q))  # q is a negative real number
        # Note that we check this branch only after ruling out cosv≥0 because this could
        # otherwise correspond to *positive* real numbers, which are treated correctly by
        # the preceding branch, but only the preceding branch will behave correctly for AD.
        return Quaternion{T}(logabs(q, a), false, false, π)
    else  # q[1] < 0 but q⃗ ≠ 0
        sinv = absvec(q)
        v′ = atan(sinv, -cosv)
        # Dividing the vector part by `sinv` first avoids overflow when `sinv` is tiny.
        return logabs(q, a) + (π - v′) * (quatvec(q) / sinv)
    end
end
function Base.log(q::Rotor)
    q = float(q)
    cosv = q[1]
    if cosv ≥ 0  # q[1] ≥ 0
        f = if isnearscalar(q)
            # Series branch required for AD; see `log(::Quaternion)` above
            logseries(abs2vec(quatvec(q) / cosv))
        else
            invsinc(atan(absvec(q), cosv))
        end
        return f * quatvec(q)
    elseif iszerovalue(vec(q))  # q is a negative real number
        # Note that we check this branch only after ruling out cosv≥0 because this could
        # otherwise correspond to *positive* real numbers, which are treated correctly by
        # the preceding branch, but only the preceding branch will behave correctly for AD.
        return QuatVec{basetype(q)}(false, false, false, π)
    else  # q[1] < 0 but q⃗ ≠ 0
        sinv = absvec(q)
        v′ = atan(sinv, -cosv)
        return (π - v′) * (quatvec(q) / sinv)
    end
end

function Base.log(q::Quaternion{<:Complex})
    q = float(q)
    T = basetype(q)
    if iszerovalue(q)
        return Quaternion{T}(T(-Inf), false, false, false)
    end
    a = abs(q)
    cosv = q[1]
    # The first branch below returns the angle v with -π/2 ≤ real(v) ≤ π/2, whose cosine has
    # a nonnegative real part, so the branch must be chosen by cos(v) = q[1]/a, not by q[1]
    # alone.
    if real(cosv / a) ≥ 0
        f = if isnearscalar(q)
            # Series branch required for AD; see `log(::Quaternion)` above
            logseries(abs2vec(quatvec(q) / cosv))
        else
            # The angle v satisfies cos(v) = cosv/a and sin(v) = sinv/a, so that exp(iv) =
            # (cosv + i sinv)/a, where `i` is the imaginary unit of the components.  The
            # principal logarithm gives the v with -π/2 ≤ real(v) ≤ π/2 in this branch.
            # Unlike atan(sinv/cosv), this is smooth where cosv = 0, so its derivatives are
            # correct there too.
            sinv = absvec(q)
            v = -im * log((cosv + im * sinv) / a)
            v * a / sinv
        end
        return logabs(q, a) + f * (quatvec(q) / a)
    elseif iszerovalue(vec(q))
        return Quaternion{T}(logabs(q, a), false, false, T(π))
    else
        sinv = absvec(q)
        v′ = atan(-sinv / cosv)
        f = -invsinc(v′) * (v′ - π) / v′
        return logabs(q, a) + f * (quatvec(q) / a)
    end
end

function Base.log(q::Rotor{<:Complex})
    q = float(q)
    T = basetype(q)
    cosv = q[1]
    if real(cosv) ≥ 0
        f = if isnearscalar(q)
            # Series branch required for AD; see `log(::Quaternion)` above
            logseries(abs2vec(quatvec(q) / cosv))
        else
            # As in `log(::Quaternion{<:Complex})` above, with a = 1
            sinv = absvec(q)
            v = -im * log(cosv + im * sinv)
            v / sinv
        end
        return f * quatvec(q)
    elseif iszerovalue(vec(q))
        return QuatVec{T}(false, false, T(π))
    else
        sinv = absvec(q)
        v′ = atan(-sinv / cosv)
        f = -invsinc(v′) * (v′ - π) / v′
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

The exponential of a `QuatVec` is a `Rotor`.  The exponential of a `Quaternion` or of a
`Rotor` is a `Quaternion`.

!!! note "Derivatives near 0"

    When the vector part of `q` is small, its contribution is evaluated as a Taylor series
    so that automatic differentiation is accurate there.  Derivatives with respect to the
    vector components *at exactly* zero vector part are therefore accurate only through
    eleventh order.  Using twelfth-order derivatives or higher is so unusual that this will
    likely not be a problem in practice.

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
if we do nothing but return the value.  Moreover, when ``v`` is small but nonzero, ``v`` is
the square root of a small number, and higher derivatives computed through that square root
and through ``\sin(v)/v`` lose accuracy to cancellation.  Therefore, whenever ``v^2`` is
small, we evaluate ``\cos(v)`` and ``\sin(v)/v`` as Taylor series in ``v^2``, which are
accurate to machine precision there, and give accurate derivatives.

With symbolic component types such as `Symbolics.Num`, the vector part cannot be tested for
being small, so the general formula is returned.  That expression has a removable
singularity where the vector part vanishes, so it should not be evaluated there by direct
substitution.
"""
function Base.exp(q::Quaternion)
    q = float(q)
    e = cexp(q[1])
    a² = abs2vec(q)
    c, s = if issmallvalue(a², trigseriestolerance)
        cosseries(a²), sincseries(a²)
    else
        a = absvec(q)
        ccos(a), _sincu(a)
    end
    ec, es = e*c, e*s
    Quaternion{typeof(ec)}(ec, es*q[2], es*q[3], es*q[4])
end
function Base.exp(v⃗::QuatVec)
    v⃗ = float(v⃗)
    a² = abs2vec(v⃗)
    c, s = if issmallvalue(a², trigseriestolerance)
        cosseries(a²), sincseries(a²)
    else
        a = absvec(v⃗)
        ccos(a), _sincu(a)
    end
    Rotor{typeof(c)}(c, s*v⃗[2], s*v⃗[3], s*v⃗[4])
end
Base.exp(q::Rotor) = exp(quaternion(q))

@doc raw"""
    sqrt(q)

Square root of a quaternion.

The square root of a `Rotor` is a `Rotor`, and the square root of a `Quaternion` or of a
`QuatVec` is a `Quaternion`.  (The square root of a pure vector is not a pure vector.)
Integer and rational components are converted to floating-point numbers first.

!!! note "Branch-cut behavior"

    As with the logarithm, the quaternionic square-root has a branch cut along the
    non-positive real axis: if the vector components of `q` are precisely zero *and* the
    scalar component is negative, the returned quaternion will be `√(-q[1]) * 𝐤`.  Unlike in
    the complex case, the choice of `𝐤` is arbitrary; the vector component of the returned
    quaternion could be any unit vector.

    For complex components, when the vector components are precisely zero and the real
    part of the scalar component is not positive, the result is instead the scalar
    `sqrt(q[1])`, computed as a complex number, with zero vector part.

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
\sqrt{q} = \frac{|q⃗|}{\sqrt{2(|q| − q[1])}} + \frac{q⃗}{|q⃗|} \sqrt{\frac{|q| − q[1]}{2}}.
```

For real floating-point components, a quaternion whose norm is so large or so small that
these expressions would overflow or underflow is first rescaled by a power of 4, and the
result is rescaled by the corresponding power of 2.  In the second expression, the unit
vector ``q⃗/|q⃗|`` is computed from the unscaled vector part, after dividing it by its
largest component, so that it keeps its precision even when the vector part is tiny compared
to the scalar part.  For complex components, the first expression is used whenever the real
part of `q[1]` is nonnegative, and the second with ``|q⃗|^2`` evaluated as the complex sum
`abs2vec(q)` otherwise.

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

Because the algorithm chooses its branch by comparing component values, it does not work
with symbolic component types such as `Symbolics.Num`.
"""
function Base.sqrt(q::Union{Quaternion{T},Rotor{T}}) where {T<:Real}
    q = float(q)
    Q = typeof(q)
    w = q[1]
    if w ≤ 0 && iszerovalue(vec(q))
        return Q(false, false, false, √(-w))
    end
    k = sqrtscale(q)
    if !iszero(k)
        return rescaledsqrt(q, k)
    end
    # `ifelse` would evaluate both arms, and the unused arm is 0/0 whenever `q` is a
    # positive real; the resulting NaN contaminates the derivative.  Branch instead, so
    # that arm is never evaluated.  As of 2026-08-19 Mooncake and ReverseDiff both still
    # require this.
    a = abs(q)
    if w ≥ 0
        c = √(2(a + w))
        return Q(c/2, q[2]/c, q[3]/c, q[4]/c)
    else
        # Here the vector part is nonzero.  Its direction q⃗/|q⃗| is computed after dividing
        # q⃗ by its largest component `m`, so that |q⃗|/m = `s` neither overflows nor loses
        # precision to underflow, even when the vector part is tiny or subnormal.  The norm
        # `s` is computed as `abs` of a `QuatVec`, which is `hypotenuse` of the components,
        # so that Enzyme uses the rule for `hypotenuse` in QuaternionicEnzymeExt: Enzyme's
        # own reverse rule for `hypot` with three arguments fails at batch widths above 1.
        m = maximum(abs, vec(q))
        u = vec(q) / m
        s = abs(quatvec(u...))
        d = √(2(a - w))
        h = d/2
        return Q((m/d)*s, (u[1]/s)*h, (u[2]/s)*h, (u[3]/s)*h)
    end
end
function Base.sqrt(q::Union{Quaternion{Float16},Rotor{Float16}})
    T = wrapper(typeof(q))
    T{Float16}(sqrt(T{Float32}(q)))
end
function Base.sqrt(q::Union{Quaternion{Complex{T}},Rotor{Complex{T}}}) where {T<:Real}
    q = float(q)
    Q = typeof(q)
    # As in the real case, only a vanishing vector part with `real(q[1]) ≤ 0` needs a
    # special branch, because the general formula below would divide by zero.  Elsewhere,
    # the general formula is correct, and taking this branch would lose the derivative with
    # respect to the vector part.
    if real(q[1]) ≤ 0 && iszerovalue(vec(q))
        return Q(csqrt(q[1]), false, false, false)
    end
    c₁ = if real(q[1]) ≥ 0
        (abs(q) + q[1])
    else
        (abs2vec(q) / (abs(q) - q[1]))
    end
    c₂ = csqrt(inv(2c₁))
    return Q(c₁*c₂, q[2]*c₂, q[3]*c₂, q[4]*c₂)
end
Base.sqrt(q::QuatVec) = sqrt(quaternion(q))

# The power `k` such that `sqrt` should rescale `q` by 4⁻ᵏ to avoid overflow and underflow,
# or 0 if no rescaling is needed.  The test uses the largest component, which is finite even
# when `abs(q)` overflows.  A `Rotor` and types other than floating-point numbers (such as
# dual numbers) are not rescaled.
function sqrtscale(q::Quaternion{T}) where {T<:AbstractFloat}
    m = maximum(abs, components(q))
    if isfinite(m) && !iszero(m) && !(floatmin(T)/eps(T) ≤ m ≤ floatmax(T)/8)
        exponent(m) ÷ 2
    else
        0
    end
end
sqrtscale(_) = 0

# The square root of a floating-point `Quaternion` that `sqrt` rescales by 4⁻ᵏ, which is
# exact, so that |q| is near 1, and then unscales by 2ᵏ.  The factor 4⁻ᵏ itself may not be
# representable, so `q` is multiplied twice by 2⁻ᵏ.  The scale factors are constants made
# with `ldexp`, and the components are only multiplied by them, because applying `ldexp` to
# the components themselves breaks Enzyme's forward mode for every call to `sqrt`.  The
# magnitudes `c` and `d` below are computed from the rescaled quaternion `p`, but the vector
# part of `p` may underflow when it is much smaller than the scalar part, so the vector part
# of the result is computed from that of `q` instead.
function rescaledsqrt(q::Quaternion{T}, k) where {T<:AbstractFloat}
    t = ldexp(one(T), -k)  # 2⁻ᵏ
    t⁻¹ = ldexp(one(T), k)  # 2ᵏ
    p = (q * t) * t
    if q[1] ≥ 0
        c = √(2(abs(p) + p[1]))  # 2⁻ᵏ times the `c` in `sqrt`
        return Quaternion{T}(t⁻¹ * (c/2), (t*q[2])/c, (t*q[3])/c, (t*q[4])/c)
    else
        # See the comments on the unscaled case in `sqrt`.
        d = √(2(abs(p) - p[1]))  # 2⁻ᵏ times the `d` in `sqrt`
        m = maximum(abs, vec(q))
        u = vec(q) / m
        s = abs(quatvec(u...))
        h = d/2
        return Quaternion{T}(
            ((t*m) / d) * s, t⁻¹ * ((u[1]/s)*h), t⁻¹ * ((u[2]/s)*h), t⁻¹ * ((u[3]/s)*h)
        )
    end
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

For real components, this is computed as `2atan(absvec(q), q[1])`, which is equivalent to
`2absvec(log(q))` but is accurate, along with its derivatives, even near `-1`.  The angle
has cone-shaped kinks where the vector part vanishes, at which derivatives are returned as
zero.

For complex components, as for a [`Lorentz`](@ref) transformation, the result is the complex
number `2absvec(log(q))`.  When the rotation and the boost share an axis, this is ``θ ±
iη``, where ``θ`` is the rotation angle, ``η`` is the rapidity of the boost, and the sign is
positive when the boost is parallel to the rotation axis.  In general, the result does not
separate so simply into a rotation angle and a rapidity, and the imaginary part may be
negative.

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
function Base.angle(q::Union{Quaternion{<:Real},Rotor{<:Real}})
    if iszerovalue(vec(q))
        # `absvec` is not differentiable here; return the same value with zero derivative.
        return 2atan(zero(q[1]), q[1])
    end
    2atan(absvec(q), q[1])
end
Base.angle(q::Union{Quaternion{<:Complex},Rotor{<:Complex}}) = 2 * absvec(log(q))


@doc raw"""
    q ^ s
    ^(q, s)

Exponentiation operator, equivalent to ``\exp(s \log(q))``.

When `s` is a real number, this is useful for natural "linear" interpolation/extrapolation
through quaternion space, starting from `q^0 = 1` and going directly to `q^1 = q`.
Specifically, when `q` is a `Rotor`, this provides a geodesic on the unit 3-sphere going
between those two points.  For more general interpolations, see the [`slerp`](@ref)
function.

The type of the result depends on the types of the arguments:

  * A `Rotor` raised to a real or complex power is a `Rotor`.  The norm of the input is
    assumed to be 1, so the result always has unit norm.  This includes `Lorentz` rotors,
    with complex components.
  * A `Quaternion` or `QuatVec` raised to a non-integer power is a `Quaternion`.
  * When `s` is itself a quaternion, the result is the `Quaternion` ``\exp(s \log(q))``,
    with the product taken in that order.
  * When `s` is an `Integer`, the power is computed by repeated multiplication.  Then a
    `Rotor` gives a `Rotor`, and a `Quaternion` or `QuatVec` gives a `Quaternion`.  With
    integer components, a nonnegative power keeps integer components, and a negative power
    gives floating-point components.

The branch cut is the same as for [`log`](@ref): a quaternion whose vector part is exactly
zero and whose scalar part is negative is treated as ``|q| \exp(π𝐤)``.  Just as for `log`,
derivatives are not defined there for non-integer `s`, and will be returned as zero with
respect to `q`.

Non-integer powers choose their branch by comparing component values, so with symbolic
component types such as `Symbolics.Num` they are evaluated as `exp(s * log(q))`, with the
general formula for `log`.  Integer powers work with them directly.
"""
function Base.:^(q::Quaternion, s::Number)
    exp(s * log(q))
end
function Base.:^(q::Rotor, s::Number)
    q = float(q)
    w = q[1]
    if real(w) > 0 && isnearscalar(q)
        # Near the identity, `absvec(q)` is the square root of a tiny number, whose
        # derivatives are huge and lose all accuracy when they cancel, or are NaN when it is
        # exactly 0.  Instead, we write log(q) = A 𝐯/q[1], where A = atan(√x)/√x and x =
        # |𝐯|²/q[1]², and q^s = exp(f 𝐯) with f = s A/q[1], using series in `x` and in y =
        # |f 𝐯|² that are smooth at the identity.  Note that `x` is computed only here, so
        # that no division by a vanishing q[1] is ever recorded by reverse-mode AD.
        v² = abs2vec(q)
        A = atanseries((v² / w) / w)
        f = s * A / w
        y = f * f * v²
        c, sincu = if issmallvalue(y, trigseriestolerance)
            cosseries(y), sincseries(y)
        else
            # Reached when |s| is large enough that √y is not small, or, for types with
            # high precision, when the series in `y` would not be accurate enough.  For a
            # complex `Rotor`, `y` is complex, so its square root goes through `csqrt`.
            a = csqrt(y)
            ccos(a), csin(a) / a
        end
        return Rotor{typeof(c)}(c, sincu*f*q[2], sincu*f*q[3], sincu*f*q[4])
    end
    power(q, s)
end
Base.:^(q::Rotor, s::AbstractQuaternion) = quaternion(q)^s

# Rational exponents would otherwise be ambiguous with Base's `^(::Number, ::Rational)`.  We
# forward the exponent unchanged, rather than converting it to a float, so that the
# precision of types like `BigFloat` is kept.
Base.:^(q::Quaternion, s::Rational) = exp(s * log(q))
Base.:^(q::Rotor, s::Rational) = invoke(^, Tuple{Rotor,Number}, q, s)
Base.:^(q::QuatVec, s::Rational) = quaternion(q)^s

# The power of a real `Rotor` away from the identity
function power(q::Rotor{<:Real}, s)
    if iszerovalue(vec(q))
        # q is -1, so log(q) = π𝐤.  We promote `s` (as in the general case below) rather than
        # converting it, so that the result has the same type as in the other branches, and
        # so that dual numbers work.
        sπ = s * one(basetype(q))
        sin_πs, cos_πs = sinpi(sπ), cospi(sπ)
        return Rotor{typeof(cos_πs)}(cos_πs, false, false, sin_πs)
    end
    absolutevec = absvec(q)
    f1 = s * atan(absolutevec, q[1])
    # `sin` and `cos` are called separately, rather than through `sincos`, because Enzyme
    # cannot take second derivatives of `sincos`.
    sin_f1, cos_f1 = sin(f1), cos(f1)
    # Dividing each component by `absolutevec` before multiplying avoids overflow when
    # `absolutevec` is tiny.
    Rotor{typeof(cos_f1)}(
        cos_f1, sin_f1*(q[2]/absolutevec), sin_f1*(q[3]/absolutevec), sin_f1*(q[4]/absolutevec)
    )
end

# The power of a complex `Rotor` (such as a `Lorentz` rotor) away from the identity
power(q::Rotor{<:Complex}, s) = exp(s * log(q))

# We need to be more specific about the quaternion types here because we had
# to be specific about the quaternion type above, and Integer<:Number
for QT ∈ (AbstractQuaternion, Quaternion, Rotor)
    @eval function Base.:^(q::$QT, s::Integer)
        s < 0 && return inv(q ^ -s)
        # This is the square-and-multiply algorithm of `Base.power_by_squaring`, written out
        # here, with the same order of operations and therefore the same results and types.
        # As of 2026-10-01, Enzyme's reverse mode aborts the Julia process when
        # differentiating `Base.power_by_squaring` with the exponents 1 and 2.  As in Base's
        # `to_power_type`, the input is first converted to the type of a product, so that
        # the type of the result does not depend on `s`.  For example, that type has
        # floating-point components for a `Rotor`, because the product of two rotors is
        # normalized, and integer components for `Bool` components, as in `imx`.
        x = convert(Base.promote_op(*, typeof(q), typeof(q)), q)
        s == 0 && return one(x)
        while iseven(s)
            x *= x
            s >>= 1
        end
        p = x
        while (s >>= 1) > 0
            x *= x
            if isodd(s)
                p *= x
            end
        end
        p
    end
end
Base.:^(q::QuatVec, s::Integer) = quaternion(q)^s

# As with the arithmetic operators in `algebra.jl`, `@fastmath q ^ s` would otherwise
# promote `s` to a quaternion, which changes the result for a `Quaternion` and throws an
# error for a `Rotor`.  (See issue #46.)  The `Integer` method is needed to resolve an
# ambiguity with `Base.FastMath`.
Base.FastMath.pow_fast(q::AbstractQuaternion, s::Number) = q ^ s
Base.FastMath.pow_fast(q::AbstractQuaternion, s::Integer) = q ^ s
