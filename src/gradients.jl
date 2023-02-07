@doc raw"""
    ∂log(Z::Rotor)

Return the gradient of `log(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even
though change of `Z` in a direction that changes its norm would not be allowed
for a `Rotor`, we measure the gradient in that direction anyway.  That is, the
elements of the returned vector of quaternions is
```math
\begin{aligned}
  \left[
    \frac{\partial} {\partial Z_w} \log(Z),
    \frac{\partial} {\partial Z_x} \log(Z),
    \frac{\partial} {\partial Z_y} \log(Z),
    \frac{\partial} {\partial Z_z} \log(Z)
  \right].
\end{aligned}
```

Note that, even though `log(::Rotor)` is a `QuatVec`, the derivative (and
therefore each element of the result) is a general `Quaternion`.

See also [`∂exp`](@ref) for a similar function, as well as [`log∂log`](@ref)
for a function to compute the value along with the gradient.

# Examples
```julia
julia> ∂log∂w, ∂log∂x, ∂log∂y, ∂log∂z = ∂log(randn(QuatVecF64));

```
"""
function ∂log(Z::Rotor)
    # log(Z::Rotor) = log(abs2) / 2 + f(Z) * vec(Z)
    # f(Z) = atan(absvec / Z[1]) / absvec = acos(Z[1]) / absvec
    a2 = abs2vec(Z)
    a = sqrt(a2)
    if a < 2eps(typeof(a))
        return Quaternion.([Quaternion(one(a)), zero(a)*imx, zero(a)*imy, one(a)*imz])
    end
    f = acos(Z[1]) / a
    fprime = (Z[1] - f) / a2
    V = 1 + fprime * QuatVec(Z)
    [
        Quaternion(conj(Z)),
        Z[2] * V + f * imx,
        Z[3] * V + f * imy,
        Z[4] * V + f * imz,
    ]
end


"""
    log∂log(Z::Rotor)

Return the value and gradient of `log(Z)` with respect to the components of `Z`.

See [`∂log`](@ref) for more explanation of the components of the gradient.

# Examples
```julia
julia> l, ∂l = log∂log(randn(RotorF64));

```
"""
function log∂log(Z::Rotor)
    # log(Z::Rotor) = log(abs2) / 2 + f(Z) * vec(Z)
    # f(Z) = atan(absvec / Z[1]) / absvec = acos(Z[1]) / absvec
    a2 = abs2vec(Z)
    a = sqrt(a2)
    if a < 2eps(typeof(a))
        return (
            QuatVec(zero(a)),
            Quaternion.([Quaternion(one(a)), zero(a)*imx, zero(a)*imy, one(a)*imz])
        )
    end
    f = acos(Z[1]) / a
    fprime = (Z[1] - f) / a2
    V = 1 + fprime * QuatVec(Z)
    ∂log = [
        Quaternion(conj(Z)),
        Z[2] * V + f * imx,
        Z[3] * V + f * imy,
        Z[4] * V + f * imz,
    ]
    (f * QuatVec(Z), ∂log)
end


@doc raw"""
    ∂exp(Z::QuatVec)

Return the gradient of `exp(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even
though a scalar component of `Z` would not be allowed for a `QuatVec`, we
measure the gradient in that direction anyway.  That is, the first element
of the returned vector of quaternions is
```math
\begin{aligned}
  \left.\frac{\partial} {\partial Z_w} \exp(Z) \right|_{Z_w=0}.
\end{aligned}
```

Note that, even though `exp(::QuatVec)` is a `Rotor`, the derivative (and
therefore each element of the result) is a general `Quaternion`.

See also [`∂log`](@ref) for a similar function, as well as [`exp∂exp`](@ref)
for a function to compute the value along with the gradient.

# Examples
```julia
julia> ∂exp∂w, ∂exp∂x, ∂exp∂y, ∂exp∂z = ∂exp(randn(QuatVecF64));

```
"""
function ∂exp(Z::QuatVec)
    # exp(Z::Quaternion) = exp(Z[1]) * cos(absvec) + g(Z) * vec(Z)
    # g(Z) = exp(Z[1]) * sin(absvec) / absvec
    a2 = abs2vec(Z)
    a = sqrt(a2)
    if a < 2eps(typeof(a))
        return Quaternion.([Quaternion(one(a)), one(a)*imx, one(a)*imy, one(a)*imz])
    end
    s, c = sincos(a)
    sinca = s / a
    g = sinca
    gprime = (c - g) / a2
    V = -sinca + gprime * Z
    [
        c + sinca * Z,
        Z[2] * V + g * imx,
        Z[3] * V + g * imy,
        Z[4] * V + g * imz,
    ]
end


"""
    exp∂exp(Z::QuatVec)

Return the value and gradient of `exp(Z)` with respect to the components of `Z`.

See [`∂exp`](@ref) for more explanation of the components of the gradient.

# Examples
```julia
julia> e, ∂e = exp∂exp(randn(QuatVecF64));

```
"""
function exp∂exp(Z::QuatVec)
    # exp(Z::Quaternion) = exp(Z[1]) * cos(absvec) + g(Z) * vec(Z)
    # g(Z) = exp(Z[1]) * sin(absvec) / absvec
    a2 = abs2vec(Z)
    a = sqrt(a2)
    if a < 2eps(typeof(a))
        return (
            Quaternion(one(a)),
            Quaternion.([Quaternion(one(a)), one(a)*imx, one(a)*imy, one(a)*imz])
        )
    end
    s, c = sincos(a)
    sinca = s / a
    g = sinca
    gprime = (c - g) / a2
    V = -sinca + gprime * Z
    ∂exp = [
        c + sinca * Z,
        Z[2] * V + g * imx,
        Z[3] * V + g * imy,
        Z[4] * V + g * imz,
    ]
    (c+g*Z, ∂exp)
end


@doc raw"""
    slerp∂slerp(q₁, q₂, τ)

Return the value and gradient of `slerp`.

The gradient is with respect to each of the input arguments in turn, with each
quaternion regarded as a series of four arguments.  That is, a total of 10
quaternions will be returned:
```math
\begin{aligned}
  &\big[\\
    &\quad \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₁.w} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₁.x} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₁.y} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₁.z} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₂.w} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₂.x} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₂.y} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial q₂.z} \mathrm{slerp}(q₁, q₂, τ), \\
    &\quad \frac{\partial}{\partial \tau} \mathrm{slerp}(q₁, q₂, τ), \\
  &\big]
\end{aligned}
```
For convenience, this will be a 4-tuple with the `slerp` as the first element,
the first four components of the derivative, followed by the next four
components of the derivative, followed by the last component of the derivative.

See also [`slerp`](@ref) for just the value, and [`slerp∂slerp∂τ`](@ref) for
just the value and the derivative with respect to `τ`.

# Examples
```julia
julia> (q1, q2), τ = randn(RotorF64, 2), rand();

julia> s, ∂s∂q₁, ∂s∂q₂, ∂s∂τ = slerp∂slerp(q₁, q₂, τ);

```
"""
function slerp∂slerp(q₁::Rotor, q₂::Rotor, τ::Real)
    r = q₂ / q₁
    basis = typeof(r)[1, imx, imy, imz]
    l, ∂l = log∂log(r)
    e, ∂e = exp∂exp(τ*l)
    s = e * q₁
    (
        s,
        [
            sum(
                ∂e[b] * τ * ∂l[c][b] * (
                    q₂ * conj(basis[a]) #/ abs2(q₁)
                    - 2q₁[a] * (q₂/q₁) #/ abs2(q₁)
                )[c] * q₁
                for b in 1:4
                for c in 1:4
            ) + e * basis[a]
            for a in 1:4
        ],
        [
            sum(
                ∂e[b] * τ * ∂l[c][b] * (
                    basis[a] / q₁
                )[c] * q₁
                for b in 1:4
                for c in 1:4
            )
            for a in 1:4
        ],
        l * s
    )
end


"""
    slerp∂slerp∂τ(q₁, q₂, τ)

Return the value and time-derivative of `slerp`.

See also [`slerp∂slerp`](@ref), which returns the value and *all* of the
derivatives of `slerp`.

"""
function slerp∂slerp∂τ(q₁::Rotor, q₂::Rotor, τ::Real)
    l = log(q₂ / q₁)
    e = exp(τ*l)
    s = e * q₁
    ∂s∂t = l * s
    (s, ∂s∂t)
end


"""
    squad∂squad∂t(qᵢ, A, B, qᵢ₊₁, ta, tb, t)

Compute the value and time-derivative of [`squad`](@ref).

This is primarily an internal helper function, taking various parameters
computed within the `squad` function.  This will be used to compute the
derivative of `squad` when the angular velocity is also requested.  To actually
obtain the derivative, simply pass the relevant keyword to the `squad`
function.

"""
function squad∂squad∂t(qᵢ, A, B, qᵢ₊₁, ta, tb, t)
    # squad = slerp(slerp(qᵢ, qᵢ₊₁, τ), slerp(A, B, τ), 2τ*(1-τ))
    # τ = (t - ta) / (tb - ta)
    # ∂τ∂t = 1 / (tb - ta)
    # ∂(2τ*(1-τ))/∂t = ∂τ/∂t * (2 - 4τ)

    τ = (t - ta) / (tb - ta)
    ∂τ∂t = 1 / (tb - ta)

    X, ∂X∂τ = slerp∂slerp∂τ(qᵢ, qᵢ₊₁, τ)
    ∂X∂t = ∂τ∂t * ∂X∂τ

    Y, ∂Y∂τ = slerp∂slerp∂τ(A, B, τ)
    ∂Y∂t = ∂τ∂t * ∂Y∂τ

    s, ∂s∂X, ∂s∂Y, ∂s∂τparabola = slerp∂slerp(X, Y, 2τ*(1-τ))
    ∂s∂t = ∂τ∂t * (2 - 4τ) * ∂s∂τparabola

    (
        s,
        sum(∂X∂t[b] * ∂s∂X[b] for b in 1:4)
        + sum(∂Y∂t[b] * ∂s∂Y[b] for b in 1:4)
        + ∂s∂t
    )
end
