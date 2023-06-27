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
