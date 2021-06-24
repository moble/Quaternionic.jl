@doc raw"""
    ∂log(Z::Rotor)

Return the gradient of `log(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even
though change of `Z` in a direction that changes its norm would not be allowed
for a `Rotor`, we measure the gradient in that direction anyway.  That is, the
elements of the returned vector of quaternions is
```math
\begin{equation}
  \left[
    \frac{\partial} {\partial Z_w} log(Z),
    \frac{\partial} {\partial Z_x} log(Z),
    \frac{\partial} {\partial Z_y} log(Z),
    \frac{\partial} {\partial Z_z} log(Z)
  \right].
\end{equation}
```

Note that, even though `log(::Rotor)` is a `QuatVec`, the derivative (and
therefore each element of the result) is a general `Quaternion`.

See also [`∂exp`](@ref).

"""
function ∂log(Z::Rotor)
    # log(Z::Rotor) = log(abs2) / 2 + f(Z) * Z.vec
    # f(Z) = atan(absvec / Z.w) / absvec = acos(Z.w) / absvec
    a2 = abs2vec(Z)
    a = sqrt(a2)
    if a < 2eps(typeof(a))
        return Quaternion.([Quaternion(one(a)), zero(a)*imx, zero(a)*imy, one(a)*imz])
    end
    f = acos(Z.w) / a
    fprime = (Z.w - f) / a2
    V = 1 + fprime * QuatVec(Z)
    [
        Quaternion(conj(Z)),
        Z.x * V + f * imx,
        Z.y * V + f * imy,
        Z.z * V + f * imz,
    ]
end


@doc raw"""
    ∂exp(Z::QuatVec)

Return the gradient of `exp(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even
though a scalar component of `Z` would not be allowed for a `QuatVec`, we
measure the gradient in that direction anyway.  That is, the first element
of the returned vector of quaternions is
```math
\begin{equation}
  \left.\frac{\partial} {\partial Z_w} exp(Z) \right|_{Z_w=0}.
\end{equation}
```

Note that, even though `exp(::QuatVec)` is a `Rotor`, the derivative (and
therefore each element of the result) is a general `Quaternion`.

See also [`∂log`](@ref).

"""
function ∂exp(Z::QuatVec)
    # exp(Z::Quaternion) = exp(Z.w) * cos(absvec) + g(Z) * Z.vec
    # g(Z) = exp(Z.w) * sin(absvec) / absvec
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
        Z.x * V + g * imx,
        Z.y * V + g * imy,
        Z.z * V + g * imz,
    ]
end


println("***NOTE***: I'm not sure about the order of the indices on `∂l` in the `∂slerp∂q` functions.")


function ∂slerp∂q₁(q₁, q₂, τ)
    basis = Rotor[Quaternion(1), imx, imy, imz]
    e = exp(τ*log(q₂ / q₁))
    ∂e = ∂exp(τ*log(q₂ / q₁))
    ∂l = ∂log(q₂ / q₁)
    [
        sum(
            ∂e[b] * τ * ∂l[c][b] * (
                q₂ * conj(basis[a]) #/ abs2(q₁)
                - 2q₁[a] * (q₂/q₁) #/ abs2(q₁)^2
            )[c] * q₁
            for b in 1:4
            for c in 1:4
        ) + e * basis[a]
        for a in 1:4
    ]
end


function ∂slerp∂q₂(q₁, q₂, τ)
    basis = Rotor[Quaternion(1), imx, imy, imz]
    e = exp(τ*log(q₂ / q₁))
    ∂e = ∂exp(τ*log(q₂ / q₁))
    ∂l = ∂log(q₂ / q₁)
    [
        sum(
            ∂e[b] * τ * ∂l[c][b] * (basis[a] / q₁)[c] * q₁
            for b in 1:4
            for c in 1:4
        )
        for a in 1:4
    ]
end


∂slerp∂τ(q₁, q₂, τ) = log(q₂/q₁) * slerp(q₁, q₂, τ)


function ∂squad∂t(qᵢ, A, B, qᵢ₊₁, ta, tb, t)
    # squad = slerp(slerp(qᵢ, qᵢ₊₁, τ), slerp(A, B, τ), 2τ*(1-τ))
    # τ = (t - ta) / (tb - ta)
    # ∂τ∂t = 1 / (tb - ta)
    # ∂(2τ*(1-τ))/∂t = ∂τ/∂t * (2 - 4τ)

    τ = (t - ta) / (tb - ta)
    ∂τ∂t = 1 / (tb - ta)

    X = slerp(qᵢ, qᵢ₊₁, τ)
    ∂X∂t = ∂τ∂t * ∂slerp∂τ(qᵢ, qᵢ₊₁, τ)

    Y = slerp(A, B, τ)
    ∂Y∂t = ∂τ∂t * ∂slerp∂τ(A, B, τ)

#     ∂X∂qᵢ = ∂slerp∂q₁(qᵢ, qᵢ₊₁, τ)
#     ∂X∂qᵢ₊₁ = ∂slerp∂q₂(qᵢ, qᵢ₊₁, τ)
#     ∂Y∂A = ∂slerp∂q₁(A, B, τ)
#     ∂Y∂B = ∂slerp∂q₂(A, B, τ)

    ∂s∂X = ∂slerp∂q₁(X, Y, 2τ*(1-τ))
    ∂s∂Y = ∂slerp∂q₂(X, Y, 2τ*(1-τ))
    ∂s∂τ = ∂slerp∂τ(X, Y, 2τ*(1-τ))

    ∂s∂t = (
        sum(∂X∂t[b] * ∂s∂X[b] for b in 1:4)
        + sum(∂Y∂t[b] * ∂s∂Y[b] for b in 1:4)
        + ∂τ∂t * (2 - 4τ) * ∂s∂τ
    )
end
;
