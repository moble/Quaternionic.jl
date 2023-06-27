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
