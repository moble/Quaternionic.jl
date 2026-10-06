@doc raw"""
    slerp∂slerp(q₁, q₂, τ)

Return the value and gradient of `slerp`.

The gradient is with respect to each of the input arguments in turn, with each quaternion
regarded as a series of four arguments.  That is, a total of 10 quaternions will be
returned:
```math
\begin{aligned}
  \big[
    &\mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₁.w} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₁.x} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₁.y} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₁.z} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₂.w} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₂.x} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₂.y} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial q₂.z} \mathrm{slerp}(q₁, q₂, τ), \\
    &\frac{\partial}{\partial \tau} \mathrm{slerp}(q₁, q₂, τ)
  \big]
\end{aligned}
```
For convenience, this will be a 4-tuple with the `slerp` as the first element, the first
four components of the derivative, followed by the next four components of the derivative,
followed by the last component of the derivative.  The first element is a `Rotor`, and each
derivative is a general `Quaternion`.  As with [`∂log`](@ref), the derivatives with respect
to the components of `q₁` and `q₂` include "off-shell" directions that would change their
norms.

See also [`slerp`](@ref) for just the value, and [`slerp∂slerp∂τ`](@ref) for just the value
and the derivative with respect to `τ`.

# Examples
```julia
julia> (q₁, q₂), τ = randn(RotorF64, 2), rand();

julia> s, ∂s∂q₁, ∂s∂q₂, ∂s∂τ = slerp∂slerp(q₁, q₂, τ);

```
"""
function slerp∂slerp(q₁::Rotor, q₂::Rotor, τ::Real)
    # slerp(q₁, q₂, τ) = exp(τ log(r)) q₁, where r = q₂ / q₁.
    r = q₂ / q₁
    l, ∂l = log∂log(r)
    e, ∂e = exp∂exp(τ*l)
    s = e * q₁
    # ∂e∂r[c] is the derivative of e with respect to r[c].
    ∂e∂r = [τ * sum(∂e[b] * ∂l[c][b] for b in 1:4) for c in 1:4]
    # The derivative of s along a change δr of r is (Σ_c ∂e∂r[c] δr[c]) q₁.
    ∂s∂r(δr) = sum(∂e∂r[c] * δr[c] for c in 1:4) * q₁
    T = basetype(r)
    basis = (
        Quaternion{T}(1, 0, 0, 0), Quaternion{T}(0, 1, 0, 0),
        Quaternion{T}(0, 0, 1, 0), Quaternion{T}(0, 0, 0, 1)
    )
    (
        s,
        [∂s∂r(q₂ * conj(basis[a]) - 2q₁[a] * r) + e * basis[a] for a in 1:4],
        [∂s∂r(basis[a] / q₁) for a in 1:4],
        l * s
    )
end


"""
    slerp∂slerp∂τ(q₁, q₂, τ)

Return the value of `slerp` and its derivative with respect to `τ`.

See also [`slerp∂slerp`](@ref), which returns the value and *all* of the derivatives of
`slerp`.

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

This is primarily an internal helper function, taking various parameters computed within the
`squad` function.  This will be used to compute the derivative of `squad` when the angular
velocity is also requested.  To actually obtain the derivative, simply pass the relevant
keyword to the `squad` function.

"""
function squad∂squad∂t(qᵢ, A, B, qᵢ₊₁, ta, tb, t)
    # squad = slerp(X, Y, σ), with
    #   X = slerp(qᵢ, qᵢ₊₁, τ),  Y = slerp(A, B, τ),  σ = 2τ(1-τ),  τ = (t - ta) / (tb - ta).
    # Writing slerp(X, Y, σ) = exp(σ log(r)) X with r = Y / X, the derivative along τ is
    #   ∂s/∂τ = (∂e/∂τ) X + e ∂X/∂τ,  where  e = exp(σ log(r)),
    # and ∂e/∂τ follows from the chain rule through `exp∂exp` and `log∂log`, using
    #   ∂r/∂τ = (∂Y/∂τ - r ∂X/∂τ) / X.
    # This directional derivative avoids building the full Jacobians of `slerp∂slerp`.

    τ = (t - ta) / (tb - ta)
    ∂τ∂t = 1 / (tb - ta)

    X, ∂X∂τ = slerp∂slerp∂τ(qᵢ, qᵢ₊₁, τ)
    Y, ∂Y∂τ = slerp∂slerp∂τ(A, B, τ)

    σ = 2τ*(1-τ)
    ∂σ∂τ = 2 - 4τ

    r = Y / X
    ∂r∂τ = (∂Y∂τ - r * ∂X∂τ) / X
    l, ∂l = log∂log(r)
    ∂l∂τ = sum(∂l[c] * ∂r∂τ[c] for c in 1:4)
    e, ∂e = exp∂exp(σ*l)
    ∂σl∂τ = ∂σ∂τ * l + σ * ∂l∂τ
    ∂e∂τ = sum(∂e[b] * ∂σl∂τ[b] for b in 1:4)

    s = e * X
    ∂s∂τ = ∂e∂τ * X + e * ∂X∂τ

    (s, ∂τ∂t * ∂s∂τ)
end
