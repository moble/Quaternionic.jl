@doc raw"""
    ∂log(Z::Rotor)

Return the gradient of `log(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even though change
of `Z` in a direction that changes its norm would not be allowed for a `Rotor`, we measure
the gradient in that direction anyway.  Specifically, the gradient is that of the extension
of the logarithm to general quaternions,
```math
\log(Z) = \log|Z| + \frac{\mathrm{atan}(|\vec{Z}|, Z_w)}{|\vec{Z}|}\, \vec{Z},
```
so that the elements of the returned vector of quaternions are
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

Note that, even though `log(::Rotor)` is a `QuatVec`, the derivative (and therefore each
element of the result) is a general `Quaternion`.

See also [`∂exp`](@ref) for a similar function, as well as [`log∂log`](@ref) for a function
to compute the value along with the gradient.

# Examples
```julia
julia> ∂log∂w, ∂log∂x, ∂log∂y, ∂log∂z = ∂log(randn(RotorF64));

```
"""
∂log(Z::Rotor) = log∂log(Z)[2]


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
    # With w = Z[1] and a = absvec(Z), the logarithm is extended off shell as
    #
    #   log(Z) = log(|Z|) + f vec(Z),  where  f = atan(a, w) / a.
    #
    # Its derivatives are
    #
    #   ∂log(Z)/∂w  = conj(Z) / |Z|²,
    #   ∂log(Z)/∂Zᵢ = Zᵢ (1/|Z|² + f′ vec(Z)) + f 𝐞ᵢ,  where  f′ = (w/|Z|² - f) / a².
    #
    # Near the identity, the difference in f′ cancels catastrophically, and the square root
    # in `a` has an infinite derivative at a = 0.  There, both f and f′ are instead
    # computed from series in x = a²/w², which carry the derivatives of every input
    # component through AD.  The series are f = A(x)/w and f′ = S(x)/w³, where
    #
    #   A(x) = Σₖ (-1)ᵏ xᵏ / (2k+1),  and  S(x) = Σₖ (-1)ᵏ⁺¹ 2(k+1) xᵏ / (2k+3).
    #
    # Below the threshold x ≤ eps^(1/5), the six terms used here truncate each series at
    # well below `eps`.
    #
    # The term Zᵢ f′ vec(Z) is computed as uᵢ k u.  In the series branch, u = vec(Z) and
    # k = f′.  Otherwise, u = vec(Z)/a and k = f′a², because f′ alone overflows near the
    # antipode, where it grows like 1/a³.
    Z = float(Z)
    w = Z[1]
    a² = abs2vec(Z)
    n² = w*w + a²
    if value(w) > 0 && value(a²) ≤ value(w*w) * eps(typeof(value(a²)))^(1//5)
        x = a² / (w*w)
        f = evalpoly(x, (1, -1//3, 1//5, -1//7, 1//9, -1//11)) / w
        k = evalpoly(x, (-2//3, 4//5, -6//7, 8//9, -10//11, 12//13)) / (w*w*w)
        u = quatvec(Z)
        l = f * quatvec(Z)
    elseif iszerovalue(quatvec(Z))
        # Z is a negative real number, where the logarithm is singular, and its derivatives
        # with respect to the vector components are infinite.  As `log` does, we return π𝐤
        # as the value.  The gradient keeps only the finite derivatives of log(|Z|), which
        # is also what AD of `log(::Quaternion)` gives at this point.  The components are
        # tested, rather than `a²`, because `a²` underflows to zero for tiny nonzero vector
        # parts, which belong in the general branch below.
        f = zero(n²)
        k = zero(n²)
        u = quatvec(Z)
        l = QuatVec{typeof(n²)}(false, false, false, π)
    else
        a = absvec(Z)
        f = atan(a, w) / a
        k = w / n² - f
        u = quatvec(Z) / a
        # Multiplying the angle by `u` avoids overflow in `f` for a tiny `a`.
        l = atan(a, w) * u
    end
    ∂l = [
        quaternion(conj(Z)) / n²,
        Z[2] / n² + u[2] * k * u + f * imx,
        Z[3] / n² + u[3] * k * u + f * imy,
        Z[4] / n² + u[4] * k * u + f * imz,
    ]
    (l, ∂l)
end


@doc raw"""
    ∂exp(Z::QuatVec)

Return the gradient of `exp(Z)` with respect to the components of `Z`.

The result includes "off-shell" components of the gradient, meaning that even though a
scalar component of `Z` would not be allowed for a `QuatVec`, we measure the gradient in
that direction anyway.  That is, the first element of the returned vector of quaternions is
```math
\begin{aligned}
  \left.\frac{\partial} {\partial Z_w} \exp(Z) \right|_{Z_w=0}.
\end{aligned}
```

Note that, even though `exp(::QuatVec)` is a `Rotor`, the derivative (and therefore each
element of the result) is a general `Quaternion`.

See also [`∂log`](@ref) for a similar function, as well as [`exp∂exp`](@ref) for a function
to compute the value along with the gradient.

# Examples
```julia
julia> ∂exp∂w, ∂exp∂x, ∂exp∂y, ∂exp∂z = ∂exp(randn(QuatVecF64));

```
"""
∂exp(Z::QuatVec) = exp∂exp(Z)[2]


"""
    exp∂exp(Z::QuatVec)

Return the value and gradient of `exp(Z)` with respect to the components of `Z`.

The value is a `Rotor`, as returned by `exp(Z)`, while the elements of the gradient are
general `Quaternion`s.  See [`∂exp`](@ref) for more explanation of the components of the
gradient.

# Examples
```julia
julia> e, ∂e = exp∂exp(randn(QuatVecF64));

```
"""
function exp∂exp(Z::QuatVec)
    # With a = absvec(Z), the exponential is extended off shell as
    #
    #   exp(Z) = exp(Z_w) (c + g vec(Z)),  where  c = cos(a)  and  g = sin(a) / a.
    #
    # Its derivatives at Z_w = 0 are
    #
    #   ∂exp(Z)/∂Z_w = exp(Z),
    #   ∂exp(Z)/∂Zᵢ  = Zᵢ (-g + g′ vec(Z)) + g 𝐞ᵢ,  where  g′ = (c - g) / a².
    #
    # Near a = 0, the difference in g′ cancels catastrophically, and the square root in `a`
    # has an infinite derivative at a = 0.  There, c, g, and g′ are instead computed from
    # their Taylor series in a², which carry the derivatives of every input component
    # through AD.  The terms of g′ are (-1)ᵏ⁺¹ 2(k+1) a²ᵏ / (2k+3)!.
    Z = float(Z)
    a² = abs2vec(Z)
    if value(a²) ≤ eps(typeof(value(a²)))^(1//5)
        c = evalpoly(a², (1, -1//2, 1//24, -1//720, 1//40320, -1//3628800))
        g = evalpoly(a², (1, -1//6, 1//120, -1//5040, 1//362880, -1//39916800))
        g′ = evalpoly(a², (-1//3, 1//30, -1//840, 1//45360, -1//3991680, 1//518918400))
    else
        a = sqrt(a²)
        # Separate `sin` and `cos` calls (rather than `sincos`) keep second-order Enzyme
        # working.
        c = cos(a)
        g = sin(a) / a
        g′ = (c - g) / a²
    end
    e = Rotor{typeof(c)}(c, g*Z[2], g*Z[3], g*Z[4])
    V = -g + g′ * Z
    ∂e = [
        quaternion(e),
        Z[2] * V + g * imx,
        Z[3] * V + g * imy,
        Z[4] * V + g * imz,
    ]
    (e, ∂e)
end
