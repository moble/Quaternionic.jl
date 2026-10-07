# The pushforwards and pullbacks of `exp` and `log` of quaternions, and the power series
# that keep them accurate near the identity.  The time derivative that `squad` computes is
# built on the pushforwards, and the ChainRulesCore extension uses all of them in its rules
# for `exp`, `log`, and powers.  The pullbacks of `log` are evaluated in forms that cannot
# overflow or underflow where the primal does not.
#
# The functions are generic in the element type, so that nested AD (such as ForwardDiff
# running through a ChainRules rule) works, and every branch inside them tests values,
# through `value` and `iszerovalue`.  Near the identity, the functions of the vector part
# are evaluated as power series, so that all of their derivatives stay accurate there.  The
# number of terms depends on the precision of the element type, so that the series are as
# accurate for `BigFloat` as for `Float64`.

########################################################################################
## Small helpers
########################################################################################

"""
    seriestype(x)

Return the real floating-point type of the value of `x`, whose precision decides how many
terms of a power series are needed.
"""
seriestype(x) = real(float(typeof(value(x))))

"""
    bilinearvec(a, b)

Return the bilinear product `Σᵢ aᵢ bᵢ` of the vector components (2 to 4) of the quaternions
`a` and `b`, without complex conjugation.
"""
bilinearvec(a, b) = a[2] * b[2] + a[3] * b[3] + a[4] * b[4]

########################################################################################
## Power series
########################################################################################

"""
    expseriesterms(T)

Return the number of terms of the series in `cosroot`, `sincroot`, and `sincrootslope` that
makes their truncation error smaller than the precision of `T` for every `|x| < 1`.  This is
the smallest `n` with `(2n)! ≥ 2^(precision(T)+3)`, which is 10 for `Float64`.  For types
without a precision, it is 10.
"""
function expseriesterms(::Type{T}) where {T<:AbstractFloat}
    bits = precision(T) + 3
    n, l = 1, 1.0  # l = log2((2n)!)
    while l < bits
        n += 1
        l += log2((2n - 1) * (2n))
    end
    return n
end
# Julia replaces a method that returns a constant with that constant, so coverage never
# counts the lines of these methods, or of the last method of `logseriesterms`, below.
expseriesterms(::Type) = 10  # COV_EXCL_LINE
# The values for the IEEE types, which the method above also returns, are given explicitly,
# so that they cost nothing at run time.
expseriesterms(::Type{Float64}) = 10  # COV_EXCL_LINE
expseriesterms(::Type{Float32}) = 6  # COV_EXCL_LINE
expseriesterms(::Type{Float16}) = 4  # COV_EXCL_LINE

"""
    logseriesterms(T)

Return the number of terms of the series in `atanroot` and `atanrootslope` that makes their
truncation error smaller than the precision of `T` for every `0 ≤ x ≤ 1/16`.  This is 14 for
`Float64`.  For types without a precision, it is 24.
"""
logseriesterms(::Type{T}) where {T<:AbstractFloat} = cld(precision(T) + 3, 4)
logseriesterms(::Type) = 24  # COV_EXCL_LINE

"""
    rootclosedforms(x)

Return `(cos(√x), sin(√x)/√x)` in closed form.  Both are even functions of √x, so the branch
of the square root does not matter.  When the real part of `x` is negative, they are
evaluated as `cosh(√(-x))` and `sinh(√(-x))/√(-x)`, so that the square root is never taken
on or near its branch cut.  (On the cut, ForwardDiff's derivative of the complex square
root of a number whose imaginary part is a negative zero belongs to the opposite branch from
its value, which spoils nested derivatives.)
"""
function rootclosedforms(x)
    if real(value(x)) < 0
        b = sqrt(-x)
        return cosh(b), sinh(b) / b
    end
    a = sqrt(x)
    return cos(a), sin(a) / a
end

"""
    cosroot(x)

Return cos(√x) = Σₖ (-x)ᵏ/(2k)!, evaluated as a power series when `abs(value(x)) < 1`.
"""
function cosroot(x)
    abs(value(x)) < 1 || return first(rootclosedforms(x))
    r = one(x)
    for k ∈ expseriesterms(seriestype(x))-1:-1:1
        r = 1 - x / ((2k - 1) * (2k)) * r
    end
    return r
end

"""
    sincroot(x)

Return sin(√x)/√x = Σₖ (-x)ᵏ/(2k+1)!, evaluated as a power series when `abs(value(x)) < 1`.
"""
function sincroot(x)
    abs(value(x)) < 1 || return last(rootclosedforms(x))
    r = one(x)
    for k ∈ expseriesterms(seriestype(x))-1:-1:1
        r = 1 - x / ((2k) * (2k + 1)) * r
    end
    return r
end

"""
    sincrootslope(x)

Return (cos(√x) - sin(√x)/√x)/x, which is twice the derivative of `sincroot(x)`, evaluated
as the power series Σₖ (-1)ᵏ⁺¹ 2(k+1) xᵏ/(2k+3)! when `abs(value(x)) < 1`.
"""
function sincrootslope(x)
    if abs(value(x)) ≥ 1
        C, S = rootclosedforms(x)
        return (C - S) / x
    end
    # The ratio of consecutive terms is -x/(2k(2k+3)), and the first term is -1/3.
    r = one(x)
    for k ∈ expseriesterms(seriestype(x))-1:-1:1
        r = 1 - x / ((2k) * (2k + 3)) * r
    end
    return -r / 3
end

"""
    atanroot(x)

Return atan(√x)/√x = Σₖ (-x)ᵏ/(2k+1) as a power series, accurate for `0 ≤ x ≤ 1/16`.
"""
function atanroot(x)
    # The ratio of consecutive terms is -(2k-1)x/(2k+1), and the first term is 1.
    r = one(x)
    for k ∈ logseriesterms(seriestype(x))-1:-1:1
        r = 1 - x * (2k - 1) / (2k + 1) * r
    end
    return r
end

"""
    atanrootslope(x)

Return (1/(1+x) - atanroot(x))/x = Σₖ₌₁ (-1)ᵏ 2k/(2k+1) xᵏ⁻¹ as a power series, accurate for
`0 ≤ x ≤ 1/16`.
"""
function atanrootslope(x)
    # The ratio of consecutive terms is -x(k+1)(2k+1)/(k(2k+3)), and the first term is -2/3.
    r = one(x)
    for k ∈ logseriesterms(seriestype(x))-1:-1:1
        r = 1 - x * ((k + 1) * (2k + 1)) / (k * (2k + 3)) * r
    end
    return -2r / 3
end

########################################################################################
## exp
########################################################################################

# exp(w + v⃗) = E (C + S v⃗), with E = exp(w), x = Σᵢ vᵢ² (bilinear for complex components),
# C = cos√x, and S = sin√x/√x.  With D = (C - S)/x = 2 dS/dx, the pushforward is
#   Ω̇ = ẇ Ω + E (-S ⟨v⃗, v⃗̇⟩ + S v⃗̇ + D ⟨v⃗, v⃗̇⟩ v⃗),
# with the bilinear ⟨⋅, ⋅⟩.  The pullback applies the conjugate transpose, which evaluates
# the same coefficients at the componentwise conjugate of `q`.  The factor Ω in the first
# term is recomputed from C and S, rather than taken from the primal, so that the rules do
# not depend on how well an AD system differentiates the source's own evaluation of `exp`
# when nested AD runs through them.  (With complex components, the source's `absvec` takes a
# complex square root that may lie on its branch cut.)

"""
    expfactors(q)

Return `(E, C, S, D)` for the quaternion `q = w + v⃗`: E = exp(w), and with x = Σᵢ vᵢ²,
C = cos√x, S = sin√x/√x, and D = (C - S)/x.
"""
function expfactors(q)
    x = q[2] * q[2] + q[3] * q[3] + q[4] * q[4]
    return exp(q[1]), cosroot(x), sincroot(x), sincrootslope(x)
end

"""
    exp_pushforward(q, q̇)

Return the pushforward of `exp` at the `Quaternion` `q`, applied to the `Quaternion` tangent
`q̇`.
"""
function exp_pushforward(q, q̇)
    E, C, S, D = expfactors(q)
    vv̇ = bilinearvec(q, q̇)
    a = E * S
    b = E * D * vv̇ + a * q̇[1]
    return quaternion(
        E * C * q̇[1] - a * vv̇,
        a * q̇[2] + b * q[2], a * q̇[3] + b * q[3], a * q̇[4] + b * q[4]
    )
end

"""
    exp_pullback(q, Δ)

Return the pullback of `exp` at the `Quaternion` `q`, applied to the `Quaternion` cotangent
`Δ`.
"""
function exp_pullback(q, Δ)
    p = quaternion(conj(q[1]), conj(q[2]), conj(q[3]), conj(q[4]))
    E, C, S, D = expfactors(p)
    vΔ = bilinearvec(p, Δ)
    a = E * S
    b = E * D * vΔ - a * Δ[1]
    return quaternion(
        E * C * Δ[1] + a * vΔ, a * Δ[2] + b * p[2], a * Δ[3] + b * p[3], a * Δ[4] + b * p[4]
    )
end

########################################################################################
## log (real components)
########################################################################################

# log(q) = log n + f v⃗ for q = w + v⃗, with n = |q|, a = |v⃗|, and f = atan(a, w)/a.  With
# g = (w/n² - f)/a², the pushforward is
#   Ω̇ = ⟨q, q̇⟩/n² + f v⃗̇ + v⃗ (g ⟨v⃗, v⃗̇⟩ - ẇ/n²),
# and the pullback is
#   ∂w = (w Δ₁ - ⟨v⃗, Δ⃗⟩)/n²,   ∂v⃗ = (Δ₁/n² + g ⟨v⃗, Δ⃗⟩) v⃗ + f Δ⃗.
# Since log(λq) = log λ + log q for λ > 0, both are evaluated at the unit quaternion p = q/n
# (where n = 1) and divided by n, which cannot overflow or underflow.  Near the positive
# real axis (w > 0 and a² ≤ w²/16), with x = a²/w², f = A(x)/w and g = B(x)/w³, where A is
# `atanroot` and B is `atanrootslope`.  Elsewhere, g ⟨v⃗, ⋅⟩ v⃗ is evaluated as
# (w - f) ⟨û, ⋅⟩ û with the unit vector û = v⃗/a, so that neither a² nor 1/a² is ever
# formed.  Exactly on the negative real axis, the source returns the constant vector part
# π𝐤, so f = g = 0 there.

"""
    logfactors(p)

Return `(f, h, isseries)` for the unit `Quaternion` `p` with real components, with `f` as in
the comment above.  When `isseries` is true, `h` is `g`, which multiplies ⟨v⃗, ⋅⟩ v⃗.
Otherwise, `h` is `w - f`, which multiplies ⟨û, ⋅⟩ û, because g = (w - f)/a² when n = 1.
"""
function logfactors(p)
    w = p[1]
    if value(w) > 0
        x = (p[2] / w)^2 + (p[3] / w)^2 + (p[4] / w)^2
        if value(x) ≤ 1//16
            return atanroot(x) / w, atanrootslope(x) / w^3, true
        end
    end
    if iszerovalue(vec(p))  # On the negative real axis
        return zero(w), zero(w), true
    end
    a = absvec(p)
    f = atan(a, w) / a
    return f, w - f, false
end

"""
    log_pullback(q, Δ)

Return the pullback of `log` at the `Quaternion` `q` with real components, applied to the
`Quaternion` cotangent `Δ`.  At `q = 0`, where the source returns a constant, the result is
zero.
"""
function log_pullback(q, Δ)
    n = abs(q)
    if iszerovalue(q)
        z = zero(Δ[1] / n)
        return quaternion(z, z, z, z)
    end
    p = q / n
    f, h, isseries = logfactors(p)
    vΔ = bilinearvec(p, Δ)
    # The term g ⟨v⃗, Δ⃗⟩ v⃗ is c u⃗, as in `log_pushforward`.
    c, u = if isseries
        h * vΔ, p
    else
        a = absvec(p)
        h * (vΔ / a), p / a
    end
    ∂ = quaternion(
        p[1] * Δ[1] - vΔ,
        Δ[1] * p[2] + c * u[2] + f * Δ[2],
        Δ[1] * p[3] + c * u[3] + f * Δ[3],
        Δ[1] * p[4] + c * u[4] + f * Δ[4]
    )
    return ∂ / n
end

"""
    log_pushforward(q, q̇)

Return the pushforward of `log` at the `Quaternion` `q` with real components, applied to the
`Quaternion` tangent `q̇`.  At `q = 0`, where the source returns a constant, the result is
zero.
"""
function log_pushforward(q, q̇)
    n = abs(q)
    if iszerovalue(q)
        z = zero(q̇[1] / n)
        return quaternion(z, z, z, z)
    end
    p = q / n
    ṗ = q̇ / n
    f, h, isseries = logfactors(p)
    vv̇ = bilinearvec(p, ṗ)
    # The term g ⟨v⃗, v⃗̇⟩ v⃗ is c u⃗.
    c, u = if isseries
        h * vv̇, p
    else
        a = absvec(p)
        h * (vv̇ / a), p / a
    end
    return quaternion(
        p[1] * ṗ[1] + vv̇,
        f * ṗ[2] + c * u[2] - ṗ[1] * p[2],
        f * ṗ[3] + c * u[3] - ṗ[1] * p[3],
        f * ṗ[4] + c * u[4] - ṗ[1] * p[4]
    )
end
