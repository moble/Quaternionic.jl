# Rules for `exp`, `log`, `sqrt`, powers, and `csqrt` (sections 3.3, 3.4, and 3.7 of the AD
# specification).
#
# Each rule is the exact derivative of what the source in `src/math.jl` computes from the
# stored components, including off the unit sphere for a `Rotor`, whose `log`, `sqrt`, and
# powers assume unit norm in ways that the rules reproduce.  The primal is always computed
# by calling the function itself.  The pullbacks and pushforwards are generic in the element
# type, so that nested AD (such as `Zygote.hessian`, which runs ForwardDiff through them)
# works, and every branch inside them tests values, through `value` and `iszerovalue`.
#
# Near the identity, the functions of the vector part are evaluated as power series, so that
# all of their derivatives stay accurate there.  The number of terms depends on the
# precision of the element type, so that the series are as accurate for `BigFloat` as for
# `Float64`.  The pullbacks of `log` and `sqrt` are evaluated in forms that cannot overflow
# or underflow where the primal does not.
#
# Only real components have rules for `log`, `sqrt`, non-integer powers, and integer powers
# of a `Rotor`; with complex components these are opted out in `optouts.jl`, so that the AD
# system differentiates the source.  `exp` and integer powers of a `Quaternion` or `QuatVec`
# have rules for any components.

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

"""
    vecpart(q)

Return the `Quaternion` with the vector part of `q` and a zero scalar part.
"""
vecpart(q) = quaternion(zero(q[2]), q[2], q[3], q[4])

"""
    pushtangent(x, ẋ)

Return the tangent `ẋ` of the argument `x` of a pushforward, normalized and projected onto
`x`.  This is a `Quaternion` for a quaternion `x` and a number for a number `x`, and it is a
zero of that type when `ẋ` is `nothing` or an `AbstractZero`, or when `x` is not
differentiable.
"""
function pushtangent(x::AbstractQuaternion, ẋ)
    z = zero(float(quaternion(x)))
    Δ = cot(ẋ)
    Δ isa AbstractZero && return z
    t = ProjectTo(x)(Δ)
    t isa AbstractZero && return z
    return quaternion(t)
end
function pushtangent(x::Number, ẋ)
    z = zero(float(x))
    Δ = scalarcot(ẋ)
    Δ isa AbstractZero && return z
    t = ProjectTo(x)(Δ)
    t isa AbstractZero && return z
    return t
end

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
expseriesterms(::Type) = 10
# The values for the IEEE types, which the method above also returns, are given explicitly,
# so that they cost nothing at run time.
expseriesterms(::Type{Float64}) = 10
expseriesterms(::Type{Float32}) = 6
expseriesterms(::Type{Float16}) = 4

"""
    logseriesterms(T)

Return the number of terms of the series in `atanroot` and `atanrootslope` that makes their
truncation error smaller than the precision of `T` for every `0 ≤ x ≤ 1/16`.  This is 14 for
`Float64`.  For types without a precision, it is 24.
"""
logseriesterms(::Type{T}) where {T<:AbstractFloat} = cld(precision(T) + 3, 4)
logseriesterms(::Type) = 24

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
# term is recomputed from C and S, rather than taken from the primal, so that the rule does
# not depend on how well an AD system differentiates the source's own evaluation of `exp`
# when nested AD runs through the rule.  (With complex components, the source's `absvec`
# takes a complex square root that may lie on its branch cut.)

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
    p = conjcomponents(q)
    E, C, S, D = expfactors(p)
    vΔ = bilinearvec(p, Δ)
    a = E * S
    b = E * D * vΔ - a * Δ[1]
    return quaternion(
        E * C * Δ[1] + a * vΔ, a * Δ[2] + b * p[2], a * Δ[3] + b * p[3], a * Δ[4] + b * p[4]
    )
end

# `exp(::QuatVec)` is the case w = 0, with an unnormalized `Rotor` result, and
# `exp(::Rotor)` is `exp(quaternion(q))`, so one rule serves every quaternion type.
for f ∈ (:exp, :(Base.FastMath.exp_fast))
    @eval begin
        function rrule(::typeof($f), q::AQ)
            Ω = $f(q)
            p = float(quaternion(q))
            proj = ProjectTo(q)
            function exp_back(ΔΩ)
                Δ = cot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(exp_pullback(p, Δ))
            end
            return Ω, exp_back
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::AQ)
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            ṗ = pushtangent(q, Δ)
            return Ω, ProjectTo(Ω)(exp_pushforward(float(quaternion(q)), ṗ))
        end
    end
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
# (w - f) ⟨û, ⋅⟩ û with the unit vector û = v⃗/a, so that a² is never formed.  Exactly on
# the negative real axis, the source returns the constant vector part π𝐤, so f = g = 0
# there.

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
    α = if isseries
        Δ[1] + h * vΔ
    else
        a = absvec(p)
        Δ[1] + (h * (vΔ / a)) / a
    end
    ∂ = quaternion(
        p[1] * Δ[1] - vΔ, α * p[2] + f * Δ[2], α * p[3] + f * Δ[3], α * p[4] + f * Δ[4]
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
    β = if isseries
        h * vv̇ - ṗ[1]
    else
        a = absvec(p)
        (h * (vv̇ / a)) / a - ṗ[1]
    end
    return quaternion(
        p[1] * ṗ[1] + vv̇, f * ṗ[2] + β * p[2], f * ṗ[3] + β * p[3], f * ṗ[4] + β * p[4]
    )
end

# The source computes log(R) of a real `Rotor` as |R| V for R[1] ≥ 0 (it evaluates θ/sin θ
# with sin θ = a/|R|) and as V for R[1] < 0, where V = vec(log(quaternion(R))) is
# homogeneous of degree 0; exactly on the negative real axis, it returns the constant π𝐤.
# The result is a `QuatVec`, so only the vector part of a cotangent matters.

"""
    rotorlog_pullback(q, L, Δ)

Return the pullback of `log` at the real `Rotor` with components `q` (a `Quaternion`),
where `L` is `log(q)`, applied to the `Quaternion` cotangent `Δ`.
"""
function rotorlog_pullback(q, L, Δ)
    Δv = vecpart(Δ)
    ∂V = log_pullback(q, Δv)
    if value(q[1]) ≥ 0
        r = abs(q)
        return r * ∂V + (bilinearvec(L, Δv) / r) * q
    elseif iszerovalue(vec(q))
        return zero(∂V)
    else
        return ∂V
    end
end

"""
    rotorlog_pushforward(q, q̇)

Return the pushforward of `log` at the real `Rotor` with components `q` (a `Quaternion`),
applied to the `Quaternion` tangent `q̇`, as a `Quaternion` with zero scalar part.
"""
function rotorlog_pushforward(q, q̇)
    V̇ = vecpart(log_pushforward(q, q̇))
    if value(q[1]) ≥ 0
        r = abs(q)
        return r * V̇ + (inner(q, q̇) / r) * vecpart(log(q))
    elseif iszerovalue(vec(q))
        return zero(V̇)
    else
        return V̇
    end
end

for f ∈ (:log, :(Base.FastMath.log_fast))
    @eval begin
        function rrule(::typeof($f), q::Quaternion{<:Real})
            Ω = $f(q)
            p = float(q)
            proj = ProjectTo(q)
            function log_back(ΔΩ)
                Δ = cot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(log_pullback(p, Δ))
            end
            return Ω, log_back
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::Quaternion{<:Real})
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            return Ω, ProjectTo(Ω)(log_pushforward(float(q), pushtangent(q, Δ)))
        end
        function rrule(::typeof($f), R::Rotor{<:Real})
            Ω = $f(R)
            p = float(quaternion(R))
            L = log(p)
            proj = ProjectTo(R)
            function rotorlog_back(ΔΩ)
                Δ = cot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(rotorlog_pullback(p, L, Δ))
            end
            return Ω, rotorlog_back
        end
        function frule((_, Ṙ)::Tuple, ::typeof($f), R::Rotor{<:Real})
            Ω = $f(R)
            Δ = cot(Ṙ)
            Δ isa AbstractZero && return Ω, Δ
            Ω̇ = rotorlog_pushforward(float(quaternion(R)), pushtangent(R, Δ))
            return Ω, ProjectTo(Ω)(Ω̇)
        end
    end
end

########################################################################################
## sqrt (real components)
########################################################################################

# For a `Quaternion` (and a `QuatVec`, through `quaternion`), s = √q = σ + u⃗ satisfies
# s² = q in every branch of the source, including its rescaling branches, which are exact.
# With r = |s|² = |q|, differentiating s² = q gives
#   σ̇ = (σ ẇ + ⟨u⃗, v⃗̇⟩)/(2r),   u⃗̇ = (v⃗̇ - 2σ̇ u⃗)/(2σ),
# and the pullback is
#   ∂w = (σ Δ₁ - ⟨u⃗, Δ⃗⟩)/(2r),   ∂v⃗ = (Δ⃗ + 2 ∂w u⃗)/(2σ).
# The factor 1/r is applied as two divisions by |s|, so that it can neither overflow nor
# underflow.
# On the non-positive real axis, the source returns √(-w) 𝐤, whose only derivative is
# ∂Ω₄/∂w = -1/(2√(-w)).

"""
    sqrt_pullback(q, s, Δ)

Return the pullback of `sqrt` at the `Quaternion` `q` with real components, where `s` is
`sqrt(q)` as a `Quaternion`, applied to the `Quaternion` cotangent `Δ`.
"""
function sqrt_pullback(q, s, Δ)
    if value(q[1]) ≤ 0 && iszerovalue(vec(q))
        ∂w = -Δ[4] / (2s[4])
        z = zero(∂w)
        return quaternion(∂w, z, z, z)
    end
    σ = s[1]
    m = abs(s)
    ∂w = ((σ / m) * Δ[1] - bilinearvec(s, Δ) / m) / (2m)
    c = 2∂w
    d = 2σ
    return quaternion(
        ∂w, (Δ[2] + c * s[2]) / d, (Δ[3] + c * s[3]) / d, (Δ[4] + c * s[4]) / d
    )
end

"""
    sqrt_pushforward(q, s, q̇)

Return the pushforward of `sqrt` at the `Quaternion` `q` with real components, where `s` is
`sqrt(q)` as a `Quaternion`, applied to the `Quaternion` tangent `q̇`.
"""
function sqrt_pushforward(q, s, q̇)
    if value(q[1]) ≤ 0 && iszerovalue(vec(q))
        ṡ₄ = -q̇[1] / (2s[4])
        z = zero(ṡ₄)
        return quaternion(z, z, z, ṡ₄)
    end
    σ = s[1]
    m = abs(s)
    σ̇ = ((σ / m) * q̇[1] + bilinearvec(s, q̇) / m) / (2m)
    c = 2σ̇
    d = 2σ
    return quaternion(
        σ̇, (q̇[2] - c * s[2]) / d, (q̇[3] - c * s[3]) / d, (q̇[4] - c * s[4]) / d
    )
end

# The source computes sqrt(R) of a real `Rotor` with abs(R) ≡ 1, so it is not √ off the unit
# sphere.  For w ≥ 0, with c = √(2(1+w)), Ω = (c/2, v⃗/c).  For w < 0, with s = |v⃗|,
# u⃗ = v⃗/s, and d = √(2(1-w)), Ω = (s/d, u⃗ d/2).  On the non-positive real axis,
# Ω = √(-w) 𝐤.

"""
    rotorsqrt_pullback(q, Δ)

Return the pullback of `sqrt` at the real `Rotor` with components `q` (a `Quaternion`),
applied to the `Quaternion` cotangent `Δ`.
"""
function rotorsqrt_pullback(q, Δ)
    w = q[1]
    if value(w) ≤ 0 && iszerovalue(vec(q))
        ∂w = -Δ[4] / (2sqrt(-w))
        z = zero(∂w)
        return quaternion(∂w, z, z, z)
    elseif value(w) ≥ 0
        c = sqrt(2(1 + w))
        return quaternion(
            Δ[1] / (2c) - bilinearvec(q, Δ) / c^3, Δ[2] / c, Δ[3] / c, Δ[4] / c
        )
    else
        s = absvec(q)
        d = sqrt(2(1 - w))
        h = d / 2
        uΔ = bilinearvec(q, Δ) / s
        ∂s = Δ[1] / d
        ∂d = -Δ[1] * s / d^2 + uΔ / 2
        # With u⃗ = v⃗/s and s = |v⃗|, ∂v⃗ = (h Δ⃗ - h ⟨u⃗, Δ⃗⟩ u⃗)/s + ∂s u⃗,
        # and ∂w = -∂d/d.
        a = h / s
        b = (∂s - a * uΔ) / s
        return quaternion(
            -∂d / d, a * Δ[2] + b * q[2], a * Δ[3] + b * q[3], a * Δ[4] + b * q[4]
        )
    end
end

"""
    rotorsqrt_pushforward(q, q̇)

Return the pushforward of `sqrt` at the real `Rotor` with components `q` (a `Quaternion`),
applied to the `Quaternion` tangent `q̇`.
"""
function rotorsqrt_pushforward(q, q̇)
    w = q[1]
    if value(w) ≤ 0 && iszerovalue(vec(q))
        ṡ₄ = -q̇[1] / (2sqrt(-w))
        z = zero(ṡ₄)
        return quaternion(z, z, z, ṡ₄)
    elseif value(w) ≥ 0
        c = sqrt(2(1 + w))
        b = q̇[1] / c^3
        return quaternion(
            q̇[1] / (2c), q̇[2] / c - b * q[2], q̇[3] / c - b * q[3], q̇[4] / c - b * q[4]
        )
    else
        s = absvec(q)
        d = sqrt(2(1 - w))
        h = d / 2
        ṡ = bilinearvec(q, q̇) / s
        ḋ = -q̇[1] / d
        # Ω⃗ = u⃗ h, with u⃗̇ = (v⃗̇ - u⃗ ṡ)/s and ḣ = ḋ/2
        a = h / s
        b = (ḋ / 2 - a * ṡ) / s
        return quaternion(
            ṡ / d - s * ḋ / d^2,
            a * q̇[2] + b * q[2], a * q̇[3] + b * q[3], a * q̇[4] + b * q[4]
        )
    end
end

# A `QuatVec` has the same rule as a `Quaternion`, through `quaternion`.
const RealQV = Union{Quaternion{<:Real},QuatVec{<:Real}}
for f ∈ (:sqrt, :(Base.FastMath.sqrt_fast))
    @eval begin
        function rrule(::typeof($f), q::RealQV)
            Ω = $f(q)
            p = float(quaternion(q))
            proj = ProjectTo(q)
            function sqrt_back(ΔΩ)
                Δ = cot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(sqrt_pullback(p, Ω, Δ))
            end
            return Ω, sqrt_back
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::RealQV)
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            ṗ = pushtangent(q, Δ)
            return Ω, ProjectTo(Ω)(sqrt_pushforward(float(quaternion(q)), Ω, ṗ))
        end
        function rrule(::typeof($f), R::Rotor{<:Real})
            Ω = $f(R)
            p = float(quaternion(R))
            proj = ProjectTo(R)
            function rotorsqrt_back(ΔΩ)
                Δ = cot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(rotorsqrt_pullback(p, Δ))
            end
            return Ω, rotorsqrt_back
        end
        function frule((_, Ṙ)::Tuple, ::typeof($f), R::Rotor{<:Real})
            Ω = $f(R)
            Δ = cot(Ṙ)
            Δ isa AbstractZero && return Ω, Δ
            ṗ = pushtangent(R, Δ)
            return Ω, ProjectTo(Ω)(rotorsqrt_pushforward(float(quaternion(R)), ṗ))
        end
    end
end

########################################################################################
## Integer powers
########################################################################################

# With Ω = qⁿ for n ≥ 1, Ω̇ = Σₖ qᵏ⁻¹ q̇ qⁿ⁻ᵏ, which is exact everywhere.  The pullback is
# accumulated by the recurrence S₁ = Δ, Sₘ₊₁ = q† Sₘ + Δ (q†)ᵐ, and the pushforward by the
# same recurrence with q in place of q†; neither mutates anything, so that nested AD works.
# A negative power is the inverse of the positive power, as in the source, and the pullback
# of Ω = inv(P) is -Ω† Δ Ω†.

"""
    intpow_pullback(p, n, Ω, Δ)

Return the pullback of `p ↦ pⁿ` at the `Quaternion` `p` (with floating-point components),
where `Ω = pⁿ` as a `Quaternion`, applied to the `Quaternion` cotangent `Δ`.
"""
function intpow_pullback(p, n::Integer, Ω, Δ)
    n == 0 && return zero(p * Δ)
    if n < 0
        Ωa = qadjoint(Ω)
        return intpow_pullback(p, -n, Ω, -(Ωa * Δ * Ωa))
    end
    a = qadjoint(p)
    S = Δ * one(a)
    P = a
    for _ ∈ 2:n
        S = a * S + Δ * P
        P = P * a
    end
    return S
end

"""
    intpow_pushforward(p, n, Ω, ṗ)

Return the pushforward of `p ↦ pⁿ` at the `Quaternion` `p` (with floating-point components),
where `Ω = pⁿ` as a `Quaternion`, applied to the `Quaternion` tangent `ṗ`.
"""
function intpow_pushforward(p, n::Integer, Ω, ṗ)
    n == 0 && return zero(p * ṗ)
    if n < 0
        return -(Ω * intpow_pushforward(p, -n, Ω, ṗ) * Ω)
    end
    S = ṗ * one(p)
    P = p
    for _ ∈ 2:n
        S = p * S + ṗ * P
        P = P * p
    end
    return S
end

"""
    intpow_rule(q, n)

Return the primal and the pullback of `q^n` for an integer `n` and a quaternion `q` that is
not a `Rotor`.  A `QuatVec` is raised to the power as `quaternion(q)`.
"""
function intpow_rule(q::AQ, n::Integer)
    Ω = q^n
    p = float(quaternion(q))
    Ωq = quaternion(Ω)
    proj = ProjectTo(q)
    function intpow_back(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, NoTangent())
        return NoTangent(), proj(intpow_pullback(p, n, Ωq, Δ)), NoTangent()
    end
    return Ω, intpow_back
end

"""
    intpow_frule(q̇, q, n)

Return the primal and the tangent of `q^n` for an integer `n` and a quaternion `q` that is
not a `Rotor`.
"""
function intpow_frule(q̇, q::AQ, n::Integer)
    Ω = q^n
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    ṗ = pushtangent(q, Δ)
    return Ω, ProjectTo(Ω)(intpow_pushforward(float(quaternion(q)), n, quaternion(Ω), ṗ))
end

# The source raises a real `Rotor` to an integer power by repeated normalized products,
# except that R⁰ = 1, R¹ = R, and R⁻¹ = conj(R) are not normalized.  Because |ab| = |a||b|,
# the result for |n| ≥ 2 is normalize(R^|n|), conjugated when n < 0.

"""
    rotorintpow_rule(R, n)

Return the primal and the pullback of `R^n` for a real `Rotor` `R` and an integer `n`.
"""
function rotorintpow_rule(R::Rotor{<:Real}, n::Integer)
    Ω = R^n
    p = float(quaternion(R))
    m = abs(n)
    P = p^m
    proj = ProjectTo(R)
    function rotorintpow_back(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, NoTangent())
        Δ′ = n < 0 ? conj(Δ) : Δ
        ∂ = if m == 0
            zero(p * Δ)
        elseif m == 1
            Δ′ * one(p)
        else
            intpow_pullback(p, m, P, normalize_pullback(P, Δ′))
        end
        return NoTangent(), proj(∂), NoTangent()
    end
    return Ω, rotorintpow_back
end

"""
    rotorintpow_frule(Ṙ, R, n)

Return the primal and the tangent of `R^n` for a real `Rotor` `R` and an integer `n`.
"""
function rotorintpow_frule(Ṙ, R::Rotor{<:Real}, n::Integer)
    Ω = R^n
    Δ = cot(Ṙ)
    Δ isa AbstractZero && return Ω, Δ
    p = float(quaternion(R))
    ṗ = pushtangent(R, Δ)
    m = abs(n)
    Ω̇ = if m == 0
        zero(p * ṗ)
    elseif m == 1
        ṗ * one(p)
    else
        P = p^m
        normalize_pushforward(P, intpow_pushforward(p, m, P, ṗ))
    end
    return Ω, ProjectTo(Ω)(n < 0 ? conj(Ω̇) : Ω̇)
end

########################################################################################
## Real and quaternionic exponents
########################################################################################

# For a base with real components that is not a `Rotor` raised to a real power, the source
# computes q^s = exp(s L) with L = log(p) and p = quaternion(q), multiplying s L in that
# order when s is a quaternion.  A real base x is promoted to quaternion(x).  The pullback
# runs through the pullbacks of exp at s L, of the product s L, and of log at p.

"""
    pow_pullback(p, s, L, Δ)

Return the cotangents `(∂p, ∂s)` of the base `p` (a `Quaternion` with real components) and
the exponent `s` (a real number or a quaternion) of `p^s = exp(s L)`, where `L = log(p)`,
given the `Quaternion` cotangent `Δ`.
"""
function pow_pullback(p, s, L, Δ)
    ∂sL = exp_pullback(q_(s * L), Δ)
    ∂s = s isa AQ ? ∂sL * qadjoint(L) : inner(L, ∂sL)
    ∂L = qadjoint(q_(s)) * ∂sL
    return log_pullback(p, ∂L), ∂s
end

"""
    pow_pushforward(p, s, L, ṗ, ṡ)

Return the tangent of `p^s = exp(s L)`, where `L = log(p)`, given the `Quaternion` tangent
`ṗ` of the base and the tangent `ṡ` of the exponent.
"""
function pow_pushforward(p, s, L, ṗ, ṡ)
    L̇ = log_pushforward(p, ṗ)
    return exp_pushforward(q_(s * L), ṡ * L + s * L̇)
end

"""
    realpow_rule(q, s)

Return the primal and the pullback of `q^s` for a quaternion base `q` with real components
and a real or quaternionic exponent `s`, or for a real base `q` and a quaternionic `s`.
"""
function realpow_rule(q::Union{Real,RealQ}, s::Union{Real,RealQ})
    Ω = q^s
    p = quaternion(float(q))
    L = log(p)
    pq, ps = ProjectTo(q), ProjectTo(s)
    function pow_back(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        ∂p, ∂s = pow_pullback(p, s, L, Δ)
        return NoTangent(), pq(∂p), ps(∂s)
    end
    return Ω, pow_back
end

"""
    realpow_frule(q̇, ṡ, q, s)

Return the primal and the tangent of `q^s` for the arguments of `realpow_rule`.
"""
function realpow_frule(q̇, ṡ, q::Union{Real,RealQ}, s::Union{Real,RealQ})
    Ω = q^s
    p = quaternion(float(q))
    ṗ = quaternion(pushtangent(q, q̇))
    ṡ′ = pushtangent(s, ṡ)
    return Ω, ProjectTo(Ω)(pow_pushforward(p, s, log(p), ṗ, ṡ′))
end

# The source computes R^s for a real `Rotor` R and a real s as exp(s V), where V is the
# vector part of log(quaternion(R)), which is homogeneous of degree 0: both its series
# branch near the identity and its general branch agree with this.  Exactly on the negative
# real axis (a zero vector part with R[1] ≤ 0), it returns cos(πs) + sin(πs) 𝐤, which does
# not depend on R and whose derivative with respect to s is π𝐤 Ω.

"""
    rotorpowlog(p)

Return `(V, onaxis)` for the real `Rotor` with components `p` (a `Quaternion`), where `V` is
the vector part of `log(p)` as a `Quaternion`, or π𝐤 when `onaxis` is true, that is, when
`p` lies on the non-positive real axis.
"""
function rotorpowlog(p)
    if !(value(p[1]) > 0) && iszerovalue(vec(p))
        z = zero(p[1])
        return quaternion(z, z, z, π * one(z)), true
    end
    return vecpart(log(p)), false
end

"""
    rotorpow_rule(R, s)

Return the primal and the pullback of `R^s` for a real `Rotor` `R` and a real number `s`.
"""
function rotorpow_rule(R::Rotor{<:Real}, s::Real)
    Ω = R^s
    p = float(quaternion(R))
    V, onaxis = rotorpowlog(p)
    pR, ps = ProjectTo(R), ProjectTo(s)
    function rotorpow_back(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        if onaxis
            ∂s = inner(V * Ω, Δ)
            return NoTangent(), pR(zero(p * ∂s)), ps(∂s)
        end
        ∂sV = vecpart(exp_pullback(s * V, Δ))
        ∂s = bilinearvec(V, ∂sV)
        return NoTangent(), pR(log_pullback(p, s * ∂sV)), ps(∂s)
    end
    return Ω, rotorpow_back
end

"""
    rotorpow_frule(Ṙ, ṡ, R, s)

Return the primal and the tangent of `R^s` for a real `Rotor` `R` and a real number `s`.
"""
function rotorpow_frule(Ṙ, ṡ, R::Rotor{<:Real}, s::Real)
    Ω = R^s
    p = float(quaternion(R))
    V, onaxis = rotorpowlog(p)
    ṡ′ = pushtangent(s, ṡ)
    Ω̇ = if onaxis
        ṡ′ * (V * quaternion(Ω))
    else
        V̇ = vecpart(log_pushforward(p, pushtangent(R, Ṙ)))
        exp_pushforward(s * V, ṡ′ * V + s * V̇)
    end
    return Ω, ProjectTo(Ω)(Ω̇)
end

########################################################################################
## Method table for `^` and `pow_fast`
########################################################################################

# The signatures are chosen so that no two of them, and none of them and the opt-outs in
# `optouts.jl` (`(AQ, Number)`, `(Number, AQ)`, `(AQ, AQ)`, and
# `(Rotor{<:Complex}, Integer)`), are ambiguous, and so that none of them coincides with an
# opt-out, which would delete it.
for f ∈ (:^, :(Base.FastMath.pow_fast))
    @eval begin
        rrule(::typeof($f), q::AQ, n::Integer) = intpow_rule(q, n)
        rrule(::typeof($f), q::RealQ, n::Integer) = intpow_rule(q, n)
        rrule(::typeof($f), R::Rotor{<:Real}, n::Integer) = rotorintpow_rule(R, n)
        rrule(::typeof($f), q::RealQ, s::Union{Real,RealQ}) = realpow_rule(q, s)
        rrule(::typeof($f), q::RealQ, s::RealQ) = realpow_rule(q, s)
        rrule(::typeof($f), R::Rotor{<:Real}, s::Real) = rotorpow_rule(R, s)
        rrule(::typeof($f), x::Real, s::RealQ) = realpow_rule(x, s)

        frule((_, q̇, _)::Tuple, ::typeof($f), q::AQ, n::Integer) = intpow_frule(q̇, q, n)
        frule((_, q̇, _)::Tuple, ::typeof($f), q::RealQ, n::Integer) =
            intpow_frule(q̇, q, n)
        frule((_, Ṙ, _)::Tuple, ::typeof($f), R::Rotor{<:Real}, n::Integer) =
            rotorintpow_frule(Ṙ, R, n)
        frule((_, q̇, ṡ)::Tuple, ::typeof($f), q::RealQ, s::Union{Real,RealQ}) =
            realpow_frule(q̇, ṡ, q, s)
        frule((_, q̇, ṡ)::Tuple, ::typeof($f), q::RealQ, s::RealQ) =
            realpow_frule(q̇, ṡ, q, s)
        frule((_, Ṙ, ṡ)::Tuple, ::typeof($f), R::Rotor{<:Real}, s::Real) =
            rotorpow_frule(Ṙ, ṡ, R, s)
        frule((_, ẋ, ṡ)::Tuple, ::typeof($f), x::Real, s::RealQ) = realpow_frule(ẋ, ṡ, x, s)
    end
end

# `literal_pow(^, q, Val(n))` computes `q^n`, so its rules forward to those of `^`.  For a
# complex `Rotor`, the rules of `^` are opted out and return `nothing`, which is then passed
# on.  A consumer that treats a returned `nothing` as the absence of a rule then
# differentiates the source.  Zygote does not do so by itself, because no `no_rrule`
# matches this signature, so it relies on the `rrule(::ZygoteRuleConfig, literal_pow, ...)`
# method in QuaternionicZygoteExt, which, given a `nothing`, differentiates `q^n` through
# `rrule_via_ad`.

"""
    literalpow_rule(rule)

Return the primal and pullback of `literal_pow(^, q, Val(n))`, given the result `rule` of
`rrule(^, q, n)`, or `nothing` if `rule` is `nothing`.
"""
literalpow_rule(::Nothing) = nothing
function literalpow_rule((Ω, back)::Tuple)
    literalpow_back(ΔΩ) = (NoTangent(), NoTangent(), back(ΔΩ)[2], NoTangent())
    return Ω, literalpow_back
end

function rrule(::typeof(Base.literal_pow), ::typeof(^), q::AQ, ::Val{n}) where {n}
    return literalpow_rule(rrule(^, q, n))
end
function frule(
    (_, _, q̇, _)::Tuple, ::typeof(Base.literal_pow), ::typeof(^), q::AQ, ::Val{n}
) where {n}
    return frule((NoTangent(), q̇, NoTangent()), ^, q, n)
end

########################################################################################
## csqrt
########################################################################################

# `Quaternionic.csqrt(z)` is `sqrt(z)`, through which every complex square root in the
# source passes.  Its rule exists for Mooncake, which imports it, because Base's
# `sqrt(::Complex)` reinterprets bits in a way that Mooncake cannot trace.  For real
# arguments the AD systems use their own rules for `sqrt`.
function rrule(::typeof(Quaternionic.csqrt), z::Complex)
    Ω = Quaternionic.csqrt(z)
    proj = ProjectTo(z)
    function csqrt_back(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(Δ / conj(2Ω))
    end
    return Ω, csqrt_back
end
function frule((_, ż)::Tuple, ::typeof(Quaternionic.csqrt), z::Complex)
    Ω = Quaternionic.csqrt(z)
    Δ = scalarcot(ż)
    Δ isa AbstractZero && return Ω, Δ
    return Ω, Δ / (2Ω)
end
