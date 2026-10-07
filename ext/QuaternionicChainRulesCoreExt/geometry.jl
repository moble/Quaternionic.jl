# Rules for the norms (`abs`, `abs2`, `absvec`, `abs2vec`, and `norm`), for the geometric
# functions (`R(v)`, `angle`, `distance2`, `distance`, and `dot`), and for the LAPACK
# eigenvector computation behind `from_rotation_matrix` and `align`.
#
# Each rule is the exact derivative of what the source computes from the stored components,
# including off the unit sphere for `Rotor`s.  In particular, `abs`, `abs2`, and (for real
# components) `norm` of a `Rotor` return the constant 1 without reading the components, so
# their derivatives are zero, whereas `absvec` and `abs2vec` of a `Rotor` read the stored
# vector part.  Every branch inside a rule tests values, as the source does, so that dual
# numbers (as in `Zygote.hessian`) take the same branch as the floats they wrap.  Each rule
# computes a single gradient quaternion that serves both the pullback and the pushforward,
# so that each formula exists once.
#
# `slerp` needs no rule of its own: it is `(q₂/q₁)^τ q₁`, and the composition of the rules
# for the `Rotor` quotient, for `Rotor^s`, and for the `Rotor` product is exact, including
# at `q₁ == q₂`.

########################################################################################
## Norms
########################################################################################

# `abs2` is the sum of the squares of the components (the holomorphic spinor norm `Σᵢ zᵢ²`
# for complex components), so its gradient is `2 conjcomponents(q)` in the conventions of
# `projection.jl`.  For a `QuatVec`, the scalar slot is zero by construction, and the
# projection discards it anyway.
function rrule(::typeof(abs2), q::AQ)
    proj = ProjectTo(q)
    function abs2_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(2 * conjcomponents(q) * Δ)
    end
    return abs2(q), abs2_pullback
end
function frule((_, q̇)::Tuple, ::typeof(abs2), q::AQ)
    Ω = abs2(q)
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    Δq = ProjectTo(q)(Δ)
    Δq isa AbstractZero && return Ω, Δq
    return Ω, 2 * bilinear(q, Δq)
end

function rrule(::typeof(abs2vec), q::AQ)
    proj = ProjectTo(q)
    function abs2vec_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(2 * conjcomponents(quatvec(q)) * Δ)
    end
    return abs2vec(q), abs2vec_pullback
end
function frule((_, q̇)::Tuple, ::typeof(abs2vec), q::AQ)
    Ω = abs2vec(q)
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    Δq = ProjectTo(q)(Δ)
    Δq isa AbstractZero && return Ω, Δq
    return Ω, 2 * bilinear(quatvec(q), Δq)
end

"""
    hypotfactor(Δ, Ω)

Return `Δ / conj(Ω)`, the factor by which the conjugated components are multiplied in the
pullback of a norm `Ω` computed by `hypot`, or zero when `Ω` is zero by value.  At that
point the norm has a cone-shaped kink, and the zero is a subgradient, as for `hypot`.
"""
function hypotfactor(Δ, Ω)
    s = Δ / conj(Ω)
    return iszerovalue(Ω) ? zero(s) : s
end

# `abs` and `absvec` are computed by `hypot` (or by its complex analogue, the principal
# square root of the spinor norm), so their gradients are the conjugated components divided
# by the conjugated norm, with a zero cotangent where the norm vanishes.  The real `norm` of
# a quaternion with real components is `abs`.
for (f, part) ∈ ((:abs, :identity), (:absvec, :quatvec), (:norm, :identity))
    @eval begin
        function rrule(::typeof($f), q::AQ)
            Ω = $f(q)
            proj = ProjectTo(q)
            function norm_pullback(ΔΩ)
                Δ = scalarcot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(conjcomponents($part(q)) * hypotfactor(Δ, Ω))
            end
            return Ω, norm_pullback
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::AQ)
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            Δq = ProjectTo(q)(Δ)
            Δq isa AbstractZero && return Ω, Δq
            s = bilinear($part(q), Δq) / Ω
            return Ω, iszerovalue(Ω) ? zero(s) : s
        end
    end
end

# For complex components, `norm` is instead the real Euclidean norm of the eight real
# numbers that make up the components.  It is a real function of complex variables, whose
# ChainRules gradient is `q / Ω` (not conjugated), applied to the real part of the
# cotangent, as for the `norm` of a complex array.  This method also covers a complex
# `Rotor`, whose `norm` reads the stored components.
function rrule(::typeof(norm), q::AbstractQuaternion{<:Complex})
    Ω = norm(q)
    proj = ProjectTo(q)
    function complexnorm_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(quaternion(q) * hypotfactor(real(Δ), Ω))
    end
    return Ω, complexnorm_pullback
end
function frule((_, q̇)::Tuple, ::typeof(norm), q::AbstractQuaternion{<:Complex})
    Ω = norm(q)
    Δ = cot(q̇)
    Δ isa AbstractZero && return Ω, Δ
    Δq = ProjectTo(q)(Δ)
    Δq isa AbstractZero && return Ω, Δq
    s = real(inner(q, Δq)) / Ω
    return Ω, iszerovalue(Ω) ? zero(s) : s
end

# `abs` and `abs2` of any `Rotor`, and `norm` of a `Rotor` with real components, return the
# constant 1 without reading the stored components, so their derivatives are zero.
for (f, RT) ∈ ((:abs, :Rotor), (:abs2, :Rotor), (:norm, :(Rotor{<:Real})))
    @eval begin
        function rrule(::typeof($f), q::$RT)
            rotornorm_pullback(ΔΩ) = (NoTangent(), ZeroTangent())
            return $f(q), rotornorm_pullback
        end
        frule(::Tuple, ::typeof($f), q::$RT) = $f(q), ZeroTangent()
    end
end

########################################################################################
## The action of a rotor on a vector
########################################################################################

# `R(v)` evaluates the vector part of the sandwich product `R v R̄` by a homogeneous
# quadratic formula, valid for any stored `R`, whatever its norm, and returns it as a
# `QuatVec`.  With `Δ⃗` the vector part of the cotangent, the adjoint of `x ↦ a x b` is
# `x ↦ a† x b†`, and quaternion conjugation has a real Jacobian, so
# `∂R = Δ⃗ (v R̄)† + conj((R v)† Δ⃗)` and `∂v = R† Δ⃗ (R̄)†`, where `(R̄)† = conjcomponents(R)`.
function rrule(R::Rotor, v::QuatVec)
    Ω = R(v)
    projR = ProjectTo(R)
    projv = ProjectTo(v)
    Rq = quaternion(R)
    vq = quaternion(v)
    function action_pullback(ΔΩ)
        Δ = cot(ΔΩ)
        Δ isa AbstractZero && return (Δ, Δ)
        Δv = quaternion(false, Δ[2], Δ[3], Δ[4])
        ∂R = Δv * qadjoint(vq * conj(Rq)) + conj(qadjoint(Rq * vq) * Δv)
        ∂v = qadjoint(Rq) * Δv * conjcomponents(Rq)
        return projR(∂R), projv(∂v)
    end
    return Ω, action_pullback
end
function frule((Ṙ, v̇)::Tuple, R::Rotor, v::QuatVec)
    Ω = R(v)
    ΔR = cot(Ṙ)
    Δv = cot(v̇)
    ΔR isa AbstractZero && Δv isa AbstractZero && return Ω, ZeroTangent()
    Rq = quaternion(R)
    vq = quaternion(v)
    ΔRq = ΔR isa AbstractZero ? zero(Rq) : ProjectTo(R)(ΔR)
    Δvq = Δv isa AbstractZero ? zero(vq) : quaternion(false, Δv[2], Δv[3], Δv[4])
    Ω̇ = ΔRq * vq * conj(Rq) + Rq * Δvq * conj(Rq) + Rq * vq * conj(ΔRq)
    return Ω, ProjectTo(Ω)(Ω̇)
end

########################################################################################
## Angle
########################################################################################

"""
    angle_gradient(q)

Return the gradient, as a `Quaternion`, of `angle(q) = 2atan(a, w)` with respect to the
components of the real quaternion `q`, where `w = q[1]` and `a = absvec(q)`.  This is
`(2/n²) (-a, (w/a) v⃗)`, with `n = abs(quaternion(q))` and `v⃗` the vector part of `q`.
Where the vector part vanishes (by value), the source returns a constant, so the gradient
is zero there.
"""
function angle_gradient(q::AbstractQuaternion{<:Real})
    w = q[1]
    zerovec = iszerovalue(vec(q))
    a₀ = absvec(q)
    a = zerovec ? one(a₀) : a₀
    n = abs(quaternion(q))
    c = 2 / n
    wa = (w / n) * c / a
    g = quaternion(-(a / n) * c, wa * q[2], wa * q[3], wa * q[4])
    return zerovec ? zero(g) : g
end

for f ∈ (:angle, :(Base.FastMath.angle_fast))
    @eval begin
        function rrule(::typeof($f), q::Union{Quaternion{<:Real},Rotor{<:Real}})
            proj = ProjectTo(q)
            function angle_pullback(ΔΩ)
                Δ = scalarcot(ΔΩ)
                Δ isa AbstractZero && return (NoTangent(), Δ)
                return NoTangent(), proj(angle_gradient(q) * Δ)
            end
            return $f(q), angle_pullback
        end
        function frule((_, q̇)::Tuple, ::typeof($f), q::Union{Quaternion{<:Real},Rotor{<:Real}})
            Ω = $f(q)
            Δ = cot(q̇)
            Δ isa AbstractZero && return Ω, Δ
            Δq = ProjectTo(q)(Δ)
            Δq isa AbstractZero && return Ω, Δq
            return Ω, bilinear(angle_gradient(q), Δq)
        end
    end
end

########################################################################################
## Distances between rotors
########################################################################################

# The source computes `distance2(R₁, R₂)` from the normalized quotient `Q = R₁ / R₂`, which
# is `R₁ R̄₂ / abs(R₁ R̄₂)`, since the source divides by `abs2(::Rotor) ≡ 1`.  The function of
# `Q` that it evaluates is homogeneous of degree zero, as is its choice of branch (up to
# rounding), so it is the same function of the raw product `P = R₁ R̄₂`, and the rules
# differentiate it there, without the normalization.

"""
    distance2_gradient(P)

Return the gradient, as a `Quaternion`, of the function that `distance2(R₁, R₂)` evaluates
on the quotient `R₁ / R₂`, evaluated at the real quaternion `P`, which may be any positive
multiple of that quotient.  The source evaluates `atan(a, abs(w))^2`, with `w = P[1]` and
`a = absvec(P)`, except where `a²` is small compared with `w²` (by value), where it
evaluates the series `x * evalpoly(x, (1, -2//3, 23//45, -44//105))` in `x = a²/w²`.  This
function differentiates whichever of the two the source evaluates.  The series branch
includes the point `a = 0`, where the gradient is zero and the function is smooth.
"""
function distance2_gradient(P::Quaternion{<:Real})
    w = P[1]
    v² = abs2vec(P)
    w² = w^2
    # These are the threshold and the comparison of the source, applied to the values.
    w₀ = value(w)
    T = w₀ isa AbstractFloat ? typeof(w₀) : Float64
    small = value(v²) ≤ (∜eps(T) / 2) * value(w²)
    if small === true
        x = v² / w²
        # The derivative of the series with respect to `x`, and the chain rule through
        # `x = a²/w²`, for which `∂x/∂w = -2x/w` and `∂x/∂v⃗ = 2v⃗/w²`.
        dfdx = evalpoly(x, (1, -4//3, 23//15, -176//105))
        c = 2dfdx / w²
        return quaternion(-2dfdx * x / w, c * P[2], c * P[3], c * P[4])
    else
        a = absvec(P)
        aw = abs(w)
        n = abs(P)
        θ = atan(a, aw)
        c = 2θ / n
        ca = (aw / n) * c / a
        return quaternion(-sign(w) * (a / n) * c, ca * P[2], ca * P[3], ca * P[4])
    end
end

"""
    rotorquotient_pullback(R₁, R₂, ∂P)

Return the cotangents of the quaternions `R₁` and `R₂` in the raw product
`P = R₁ * conj(R₂)` of their stored components, given the cotangent `∂P` of `P`, for real
components.
"""
function rotorquotient_pullback(R₁::AbstractQuaternion, R₂::AbstractQuaternion, ∂P)
    return ∂P * quaternion(R₂), conj(∂P) * quaternion(R₁)
end

"""
    rotorquotient_pushforward(R₁, R₂, Ṙ₁, Ṙ₂)

Return the tangent of the raw product `P = R₁ * conj(R₂)` of the stored components of the
rotors `R₁` and `R₂`, given their tangents `Ṙ₁` and `Ṙ₂`, either of which may be an
`AbstractZero`, or `ZeroTangent()` if both are.
"""
function rotorquotient_pushforward(R₁::Rotor, R₂::Rotor, Ṙ₁, Ṙ₂)
    Δ₁ = cot(Ṙ₁)
    Δ₂ = cot(Ṙ₂)
    Δ₁ isa AbstractZero && Δ₂ isa AbstractZero && return ZeroTangent()
    q₁ = quaternion(R₁)
    q₂ = quaternion(R₂)
    Δ₁q = Δ₁ isa AbstractZero ? zero(q₁) : ProjectTo(R₁)(Δ₁)
    Δ₂q = Δ₂ isa AbstractZero ? zero(q₂) : ProjectTo(R₂)(Δ₂)
    return Δ₁q * conj(q₂) + q₁ * conj(Δ₂q)
end

function rrule(::typeof(distance2), R₁::Rotor{<:Real}, R₂::Rotor{<:Real})
    Ω = distance2(R₁, R₂)
    p₁, p₂ = ProjectTo(R₁), ProjectTo(R₂)
    function distance2_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        P = quaternion(R₁) * conj(quaternion(R₂))
        ∂₁, ∂₂ = rotorquotient_pullback(R₁, R₂, distance2_gradient(P) * real(Δ))
        return NoTangent(), p₁(∂₁), p₂(∂₂)
    end
    return Ω, distance2_pullback
end
function frule((_, Ṙ₁, Ṙ₂)::Tuple, ::typeof(distance2), R₁::Rotor{<:Real}, R₂::Rotor{<:Real})
    Ω = distance2(R₁, R₂)
    Ṗ = rotorquotient_pushforward(R₁, R₂, Ṙ₁, Ṙ₂)
    Ṗ isa AbstractZero && return Ω, Ṗ
    P = quaternion(R₁) * conj(quaternion(R₂))
    return Ω, distance2_gradient(P) ⋅ Ṗ
end

# `distance(R₁, R₂)` is `√distance2(R₁, R₂)`.  Where it vanishes (by value), as at
# `R₁ == R₂`, it has a cone-shaped kink, and the cotangent is zero, rather than the `0 · ∞`
# that the chain rule through `√` would give.
function rrule(::typeof(distance), R₁::Rotor{<:Real}, R₂::Rotor{<:Real})
    Ω = distance(R₁, R₂)
    p₁, p₂ = ProjectTo(R₁), ProjectTo(R₂)
    function distance_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        P = quaternion(R₁) * conj(quaternion(R₂))
        ∂₁, ∂₂ = rotorquotient_pullback(R₁, R₂, distance2_gradient(P) * hypotfactor(real(Δ), 2Ω))
        return NoTangent(), p₁(∂₁), p₂(∂₂)
    end
    return Ω, distance_pullback
end
function frule((_, Ṙ₁, Ṙ₂)::Tuple, ::typeof(distance), R₁::Rotor{<:Real}, R₂::Rotor{<:Real})
    Ω = distance(R₁, R₂)
    Ṗ = rotorquotient_pushforward(R₁, R₂, Ṙ₁, Ṙ₂)
    Ṗ isa AbstractZero && return Ω, Ṗ
    P = quaternion(R₁) * conj(quaternion(R₂))
    return Ω, hypotfactor(distance2_gradient(P) ⋅ Ṗ, 2Ω)
end

########################################################################################
## The componentwise dot product
########################################################################################

# `p ⋅ q` is the bilinear sum `Σᵢ pᵢ qᵢ` of the components (without complex conjugation), so
# its partial derivatives are `qᵢ` and `pᵢ`, and the pullback multiplies the conjugated
# components by the cotangent.
function rrule(::typeof(dot), p::AQ, q::AQ)
    projp, projq = ProjectTo(p), ProjectTo(q)
    function dot_pullback(ΔΩ)
        Δ = scalarcot(ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        return NoTangent(), projp(conjcomponents(q) * Δ), projq(conjcomponents(p) * Δ)
    end
    return dot(p, q), dot_pullback
end
function frule((_, ṗ, q̇)::Tuple, ::typeof(dot), p::AQ, q::AQ)
    Ω = dot(p, q)
    Δp = cot(ṗ)
    Δq = cot(q̇)
    Δp isa AbstractZero && Δq isa AbstractZero && return Ω, ZeroTangent()
    Δpq = Δp isa AbstractZero ? zero(quaternion(p)) : ProjectTo(p)(Δp)
    Δqq = Δq isa AbstractZero ? zero(quaternion(q)) : ProjectTo(q)(Δq)
    return Ω, bilinear(Δpq, q) + bilinear(p, Δqq)
end

########################################################################################
## The dominant eigenvector of a symmetric matrix, computed by LAPACK
########################################################################################

# `Quaternionic.dominant_eigenvector_lapack(A, uplo)` is the only place where
# `from_rotation_matrix` and `align` call LAPACK, and backends cannot differentiate LAPACK.
# The derivative formulas are plain functions in `src/conversion.jl`, shared with the
# Enzyme extension, and they read only the triangle of `A` (or of its tangent) named by
# `uplo`.  The pullback returns a dense matrix that is zero in the other triangle.
function rrule(
    ::typeof(Quaternionic.dominant_eigenvector_lapack), A::Matrix{<:Union{Float32,Float64}},
    uplo::Char
)
    v = Quaternionic.dominant_eigenvector_lapack(A, uplo)
    proj = ProjectTo(A)
    function dominant_eigenvector_pullback(Δv)
        v̄ = unthunk(Δv)
        v̄ isa Union{Nothing,AbstractZero} && return (NoTangent(), ZeroTangent(), NoTangent())
        Ā = Quaternionic.dominant_eigenvector_pullback(A, uplo, v, collect(v̄))
        return NoTangent(), proj(Ā), NoTangent()
    end
    return v, dominant_eigenvector_pullback
end
function frule(
    (_, ΔA, _)::Tuple, ::typeof(Quaternionic.dominant_eigenvector_lapack),
    A::Matrix{<:Union{Float32,Float64}}, uplo::Char
)
    v = Quaternionic.dominant_eigenvector_lapack(A, uplo)
    Ȧ = unthunk(ΔA)
    Ȧ isa Union{Nothing,AbstractZero} && return v, ZeroTangent()
    return v, Quaternionic.dominant_eigenvector_pushforward(A, uplo, v, collect(Ȧ))
end

# `Quaternionic.value` strips the derivative parts of a (possibly nested) dual number, and
# the source uses its result only in branch decisions, in thresholds, and as the starting
# point of the Newton refinement in `refined_dominant_eigenvector`, whose derivatives the
# Newton steps supply.  Declaring it non-differentiable keeps Zygote from differentiating
# the fields of a dual number, which it cannot do, so that `Zygote.hessian` (ForwardDiff
# over Zygote) works through `from_rotation_matrix` and `align`.  The rules are written out,
# rather than generated by `@non_differentiable`, because the `frule` generated by
# ChainRulesCore 1.15 is ambiguous with its own `frule(::RuleConfig, args...)`.
function rrule(::typeof(Quaternionic.value), x)
    value_pullback(ΔΩ) = (NoTangent(), NoTangent())
    return Quaternionic.value(x), value_pullback
end
frule(::Tuple, ::typeof(Quaternionic.value), x) = Quaternionic.value(x), NoTangent()
