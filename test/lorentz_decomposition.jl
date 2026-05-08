# Tests for the RB / BR polar decomposition of Lorentz{T} rotors.
#
# Strategy: metamorphic testing.
#   - Round-trip:   Lorentz(R) * B ≈ ±Λ  (and B * Lorentz(R) ≈ ±Λ for BR)
#   - Structural:   R has zero imaginary components; B has zero real vector components
#   - Special cases: identity, pure rotation, pure boost, minus-identity
#   - Multi-precision: error shrinks as precision increases

@testmodule LorentzDecompData begin
    using Random: seed!
    using Quaternionic
    import Quaternionic: RB, BR, ℂconj, ℂreal, ℂimag, ℂreim

    seed!(42)
    const n = 20
    const rand_rotations = [Lorentz(randn(Rotor{Float64})) for _ ∈ 1:n]

    seed!(123)
    const rand_boosts = [
        Boost(abs(randn()) + 0.1, QuatVec(normalize(randn(3))))
        for _ ∈ 1:n
    ]

    # Accumulated products: alternating rotation × boost
    const mixed = accumulate(*, [x for pair ∈ zip(rand_rotations, rand_boosts) for x ∈ pair])

    # ── Predicate helpers ─────────────────────────────────────────────────────

    # True if A and B represent the same element of SO⁺(3,1) (±1 spinor cover)
    function same_Λ(A, B; atol)
        # `norm` acting on a complex quaternion is the spinor norm, which will return a
        # complex number.  We just need to take the absolute value of that to get a
        # real-valued distance metric.
        let norm = abs ∘ norm
            isapprox(A,  B; atol, norm) || isapprox(A, -B; atol, norm)
        end
    end

    # True if R (as a Lorentz rotor) has negligible imaginary parts
    is_real_rotor(R::Rotor{T}; atol) where {T} =
        all(c -> abs(imag(c)) < atol, getfield(Lorentz(R), :components))

    # True if B is a pure boost: imaginary scalar ≈ 0 and real vector ≈ 0
    function is_pure_boost(B::Lorentz{T}; atol) where {T}
        cs = getfield(B, :components)
        abs(imag(cs[1])) < atol &&
        all(c -> abs(real(c)) < atol, cs[2:4])
    end
end

# ── Special cases ─────────────────────────────────────────────────────────────

@testitem "RB/BR: identity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: same_Λ, is_real_rotor, is_pure_boost
    for T ∈ (Float32, Float64)
        ε = 10eps(T)
        Λ = one(Lorentz{T})

        R, B = RB(Λ)
        @test same_Λ(Lorentz(R) * B, Λ; atol=ε)
        @test is_real_rotor(R; atol=ε)
        @test is_pure_boost(B; atol=ε)

        B2, R2 = BR(Λ)
        @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
        @test is_real_rotor(R2; atol=ε)
        @test is_pure_boost(B2; atol=ε)
    end
end

@testitem "RB/BR: minus identity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: same_Λ, is_real_rotor, is_pure_boost
    for T ∈ (Float32, Float64)
        ε = 10eps(T)
        Λ = -one(Lorentz{T})

        R, B = RB(Λ)
        @test same_Λ(Lorentz(R) * B, Λ; atol=ε)
        @test is_real_rotor(R; atol=ε)
        @test is_pure_boost(B; atol=ε)

        B2, R2 = BR(Λ)
        @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
    end
end

@testitem "RB/BR: pure rotation" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: rand_rotations, same_Λ, is_pure_boost
    ε = 64eps(Float64)
    for Λ ∈ rand_rotations
        R, B = RB(Λ)
        @test is_pure_boost(B; atol=ε)            # B ≈ identity boost
        @test same_Λ(Lorentz(R) * B, Λ; atol=ε)

        B2, R2 = BR(Λ)
        @test is_pure_boost(B2; atol=ε)
        @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
    end
end

@testitem "RB/BR: pure boost" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: rand_boosts, same_Λ, is_real_rotor
    ε = 64eps(Float64)
    for Λ ∈ rand_boosts
        R, B = RB(Λ)
        @test is_real_rotor(R; atol=ε)            # R ≈ identity rotation
        @test same_Λ(Lorentz(R) * B, Λ; atol=ε)

        B2, R2 = BR(Λ)
        @test is_real_rotor(R2; atol=ε)
        @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
    end
end

# ── Random round-trips ────────────────────────────────────────────────────────

@testitem "RB round-trip: random mixed" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: RB
    using .LorentzDecompData: mixed, same_Λ, is_real_rotor, is_pure_boost
    ε = 4096eps(Float64)  # products of ~40 elements accumulate ~40× eps error
    for Λ ∈ mixed
        R, B = RB(Λ)
        @test is_real_rotor(R; atol=ε)
        @test is_pure_boost(B; atol=ε)
        @test same_Λ(Lorentz(R) * B, Λ; atol=ε)
    end
end

@testitem "BR round-trip: random mixed" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: BR
    using .LorentzDecompData: mixed, same_Λ, is_real_rotor, is_pure_boost
    ε = 200eps(Float64)
    for Λ ∈ mixed
        B, R = BR(Λ)
        @test is_real_rotor(R; atol=ε)
        @test is_pure_boost(B; atol=ε)
        @test same_Λ(B * Lorentz(R), Λ; atol=ε)
    end
end

@testitem "Float32 round-trip: explicit boost + rotation" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: same_Λ, is_real_rotor, is_pure_boost
    T = Float32
    ε = 128eps(T)
    n̂ = QuatVec{T}(1, 1, 1)/√T(3)
    B0 = Boost(T(1.2), n̂)
    R0 = Lorentz(rotor(T(0.6), T(0.5), T(0.3), T(0.4)))
    Λ  = R0 * B0

    R, B = RB(Λ)
    @test is_real_rotor(R; atol=ε)
    @test is_pure_boost(B; atol=ε)
    @test same_Λ(Lorentz(R) * B, Λ; atol=ε)

    B2, R2 = BR(Λ)
    @test is_real_rotor(R2; atol=ε)
    @test is_pure_boost(B2; atol=ε)
    @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
end

# ── RB and BR equivariance ────────────────────────────────────────────────────

@testitem "equivariance: B_br = R·B_rb·R⁻¹" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: mixed, same_Λ
    using DoubleFloats: Double64
    # Both decompositions share R = rotor(ℜΛ).  The boost factors are related by
    # conjugation: B_br = R * B_rb * R⁻¹ (boost direction rotated to lab frame).
    # Error scales with the boost magnitude cosh(η/2), which bounds it to ~36 eps
    # for simple inputs and ~139 eps for the accumulated `mixed` products.
    for T ∈ (Float32, Float64, Double64)
        ε = 100eps(T)
        for η ∈ T[0.3, 1.0, 2.5, 5.0]
            for n̂ ∈ [QuatVec{T}(1,0,0), QuatVec{T}(0,1,0), QuatVec(T(1)/√T(3), T(1)/√T(3), T(1)/√T(3))]
                for θ ∈ T[π/6, π/3, 2π/3]
                    Λ = Lorentz(rotor(cos(θ/2), 0, 0, sin(θ/2))) * Boost(η, n̂)
                    R, B_rb = RB(Λ)
                    B_br, _ = BR(Λ)
                    @test same_Λ(B_br, Lorentz(R) * B_rb * Lorentz(conj(R)); atol=ε)
                end
            end
        end
    end
    # Float64 accumulated products accumulate more error (~139 eps max on this data)
    for Λ ∈ mixed
        R, B_rb = RB(Λ)
        B_br, _ = BR(Λ)
        @test same_Λ(B_br, Lorentz(R) * B_rb * Lorentz(conj(R)); atol=256eps(Float64))
    end
end

# ── Multi-precision convergence ───────────────────────────────────────────────

@testitem "multi-precision convergence" tags=[:validation, :slow] begin
    import Quaternionic: RB
    using DoubleFloats: Double64

    function round_trip_error(T)
        θ, η = T(1.3), T(0.8)
        n̂T = [1/√T(3), 1/√T(3), 1/√T(3)]
        B0 = Boost(η, n̂T)
        R0 = Lorentz(rotor(cos(θ/2), sin(θ/2),0, 0))
        Λ  = R0 * B0
        R, B = RB(Λ)
        maximum(abs, components(Lorentz(R) * B) - components(Λ))
    end

    err32  = round_trip_error(Float32)
    err64  = round_trip_error(Float64)
    errD64 = round_trip_error(Double64)
    errBig = round_trip_error(BigFloat)

    # Algorithm is algebraically direct; errors should scale with machine epsilon.
    # We do not compare across precisions because Float32 can achieve zero error
    # on simple inputs (exact floating-point results at each step).
    @test err32  ≤ 10eps(Float32)
    @test err64  ≤ 10eps(Float64)
    @test errD64 ≤ 10eps(Double64)
    @test errBig ≤ 10eps(BigFloat)
end

# ── Rv and vR ─────────────────────────────────────────────────────────────────

@testitem "Rv/vR: identity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: Rv, vR
    using .LorentzDecompData: same_Λ
    for T ∈ (Float32, Float64)
        ε = 10eps(T)
        Λ = one(Lorentz{T})

        R, v⃗ = Rv(Λ)
        @test same_Λ(Lorentz(R), Λ; atol=ε)
        @test norm(v⃗) ≤ ε

        v⃗2, R2 = vR(Λ)
        @test same_Λ(Lorentz(R2), Λ; atol=ε)
        @test norm(v⃗2) ≤ ε
    end
end

@testitem "Rv/vR: pure rotation → zero velocity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: Rv, vR
    using .LorentzDecompData: rand_rotations, same_Λ
    ε = 64eps(Float64)
    for Λ ∈ rand_rotations
        R, v⃗ = Rv(Λ)
        @test norm(v⃗) ≤ ε
        @test same_Λ(Lorentz(R), Λ; atol=ε)

        v⃗2, R2 = vR(Λ)
        @test norm(v⃗2) ≤ ε
        @test same_Λ(Lorentz(R2), Λ; atol=ε)
    end
end

@testitem "Rv/vR: pure boost → identity rotation and known velocity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: Rv, vR
    using .LorentzDecompData: is_real_rotor
    for T ∈ (Float32, Float64)
        ε = 64eps(T)
        for η ∈ T[0.375, 0.75, 1.5, 2.0]
            β = tanh(η)
            for n̂ ∈ [
                QuatVec{T}(1, 0, 0),
                QuatVec{T}(0, 1, 0),
                QuatVec{T}(0, 0, 1),
                QuatVec{T}(1, 1, 1)/√T(3),
            ]
                Λ = Boost(η, n̂)

                R, v⃗ = Rv(Λ)
                @test is_real_rotor(R; atol=ε)
                @test norm(v⃗ - β * n̂) ≤ ε

                v⃗2, R2 = vR(Λ)
                @test is_real_rotor(R2; atol=ε)
                @test norm(v⃗2 - β * n̂) ≤ ε
            end
        end
    end
end

@testitem "equivariance: R·v_Rv·R⁻¹ = v_vR" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: Rv, vR
    using .LorentzDecompData: mixed, rand_rotations, rand_boosts
    using DoubleFloats: Double64
    # Both Rv and vR share R = rotor(ℜΛ) (bit-for-bit identical).  The velocities
    # are related by the adjoint action of R: v_vR = R * v_Rv * R⁻¹.
    # Holds to ~2 eps: velocities lie in [0,1), so there is no magnitude amplification.
    for T ∈ (Float32, Float64, Double64)
        ε = 10eps(T)
        for η ∈ T[0.3, 1.0, 2.5, 5.0]
            for n̂ ∈ [QuatVec{T}(1,0,0), QuatVec{T}(0,1,0), QuatVec{T}(1, 1, 1)/√T(3)]
                for θ ∈ T[π/6, π/3, 2π/3]
                    Λ = Lorentz(rotor(cos(θ/2), zero(T), zero(T), sin(θ/2))) * Boost(η, n̂)
                    R, v_rv = Rv(Λ)
                    v_vr, R2 = vR(Λ)
                    @test R == R2
                    @test norm(v_vr - QuatVec(R * v_rv * conj(R))) ≤ ε
                end
            end
        end
    end
    # Float64 accumulated products and random inputs: equivariance still holds to ~2 eps
    for Λ ∈ [mixed; rand_rotations; rand_boosts]
        R, v_rv = Rv(Λ)
        v_vr, R2 = vR(Λ)
        @test R == R2
        @test norm(v_vr - QuatVec(R * v_rv * conj(R))) ≤ 10eps(Float64)
    end
end

@testitem "Rv/vR multi-precision convergence" tags=[:validation, :slow] begin
    import Quaternionic: Rv, vR
    using DoubleFloats: Double64

    function rv_error(T)
        θ, η = T(1.3), T(0.8)
        n̂ = QuatVec{T}(1, 1, 1)/√T(3)
        Λ = Lorentz(rotor(cos(θ/2), sin(θ/2), 0, 0)) * Boost(η, n̂)
        R, v⃗ = Rv(Λ)
        β = norm(v⃗)
        B = Boost(atanh(β), v⃗ * inv(β))
        maximum(abs, components(Lorentz(R) * B) - components(Λ))
    end

    function vr_error(T)
        θ, η = T(1.3), T(0.8)
        n̂ = QuatVec{T}(1, 1, 1)/√T(3)
        Λ = Lorentz(rotor(cos(θ/2), sin(θ/2), 0, 0)) * Boost(η, n̂)
        v⃗, R = vR(Λ)
        β = norm(v⃗)
        B = Boost(atanh(β), v⃗ * inv(β))
        maximum(abs, components(B * Lorentz(R)) - components(Λ))
    end

    @test rv_error(Float32)  ≤ 10eps(Float32)
    @test rv_error(Float64)  ≤ 10eps(Float64)
    @test rv_error(Double64) ≤ 10eps(Double64)
    @test rv_error(BigFloat) ≤ 10eps(BigFloat)

    @test vr_error(Float32)  ≤ 10eps(Float32)
    @test vr_error(Float64)  ≤ 10eps(Float64)
    @test vr_error(Double64) ≤ 10eps(Double64)
    @test vr_error(BigFloat) ≤ 10eps(BigFloat)
end
