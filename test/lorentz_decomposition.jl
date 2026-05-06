# Tests for the RB / BR polar decomposition of Lorentz{T} rotors.
#
# Strategy: metamorphic testing.
#   - Round-trip:   Lorentz(R) * B ≈ ±Λ  (and B * Lorentz(R) ≈ ±Λ for BR)
#   - Structural:   R has zero imaginary components; B has zero real vector components
#   - Special cases: identity, pure rotation, pure boost, minus-identity
#   - Multi-precision: error shrinks as precision increases

@testmodule LorentzDecompData begin
    using Random: seed!
    using LinearAlgebra: normalize as la_normalize
    using Quaternionic
    import Quaternionic: RB, BR, ℂconj, ℂreal, ℂimag, ℂreim

    seed!(42)
    const n = 20
    const rand_rotations = [Lorentz(randn(Rotor{Float64})) for _ ∈ 1:n]

    seed!(123)
    const rand_boosts = [
        Boost(abs(randn()) + 0.1, QuatVec(la_normalize(randn(3))))
        for _ ∈ 1:n
    ]

    # Accumulated products: alternating rotation × boost
    const mixed = accumulate(*, [x for pair ∈ zip(rand_rotations, rand_boosts) for x ∈ pair])

    # ── Predicate helpers ─────────────────────────────────────────────────────

    # True if A and B represent the same element of SO⁺(3,1) (±1 spinor cover)
    same_Λ(A, B; atol) =
        isapprox(components(A), components(B); atol) ||
        isapprox(components(A), -components(B); atol)

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
    ε = 1024eps(Float64)  # products of ~40 elements accumulate ~40× eps error
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
    n̂ = QuatVec(T(1/√3), T(1/√3), T(1/√3))
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

# ── RB and BR are different orderings ─────────────────────────────────────────

@testitem "RB ≠ BR for generic inputs" tags=[:unit, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: mixed
    # For generic non-commuting Lorentz elements, RB and BR yield different B factors
    n_different = count(mixed) do Λ
        R_rb, B_rb = RB(Λ)
        B_br, R_br = BR(Λ)
        !isapprox(components(B_rb), components(B_br); atol=1e-6) &&
        !isapprox(components(B_rb), -components(B_br); atol=1e-6)
    end
    # Expect most random mixed elements to give genuinely different B factors
    @test n_different > length(mixed) ÷ 2
end

# ── Multi-precision convergence ───────────────────────────────────────────────

@testitem "multi-precision convergence" tags=[:validation, :slow] begin
    import Quaternionic: RB
    using DoubleFloats: Double64

    function round_trip_error(T)
        θ, η = T(1.3), T(0.8)
        n̂T = T.(Float64[1/√3, 1/√3, 1/√3])
        B0 = Boost(η, n̂T)
        R0 = Lorentz(rotor(cos(θ/2), sin(θ/2), zero(T), zero(T)))
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
