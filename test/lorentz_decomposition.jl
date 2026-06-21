# Tests for the RB / BR polar decomposition of Lorentz{T} rotors.
#
# Strategy: metamorphic testing.
#   - Round-trip:   Lorentz(R) * B ≈ ±Λ  (and B * Lorentz(R) ≈ ±Λ for BR)
#   - Structural:   R has zero imaginary components; B has zero real vector components
#   - Special cases: identity, pure rotation, pure boost, minus-identity
#   - Multi-precision: error shrinks as precision increases

@testmodule LorentzDecompData begin
    using Random: Xoshiro
    using Quaternionic
    using DoubleFloats: Double64
    import Quaternionic: RB, BR, KAN, ℂconj, ℂreal, ℂimag, ℂreim

    const n = 20

    # Raw data generated at Double64 precision so all float types see the same
    # rapidity values (different types produce different randn sequences).
    const _rot_spinors = let rng = Xoshiro(42)
        [randn(rng, Rotor{Double64}) for _ ∈ 1:n]
    end

    const _boost_data = let rng = Xoshiro(123)
        [(min(abs(randn(rng, Double64)), Double64(2)) + Double64(0.1),
          QuatVec(normalize(randn(rng, Double64, 3))))
         for _ ∈ 1:n]
    end

    rand_rotations(T) = [Lorentz(Rotor{T}(T.(components(r))...)) for r ∈ _rot_spinors]
    rand_boosts(T)    = [Boost(T(η), QuatVec{T}(T.(components(n̂))...)) for (η, n̂) ∈ _boost_data]

    # Accumulated products: alternating rotation × boost
    function mixed(T)
        seq = [x for pair ∈ zip(rand_rotations(T), rand_boosts(T)) for x ∈ pair]
        accumulate(*, seq)
    end

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

    # ── KAN-specific helpers ───────────────────────────────────────────────────

    # A null rotation in N (fixing ℓ = (𝐭+𝐳)/√2), built from screw parameters ξˣ, ξʸ
    # using the same encoding as `KAN` returns for its Rₙ factor.
    function null_rotation(::Type{T}, ξˣ, ξʸ) where {T}
        a, b = T(ξˣ)/sqrt(T(2)), T(ξʸ)/sqrt(T(2))
        Lorentz{T}(1, (im*a + b)/2, (im*b - a)/2, 0)
    end

    # A boost along the z-axis with rapidity φ — an element of A.
    z_boost(::Type{T}, φ) where {T} = Boost(T(φ), QuatVec{T}(0, 0, 1))

    # The null vector ℓ = (𝐭+𝐳)/√2 fixed by N, as a Minkowski 4-vector.
    ℓ_vec(::Type{T}) where {T} = T[1, 0, 0, 1] / sqrt(T(2))

    # Structural predicate for A: pure z-boost (w real, x=y=0, z pure imaginary).
    # Accepts any AbstractQuaternion (KAN returns Rₐ as a bare `Quaternion`).
    function is_z_boost(R; atol)
        w, x, y, z = components(R)
        abs(imag(w)) < atol && abs(x) < atol && abs(y) < atol && abs(real(z)) < atol
    end

    # Structural predicate for N: null rotation (w=1, z=0, Re(x)=Im(y), Im(x)=−Re(y)).
    function is_null_rotation(R; atol)
        w, x, y, z = components(R)
        abs(w - 1) < atol && abs(z) < atol &&
        abs(real(x) - imag(y)) < atol && abs(imag(x) + real(y)) < atol
    end
end

# ── Special cases ─────────────────────────────────────────────────────────────

@testitem "RB/BR: identity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: same_Λ, is_real_rotor, is_pure_boost
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
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
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
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
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        ε = 64eps(T)
        for Λ ∈ rand_rotations(T)
            R, B = RB(Λ)
            @test is_pure_boost(B; atol=ε)            # B ≈ identity boost
            @test same_Λ(Lorentz(R) * B, Λ; atol=ε)

            B2, R2 = BR(Λ)
            @test is_pure_boost(B2; atol=ε)
            @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
        end
    end
end

@testitem "RB/BR: pure boost" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: RB, BR
    using .LorentzDecompData: rand_boosts, same_Λ, is_real_rotor
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        ε = 64eps(T)
        for Λ ∈ rand_boosts(T)
            R, B = RB(Λ)
            @test is_real_rotor(R; atol=ε)            # R ≈ identity rotation
            @test same_Λ(Lorentz(R) * B, Λ; atol=ε)

            B2, R2 = BR(Λ)
            @test is_real_rotor(R2; atol=ε)
            @test same_Λ(B2 * Lorentz(R2), Λ; atol=ε)
        end
    end
end

# ── Random round-trips ────────────────────────────────────────────────────────

@testitem "RB round-trip: random mixed" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: RB
    using .LorentzDecompData: mixed, same_Λ, is_real_rotor, is_pure_boost
    using DoubleFloats: Double64
    # Errors scale as N×eps(T)×cosh(η/2): N≤109 for structural checks, N≤23 for round-trip.
    for T ∈ (Float32, Float64, Double64)
        for Λ ∈ mixed(T)
            R, B = RB(Λ)
            cosh_half = abs(real(getfield(B, :components)[1]))
            @test is_real_rotor(R; atol=128eps(T)*cosh_half)
            @test is_pure_boost(B; atol=128eps(T)*cosh_half)
            @test same_Λ(Lorentz(R) * B, Λ; atol=32eps(T)*cosh_half)
        end
    end
end

@testitem "BR round-trip: random mixed" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: BR
    using .LorentzDecompData: mixed, same_Λ, is_real_rotor, is_pure_boost
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        for Λ ∈ mixed(T)
            B, R = BR(Λ)
            cosh_half = abs(real(getfield(B, :components)[1]))
            @test is_real_rotor(R; atol=128eps(T)*cosh_half)
            @test is_pure_boost(B; atol=128eps(T)*cosh_half)
            @test same_Λ(B * Lorentz(R), Λ; atol=32eps(T)*cosh_half)
        end
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
            for n̂ ∈ [QuatVec{T}(1,0,0), QuatVec{T}(0,1,0), QuatVec{T}(1,1,1)/√T(3)]
                for θ ∈ T[π/6, π/3, 2π/3]
                    Λ = Lorentz(rotor(cos(θ/2), 0, 0, sin(θ/2))) * Boost(η, n̂)
                    R, B_rb = RB(Λ)
                    B_br, _ = BR(Λ)
                    @test same_Λ(B_br, Lorentz(R) * B_rb * Lorentz(conj(R)); atol=ε)
                end
            end
        end
        # Accumulated products: error scales as N×eps(T)×cosh(η/2), N≤16
        for Λ ∈ mixed(T)
            R, B_rb = RB(Λ)
            B_br, _ = BR(Λ)
            cosh_half = abs(real(getfield(B_rb, :components)[1]))
            @test same_Λ(B_br, Lorentz(R) * B_rb * Lorentz(conj(R)); atol=32eps(T)*cosh_half)
        end
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
        R0 = Lorentz(rotor(cos(θ/2), sin(θ/2), 0, 0))
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
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
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
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        ε = 64eps(T)
        for Λ ∈ rand_rotations(T)
            R, v⃗ = Rv(Λ)
            @test norm(v⃗) ≤ ε
            @test same_Λ(Lorentz(R), Λ; atol=ε)

            v⃗2, R2 = vR(Λ)
            @test norm(v⃗2) ≤ ε
            @test same_Λ(Lorentz(R2), Λ; atol=ε)
        end
    end
end

@testitem "Rv/vR: pure boost → identity rotation and known velocity" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: Rv, vR
    using .LorentzDecompData: is_real_rotor
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
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
            for n̂ ∈ [QuatVec{T}(1,0,0), QuatVec{T}(0,1,0), QuatVec{T}(1,1,1)/√T(3)]
                for θ ∈ T[π/6, π/3, 2π/3]
                    Λ = Lorentz(rotor(cos(θ/2), 0, 0, sin(θ/2))) * Boost(η, n̂)
                    R, v_rv = Rv(Λ)
                    v_vr, R2 = vR(Λ)
                    @test R == R2
                    @test norm(v_vr - QuatVec(R * v_rv * conj(R))) ≤ ε
                end
            end
        end
        # Accumulated products and random inputs: equivariance still holds to ~2 eps
        for Λ ∈ [mixed(T); rand_rotations(T); rand_boosts(T)]
            R, v_rv = Rv(Λ)
            v_vr, R2 = vR(Λ)
            @test R == R2
            @test norm(v_vr - QuatVec(R * v_rv * conj(R))) ≤ 10eps(T)
        end
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

# ── KAN (Iwasawa) decomposition ────────────────────────────────────────────────
#
# KAN(Λ) → (Rₖ, Rₐ, Rₙ) with Λ = Rₖ·Rₐ·Rₙ, where Rₖ ∈ K is a rotation, Rₐ ∈ A is a
# boost along 𝐳, and Rₙ ∈ N is a null rotation fixing ℓ = (𝐭+𝐳)/√2.  Errors scale
# with the overall transformation magnitude, captured by `maximum(abs, components(Λ))`.

@testitem "KAN: reconstruction round-trip" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: KAN, components
    using .LorentzDecompData: mixed, rand_rotations, rand_boosts, same_Λ
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        for Λ ∈ [mixed(T); rand_rotations(T); rand_boosts(T)]
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            ε = 1024eps(T) * maximum(abs, components(Λ))
            @test same_Λ(Lorentz(Rₖ) * Lorentz(Rₐ) * Lorentz(Rₙ), Λ; atol=ε)
        end
    end
end

@testitem "KAN: structural membership (K, A, N)" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: KAN, components
    using .LorentzDecompData: mixed, rand_rotations, rand_boosts, is_real_rotor, is_z_boost, is_null_rotation
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        for Λ ∈ [mixed(T); rand_rotations(T); rand_boosts(T)]
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            ε = 1024eps(T) * maximum(abs, components(Λ))
            @test is_real_rotor(Rₖ; atol=ε)        # Rₖ ∈ K: pure rotation (no boost part)
            @test is_z_boost(Rₐ; atol=ε)           # Rₐ ∈ A: boost along 𝐳 only
            @test is_null_rotation(Rₙ; atol=ε)     # Rₙ ∈ N: null rotation, w=1, z=0
        end
    end
end

@testitem "KAN: fixed-point oracles" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: KAN, components
    using .LorentzDecompData: mixed, rand_rotations, rand_boosts, ℓ_vec
    using DoubleFloats: Double64
    # K fixes the time axis 𝐭; A (boost along 𝐳) fixes the transverse axes 𝐱, 𝐲;
    # N (null rotation) fixes the null vector ℓ = (𝐭+𝐳)/√2.
    for T ∈ (Float32, Float64, Double64)
        𝐭, 𝐱, 𝐲 = T[1,0,0,0], T[0,1,0,0], T[0,0,1,0]
        ℓ = ℓ_vec(T)
        for Λ ∈ [mixed(T); rand_rotations(T); rand_boosts(T)]
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            ε = 1024eps(T) * maximum(abs, components(Λ))
            @test isapprox(Lorentz(Rₖ)(𝐭), 𝐭; atol=ε)
            @test isapprox(Lorentz(Rₐ)(𝐱), 𝐱; atol=ε)
            @test isapprox(Lorentz(Rₐ)(𝐲), 𝐲; atol=ε)
            @test isapprox(Lorentz(Rₙ)(ℓ), ℓ; atol=ε)
        end
    end
end

@testitem "KAN: special cases" tags=[:unit, :fast, :validation] setup=[LorentzDecompData] begin
    import Quaternionic: KAN
    using .LorentzDecompData: rand_rotations, null_rotation, z_boost, same_Λ
    using DoubleFloats: Double64
    for T ∈ (Float32, Float64, Double64)
        ε = 64eps(T)
        𝟙 = one(Lorentz{T})

        # Identity and minus-identity → all three factors trivial
        for Λ ∈ (𝟙, -𝟙)
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            @test same_Λ(Lorentz(Rₖ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₐ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₙ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₖ) * Lorentz(Rₐ) * Lorentz(Rₙ), Λ; atol=ε)
        end

        # Pure rotation → Rₖ ≈ Λ, Rₐ ≈ Rₙ ≈ 𝟙
        for Λ ∈ rand_rotations(T)
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            @test same_Λ(Lorentz(Rₖ), Λ; atol=ε)
            @test same_Λ(Lorentz(Rₐ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₙ), 𝟙; atol=ε)
        end

        # Pure z-boost → Rₐ ≈ Λ, Rₖ ≈ Rₙ ≈ 𝟙
        for φ ∈ T[0.3, 1.0, 2.5]
            Λ = z_boost(T, φ)
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            @test same_Λ(Lorentz(Rₖ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₐ), Λ; atol=ε)
            @test same_Λ(Lorentz(Rₙ), 𝟙; atol=ε)
        end

        # Pure null rotation → Rₙ ≈ Λ, Rₖ ≈ Rₐ ≈ 𝟙
        for (ξˣ, ξʸ) ∈ [(T(0.4), T(-0.7)), (T(1.5), T(0.0)), (T(0.0), T(2.0))]
            Λ = null_rotation(T, ξˣ, ξʸ)
            εₙ = 64eps(T) * maximum(abs, Quaternionic.components(Λ))
            Rₖ, Rₐ, Rₙ = KAN(Λ)
            @test same_Λ(Lorentz(Rₖ), 𝟙; atol=εₙ)
            @test same_Λ(Lorentz(Rₐ), 𝟙; atol=εₙ)
            @test same_Λ(Lorentz(Rₙ), Λ; atol=εₙ)
        end
    end
end

@testitem "KAN: factor recovery (uniqueness)" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: KAN, components
    using .LorentzDecompData: rand_rotations, null_rotation, z_boost, same_Λ
    using DoubleFloats: Double64
    # Build Λ from known canonical factors K₀·A₀·N₀ and verify KAN recovers each.
    # The Iwasawa decomposition is unique, so recovery is an oracle test.
    for T ∈ (Float32, Float64, Double64)
        for K₀ ∈ rand_rotations(T)
            for φ ∈ T[0.5, 2.0]
                A₀ = z_boost(T, φ)
                for (ξˣ, ξʸ) ∈ [(T(0.4), T(-0.7)), (T(1.3), T(0.9))]
                    N₀ = null_rotation(T, ξˣ, ξʸ)
                    Λ = K₀ * A₀ * N₀
                    ε = 1024eps(T) * maximum(abs, components(Λ))
                    Rₖ, Rₐ, Rₙ = KAN(Λ)
                    @test same_Λ(Lorentz(Rₖ), K₀; atol=ε)
                    @test same_Λ(Lorentz(Rₐ), A₀; atol=ε)
                    @test same_Λ(Lorentz(Rₙ), N₀; atol=ε)
                end
            end
        end
    end
end

@testitem "KAN: idempotency (decomposing a factor)" tags=[:validation, :fast] setup=[LorentzDecompData] begin
    import Quaternionic: KAN, components
    using .LorentzDecompData: rand_rotations, null_rotation, z_boost, same_Λ
    using DoubleFloats: Double64
    # Decomposing an element that already lives in a single subgroup returns it in
    # the matching slot and the identity in the other two.
    for T ∈ (Float32, Float64, Double64)
        𝟙 = one(Lorentz{T})
        ε = 64eps(T)

        for K₀ ∈ rand_rotations(T)
            Rₖ, Rₐ, Rₙ = KAN(K₀)
            @test same_Λ(Lorentz(Rₖ), K₀; atol=ε)
            @test same_Λ(Lorentz(Rₐ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₙ), 𝟙; atol=ε)
        end

        for φ ∈ T[0.5, 2.0]
            A₀ = z_boost(T, φ)
            Rₖ, Rₐ, Rₙ = KAN(A₀)
            @test same_Λ(Lorentz(Rₖ), 𝟙; atol=ε)
            @test same_Λ(Lorentz(Rₐ), A₀; atol=ε)
            @test same_Λ(Lorentz(Rₙ), 𝟙; atol=ε)
        end

        for (ξˣ, ξʸ) ∈ [(T(0.4), T(-0.7)), (T(1.3), T(0.9))]
            N₀ = null_rotation(T, ξˣ, ξʸ)
            εₙ = 64eps(T) * maximum(abs, components(N₀))
            Rₖ, Rₐ, Rₙ = KAN(N₀)
            @test same_Λ(Lorentz(Rₖ), 𝟙; atol=εₙ)
            @test same_Λ(Lorentz(Rₐ), 𝟙; atol=εₙ)
            @test same_Λ(Lorentz(Rₙ), N₀; atol=εₙ)
        end
    end
end

@testitem "KAN: multi-precision convergence" tags=[:validation, :slow] begin
    import Quaternionic: KAN, components
    using DoubleFloats: Double64
    # Algebraically direct algorithm: reconstruction error stays at the machine-epsilon
    # floor across precisions (no Float64-specific constants or thresholds).
    function recon_error(T)
        c = sin(T(13)/20) / √T(3)
        K₀ = Lorentz(rotor(cos(T(13)/20), c, c, c))
        A₀ = Boost(T(8)/10, QuatVec{T}(0, 0, 1))
        a, b = T(4)/10/√T(2), T(-7)/10/√T(2)
        N₀ = Lorentz{T}(1, (im*a + b)/2, (im*b - a)/2, 0)
        Λ = K₀ * A₀ * N₀
        Rₖ, Rₐ, Rₙ = KAN(Λ)
        recon = Lorentz(Rₖ) * Lorentz(Rₐ) * Lorentz(Rₙ)
        maximum(abs, components(recon) - components(Λ))
    end
    @test recon_error(Float32)  ≤ 10eps(Float32)
    @test recon_error(Float64)  ≤ 10eps(Float64)
    @test recon_error(Double64) ≤ 10eps(Double64)
    @test recon_error(BigFloat) ≤ 10eps(BigFloat)
end
