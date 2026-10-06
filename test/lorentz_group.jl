# Tests for Lorentz{T} — ported from Scri.jl/test/test-lorentz.jl
#
# Strategy: metamorphic testing.  Mathematical properties are factored into
# predicate functions, then called for many random inputs.

"""Minkowski inner product with signature −+++."""
minkowski_inner(v, w) = -v[1]*w[1] + v[2]*w[2] + v[3]*w[3] + v[4]*w[4]

"""Minkowski norm squared (signed; negative for timelike vectors)."""
minkowski_norm²(v) = minkowski_inner(v, v)

metric_preserved(Λ, v, w; atol=1e-12) =
    abs(minkowski_inner(Λ(v), Λ(w)) - minkowski_inner(v, w)) ≤ atol

composition_consistent(Λ₁, Λ₂, v; atol=1e-12) =
    norm((Λ₁ * Λ₂)(v) - Λ₁(Λ₂(v))) ≤ atol

inverse_law(Λ, v; atol=1e-12) =
    norm((Λ * inv(Λ))(v) - v) ≤ atol && norm((inv(Λ) * Λ)(v) - v) ≤ atol

identity_law(Λ, v; atol=1e-12) = norm(one(Λ)(v) - v) ≤ atol

associativity_law(Λ₁, Λ₂, Λ₃, v; atol=1e-12) =
    norm(((Λ₁ * Λ₂) * Λ₃)(v) - (Λ₁ * (Λ₂ * Λ₃))(v)) ≤ atol

preserves_time(Λ, v; atol=1e-12) = abs(Λ(v)[1] - v[1]) ≤ atol

preserves_spatial_norm(Λ, v; atol=1e-12) =
    abs(norm(Λ(v)[2:4]) - norm(v[2:4])) ≤ atol

double_cover(Λ_pos, Λ_neg, v; atol=1e-12) =
    norm(Λ_pos(v) - Λ_neg(v)) ≤ atol

rotor_homomorphism(Λ₁₂, Λ₁, Λ₂, v; atol=1e-12) =
    norm(Λ₁₂(v) - (Λ₁ * Λ₂)(v)) ≤ atol

function ga_norm_conditions(Λ; atol=1e-12)
    R¹, Rᶻʸ, Rˣᶻ, Rʸˣ, Rᵗˣ, Rᵗʸ, Rᵗᶻ, Rᵗˣʸᶻ = ga_components(Λ)
    quad  = R¹^2 + Rᶻʸ^2 + Rˣᶻ^2 + Rʸˣ^2 - Rᵗˣ^2 - Rᵗʸ^2 - Rᵗᶻ^2 - Rᵗˣʸᶻ^2
    cross = R¹ * Rᵗˣʸᶻ + Rᶻʸ * Rᵗˣ + Rˣᶻ * Rᵗʸ + Rʸˣ * Rᵗᶻ
    return abs(quad - 1) ≤ atol && abs(cross) ≤ atol
end

function ga_reverse_components(Λ; atol=1e-12)
    c  = ga_components(Λ)
    ci = ga_components(inv(Λ))
    abs(ci[1] - c[1]) ≤ atol &&   # R¹     unchanged  (grade 0)
    abs(ci[2] + c[2]) ≤ atol &&   # Rᶻʸ    negated    (grade 2)
    abs(ci[3] + c[3]) ≤ atol &&   # Rˣᶻ    negated    (grade 2)
    abs(ci[4] + c[4]) ≤ atol &&   # Rʸˣ    negated    (grade 2)
    abs(ci[5] + c[5]) ≤ atol &&   # Rᵗˣ    negated    (grade 2)
    abs(ci[6] + c[6]) ≤ atol &&   # Rᵗʸ    negated    (grade 2)
    abs(ci[7] + c[7]) ≤ atol &&   # Rᵗᶻ    negated    (grade 2)
    abs(ci[8] - c[8]) ≤ atol      # Rᵗˣʸᶻ  unchanged  (grade 4)
end

# ─── Shared random data ────────────────────────────────────────────────────────

Random.seed!(42)
const lg_T = Float64
const lg_n = 20

lg_rot_rotors  = [randn(Rotor{lg_T}) for _ ∈ 1:lg_n]
lg_rot_Λs      = [Lorentz{lg_T}(R) for R ∈ lg_rot_rotors]
lg_gen_vecs    = [randn(lg_T, 4) for _ ∈ 1:lg_n]
lg_spatial_vecs = [[zero(lg_T); randn(lg_T, 3)] for _ ∈ 1:lg_n]

Random.seed!(123)
lg_boost_rapidities  = abs.(randn(lg_T, lg_n)) .+ lg_T(0.1)
lg_boost_directions  = [QuatVec(normalize(randn(lg_T, 3))) for _ ∈ 1:lg_n]
lg_boost_Λs          = [Boost(η, n̂) for (η, n̂) ∈ zip(lg_boost_rapidities, lg_boost_directions)]

lg_mixed_seq  = [x for pair ∈ zip(lg_rot_Λs, lg_boost_Λs) for x ∈ pair]
lg_composed   = accumulate(*, lg_mixed_seq)

@testset "Lorentz group" begin

    # ── Group structure: rotations ─────────────────────────────────────────────

    @testset "Rotation: identity element" begin
        for v ∈ lg_gen_vecs, Λ ∈ lg_rot_Λs
            @test identity_law(Λ, v)
        end
    end

    @testset "Rotation: composition" begin
        for i ∈ 1:(lg_n-1), v ∈ lg_gen_vecs
            @test composition_consistent(lg_rot_Λs[i], lg_rot_Λs[i+1], v)
        end
    end

    @testset "Rotation: inverse" begin
        for Λ ∈ lg_rot_Λs, v ∈ lg_gen_vecs
            @test inverse_law(Λ, v)
        end
    end

    @testset "Rotation: associativity" begin
        for i ∈ 1:(lg_n-2), v ∈ lg_gen_vecs
            @test associativity_law(lg_rot_Λs[i], lg_rot_Λs[i+1], lg_rot_Λs[i+2], v)
        end
    end

    # ── Minkowski isometry: rotations ──────────────────────────────────────────

    @testset "Rotation: preserves Minkowski metric" begin
        for Λ ∈ lg_rot_Λs, i ∈ 1:(lg_n-1)
            @test metric_preserved(Λ, lg_gen_vecs[i], lg_gen_vecs[i+1])
        end
    end

    @testset "Rotation: null vectors remain null" begin
        for Λ ∈ lg_rot_Λs, v ∈ lg_spatial_vecs
            v_sp = v[2:4]
            ℓ = [norm(v_sp); v_sp]
            ℓ′ = Λ(ℓ)
            @test abs(minkowski_norm²(ℓ′)) ≤ 1e-12
        end
    end

    # ── Rotation-specific invariants ───────────────────────────────────────────

    @testset "Rotation: preserves time component" begin
        for Λ ∈ lg_rot_Λs, v ∈ lg_gen_vecs
            @test preserves_time(Λ, v)
        end
    end

    @testset "Rotation: preserves spatial norm" begin
        for Λ ∈ lg_rot_Λs, v ∈ lg_spatial_vecs
            @test preserves_spatial_norm(Λ, v)
        end
    end

    # ── Double cover ───────────────────────────────────────────────────────────

    @testset "Rotation: double cover — R and −R give the same transformation" begin
        for (R, v) ∈ zip(lg_rot_rotors, lg_gen_vecs)
            @test double_cover(Lorentz{lg_T}(R), Lorentz{lg_T}(-R), v)
        end
    end

    # ── Spin(3) → Lorentz homomorphism ────────────────────────────────────────

    @testset "Rotation: Spin(3) → Lorentz is a group homomorphism" begin
        for i ∈ 1:(lg_n-1)
            R₁, R₂ = lg_rot_rotors[i], lg_rot_rotors[i+1]
            Λ₁₂ = Lorentz{lg_T}(R₁ * R₂)
            Λ₁  = Lorentz{lg_T}(R₁)
            Λ₂  = Lorentz{lg_T}(R₂)
            for v ∈ lg_gen_vecs
                @test rotor_homomorphism(Λ₁₂, Λ₁, Λ₂, v)
            end
        end
    end

    # ── Axis-fixing ────────────────────────────────────────────────────────────

    @testset "Rotation: rotation about an axis fixes that axis direction" begin
        for θ ∈ [0.3, 1.0, π/2, π, 2π - 0.1]
            R = Rotor(cos(θ/2), 0.0, 0.0, sin(θ/2))
            Λ = Lorentz{Float64}(R)
            @test Λ([0.0, 0.0, 0.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0] atol=1e-12
            v_xy = [0.0, 1.0, 0.5, 0.0]
            v_xy′ = Λ(v_xy)
            @test abs(v_xy′[4]) ≤ 1e-12
            @test abs(v_xy′[1]) ≤ 1e-12
        end
        for θ ∈ [0.7, π/3]
            R = Rotor(cos(θ/2), sin(θ/2), 0.0, 0.0)
            Λ = Lorentz{Float64}(R)
            @test Λ([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol=1e-12
        end
        for θ ∈ [1.2, π/4]
            R = Rotor(cos(θ/2), 0.0, sin(θ/2), 0.0)
            Λ = Lorentz{Float64}(R)
            @test Λ([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, 0.0, 1.0, 0.0] atol=1e-12
        end
    end

    @testset "Rotation: composing with itself doubles the rotation angle" begin
        for θ ∈ [0.3, 0.9, 1.5, π/3]
            R_θ  = Rotor(Quaternion(cos(θ/2), 0.0, 0.0, sin(θ/2)))
            R_2θ = Rotor(Quaternion(cos(θ),   0.0, 0.0, sin(θ)))
            Λ_θ  = Lorentz{Float64}(R_θ)
            Λ_2θ = Lorentz{Float64}(R_2θ)
            v = [0.0, 1.0, 0.0, 0.0]
            @test Λ_2θ(v) ≈ (Λ_θ * Λ_θ)(v) atol=1e-12
        end
    end

    @testset "Rotation: 2π rotation acts as the identity on 4-vectors" begin
        Random.seed!(7)
        R_2π = Rotor(-1.0, 0.0, 0.0, 0.0)
        Λ_2π = Lorentz{Float64}(R_2π)
        for _ ∈ 1:10
            v = randn(4)
            @test Λ_2π(v) ≈ v atol=1e-12
        end
    end

    # ── Boost constructors and GA algebraic properties ─────────────────────────

    @testset "Boost: explicit GA component values" begin
        for η ∈ [0.3, 0.7, 1.2, 1.8, 2.5]
            ch, sh = cosh(η/2), sinh(η/2)

            c = ga_components(Boost(η, [1.0, 0.0, 0.0]))
            @test c[1] ≈ ch  atol=1e-14   # R¹
            @test c[2] ≈ 0.0 atol=1e-14   # Rᶻʸ
            @test c[3] ≈ 0.0 atol=1e-14   # Rˣᶻ
            @test c[4] ≈ 0.0 atol=1e-14   # Rʸˣ
            @test c[5] ≈ sh  atol=1e-14   # Rᵗˣ
            @test c[6] ≈ 0.0 atol=1e-14   # Rᵗʸ
            @test c[7] ≈ 0.0 atol=1e-14   # Rᵗᶻ
            @test c[8] ≈ 0.0 atol=1e-14   # Rᵗˣʸᶻ

            c = ga_components(Boost(η, [0.0, 1.0, 0.0]))
            @test c[1] ≈ ch  atol=1e-14   # R¹
            @test c[2] ≈ 0.0 atol=1e-14   # Rᶻʸ
            @test c[3] ≈ 0.0 atol=1e-14   # Rˣᶻ
            @test c[4] ≈ 0.0 atol=1e-14   # Rʸˣ
            @test c[5] ≈ 0.0 atol=1e-14   # Rᵗˣ
            @test c[6] ≈ sh  atol=1e-14   # Rᵗʸ
            @test c[7] ≈ 0.0 atol=1e-14   # Rᵗᶻ
            @test c[8] ≈ 0.0 atol=1e-14   # Rᵗˣʸᶻ

            c = ga_components(Boost(η, [0.0, 0.0, 1.0]))
            @test c[1] ≈ ch  atol=1e-14   # R¹
            @test c[2] ≈ 0.0 atol=1e-14   # Rᶻʸ
            @test c[3] ≈ 0.0 atol=1e-14   # Rˣᶻ
            @test c[4] ≈ 0.0 atol=1e-14   # Rʸˣ
            @test c[5] ≈ 0.0 atol=1e-14   # Rᵗˣ
            @test c[6] ≈ 0.0 atol=1e-14   # Rᵗʸ
            @test c[7] ≈ sh  atol=1e-14   # Rᵗᶻ
            @test c[8] ≈ 0.0 atol=1e-14   # Rᵗˣʸᶻ
        end
    end

    @testset "Boost: velocity-vector constructors" begin
        # QuatVec form matches explicit rapidity+direction form
        for (vx, vy, vz) ∈ [(0.5, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.0, 0.4), (0.2, 0.2, 0.2)]
            v⃗ = QuatVec(vx, vy, vz)
            β = abs(v⃗)
            @test Boost(v⃗) ≈ Boost(atanh(β), v⃗ / β)
            @test Boost([vx, vy, vz]) == Boost(v⃗)
        end
        # zero velocity → identity
        @test Boost(QuatVec(0.0, 0.0, 0.0)) == one(Lorentz{Float64})
    end

    @testset "Boost: known action on 4-vectors" begin
        for η ∈ [0.3, 0.7, 1.2, 1.8, 2.5]
            Λ = Boost(η, [0.0, 0.0, 1.0])
            ch, sh = cosh(η), sinh(η)

            @test Λ([1.0, 0.0, 0.0, 0.0]) ≈ [ch,  0.0, 0.0, sh]  atol=1e-13
            @test Λ([0.0, 0.0, 0.0, 1.0]) ≈ [sh,  0.0, 0.0, ch]  atol=1e-13
            @test Λ([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol=1e-13

            eη  = exp(η)
            eη⁻¹ = exp(-η)
            @test Λ([1.0, 0.0, 0.0,  1.0]) ≈ eη   .* [1.0, 0.0, 0.0,  1.0] atol=1e-13
            @test Λ([1.0, 0.0, 0.0, -1.0]) ≈ eη⁻¹ .* [1.0, 0.0, 0.0, -1.0] atol=1e-13
        end
    end

    @testset "Rotation and boost: known action on unit 4-vectors" begin
        # Rotation by π/3 about each coordinate axis.
        # Rotor(cos(θ/2), sin(θ/2), 0, 0) rotates about x (fixes [0,1,0,0])
        # Rotor(cos(θ/2), 0, sin(θ/2), 0) rotates about y (fixes [0,0,1,0])
        # Rotor(cos(θ/2), 0, 0, sin(θ/2)) rotates about z (fixes [0,0,0,1])
        θ = π / 3
        c, s = cos(θ), sin(θ)   # c = 1/2, s = √3/2

        Λ_Rx = Lorentz{Float64}(Rotor(cos(θ/2), sin(θ/2), 0.0, 0.0))
        Λ_Ry = Lorentz{Float64}(Rotor(cos(θ/2), 0.0, sin(θ/2), 0.0))
        Λ_Rz = Lorentz{Float64}(Rotor(cos(θ/2), 0.0, 0.0, sin(θ/2)))

        # Rotation about x: fixes t and x; rotates y→z plane
        @test Λ_Rx([1.0, 0.0, 0.0, 0.0]) ≈ [1.0, 0.0, 0.0, 0.0] atol=1e-14
        @test Λ_Rx([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol=1e-14
        @test Λ_Rx([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, 0.0,  c,  s]   atol=1e-14
        @test Λ_Rx([0.0, 0.0, 0.0, 1.0]) ≈ [0.0, 0.0, -s,  c]   atol=1e-14

        # Rotation about y: fixes t and y; rotates z→x plane
        @test Λ_Ry([1.0, 0.0, 0.0, 0.0]) ≈ [1.0, 0.0, 0.0, 0.0] atol=1e-14
        @test Λ_Ry([0.0, 1.0, 0.0, 0.0]) ≈ [0.0,  c, 0.0, -s]   atol=1e-14
        @test Λ_Ry([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, 0.0, 1.0, 0.0] atol=1e-14
        @test Λ_Ry([0.0, 0.0, 0.0, 1.0]) ≈ [0.0,  s, 0.0,  c]   atol=1e-14

        # Rotation about z: fixes t and z; rotates x→y plane
        @test Λ_Rz([1.0, 0.0, 0.0, 0.0]) ≈ [1.0, 0.0, 0.0, 0.0] atol=1e-14
        @test Λ_Rz([0.0, 1.0, 0.0, 0.0]) ≈ [0.0,  c,  s, 0.0]   atol=1e-14
        @test Λ_Rz([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, -s,  c, 0.0]   atol=1e-14
        @test Λ_Rz([0.0, 0.0, 0.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0] atol=1e-14

        # Boost with β = 1/3 (rapidity η = atanh(1/3)) in each direction.
        η  = atanh(1.0 / 3.0)
        ch, sh = cosh(η), sinh(η)

        Λ_Bx = Boost(η, [1.0, 0.0, 0.0])
        Λ_By = Boost(η, [0.0, 1.0, 0.0])
        Λ_Bz = Boost(η, [0.0, 0.0, 1.0])

        # Boost in x: mixes t and x; fixes y and z
        @test Λ_Bx([1.0, 0.0, 0.0, 0.0]) ≈ [ch,  sh,  0.0, 0.0] atol=1e-14
        @test Λ_Bx([0.0, 1.0, 0.0, 0.0]) ≈ [sh,  ch,  0.0, 0.0] atol=1e-14
        @test Λ_Bx([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, 0.0, 1.0, 0.0] atol=1e-14
        @test Λ_Bx([0.0, 0.0, 0.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0] atol=1e-14

        # Boost in y: mixes t and y; fixes x and z
        @test Λ_By([1.0, 0.0, 0.0, 0.0]) ≈ [ch,  0.0, sh,  0.0] atol=1e-14
        @test Λ_By([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol=1e-14
        @test Λ_By([0.0, 0.0, 1.0, 0.0]) ≈ [sh,  0.0, ch,  0.0] atol=1e-14
        @test Λ_By([0.0, 0.0, 0.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0] atol=1e-14

        # Boost in z: mixes t and z; fixes x and y
        @test Λ_Bz([1.0, 0.0, 0.0, 0.0]) ≈ [ch,  0.0, 0.0, sh]  atol=1e-14
        @test Λ_Bz([0.0, 1.0, 0.0, 0.0]) ≈ [0.0, 1.0, 0.0, 0.0] atol=1e-14
        @test Λ_Bz([0.0, 0.0, 1.0, 0.0]) ≈ [0.0, 0.0, 1.0, 0.0] atol=1e-14
        @test Λ_Bz([0.0, 0.0, 0.0, 1.0]) ≈ [sh,  0.0, 0.0, ch]  atol=1e-14
    end

    @testset "Boost: collinear rapidities add" begin
        for (η₁, η₂) ∈ [(0.3, 0.5), (1.0, 1.5), (0.1, 2.0), (0.7, 0.7)]
            for n̂ ∈ [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
                Λ_prod = Boost(η₁, n̂) * Boost(η₂, n̂)
                Λ_sum  = Boost(η₁ + η₂, n̂)
                for v ∈ [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.5, -0.3], [1.0, 0.2, 0.3, 0.9]]
                    @test norm(Λ_prod(v) - Λ_sum(v)) ≤ 1e-12
                end
            end
        end
    end

    # ── GA norm conditions ─────────────────────────────────────────────────────

    @testset "GA norm conditions hold for pure transformations" begin
        for Λ ∈ lg_rot_Λs
            @test ga_norm_conditions(Λ)
        end
        for Λ ∈ lg_boost_Λs
            @test ga_norm_conditions(Λ)
        end
    end

    @testset "GA norm conditions preserved under composition" begin
        for Λ ∈ lg_composed
            @test ga_norm_conditions(Λ)
        end
    end

    @testset "GA reverse is the group inverse (component test)" begin
        for Λ ∈ lg_rot_Λs
            @test ga_reverse_components(Λ)
        end
        for Λ ∈ lg_boost_Λs
            @test ga_reverse_components(Λ)
        end
        for Λ ∈ lg_composed
            @test ga_reverse_components(Λ)
        end
    end

    # ── Group structure: boosts ────────────────────────────────────────────────

    @testset "Boost: group properties" begin
        for Λ ∈ lg_boost_Λs, v ∈ lg_gen_vecs
            @test identity_law(Λ, v)
        end
        for i ∈ 1:(lg_n-1), v ∈ lg_gen_vecs
            @test composition_consistent(lg_boost_Λs[i], lg_boost_Λs[i+1], v)
        end
        for Λ ∈ lg_boost_Λs, v ∈ lg_gen_vecs
            @test inverse_law(Λ, v)
        end
        for i ∈ 1:(lg_n-2), v ∈ lg_gen_vecs
            @test associativity_law(lg_boost_Λs[i], lg_boost_Λs[i+1], lg_boost_Λs[i+2], v)
        end
        for Λ ∈ lg_boost_Λs, i ∈ 1:(lg_n-1)
            @test metric_preserved(Λ, lg_gen_vecs[i], lg_gen_vecs[i+1])
        end
    end

end  # @testset "Lorentz group"
