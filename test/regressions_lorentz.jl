# Regression tests for fixes to src/Lorentz.jl.
#
# Reference values are computed independently of the package, from explicit 4×4 matrices or
# closed-form expressions evaluated in BigFloat.

@testmodule LorentzReference begin
    using LinearAlgebra

    # The 4×4 matrix of an active boost with rapidity η along the unit direction n
    function boost_matrix(η, n)
        ch, sh = cosh(η), sinh(η)
        [ch  sh*n'; sh*n  I + (ch - 1)*n*n']
    end

    # The derivative of `boost_matrix` with respect to η
    function boost_matrix_derivative(η, n)
        ch, sh = cosh(η), sinh(η)
        [sh  ch*n'; ch*n  sh*n*n']
    end

    # The 4×4 matrix of the rotation given by the unit quaternion with components (w, x, y, z)
    function rotation_matrix(w, x, y, z)
        R = [
            1-2(y^2+z^2)  2(x*y-w*z)    2(x*z+w*y);
            2(x*y+w*z)    1-2(x^2+z^2)  2(y*z-w*x);
            2(x*z-w*y)    2(y*z+w*x)    1-2(x^2+y^2)
        ]
        [one(w) zeros(typeof(w), 1, 3); zeros(typeof(w), 3, 1) R]
    end

    # The 4×4 matrix of the action of Λ on 4-vectors, column by column
    function action_matrix(Λ, ::Type{T}) where {T}
        e(i) = T[j == i for j ∈ 1:4]
        reduce(hcat, [Λ(e(i)) for i ∈ 1:4])
    end

    # The components (cosh(η/2), sinh(η/2) v⃗/β) of `Boost(v⃗)`, from the rapidity η = atanh(β)
    function boost_from_velocity(v)
        β = norm(v)
        η = atanh(β)
        shc = iszero(β) ? one(β)/2 : sinh(η/2)/β
        [cosh(η/2); shc .* v]
    end
end

@testitem "Lorentz regressions: the action promotes the element type" tags=[:unit, :fast] setup=[LorentzReference] begin
    using StaticArrays, ForwardDiff, LinearAlgebra
    using .LorentzReference: boost_matrix, boost_matrix_derivative

    n = [0.6, -0.48, 0.64]
    Λ = Boost(0.5, n)
    M = boost_matrix(big(0.5), big.(n))

    v = [1, 0, 0, 0]
    @test Λ(v) isa Vector{Float64}
    @test Λ(v) ≈ Float64.(M * v) rtol=4eps()
    @test Λ(Float32[1, 0, 0, 0]) isa Vector{Float64}
    @test Λ(Float32[1, 0, 0, 0]) ≈ Float64.(M * v) rtol=4eps()
    w = view([0.3, 1.2, -0.5, 0.7, 9.0], 1:4)
    @test Λ(w) isa Vector{Float64}
    @test Λ(w) ≈ Float64.(M * big.(w)) rtol=4eps()
    @test Λ(SA[1, 0, 0, 0]) isa SVector{4, Float64}
    @test Λ(SA[1, 0, 0, 0]) ≈ Λ(v) rtol=4eps()
    @test Λ(1:4) isa Vector{Float64}
    @test Λ(1:4) ≈ Float64.(M * (1:4)) rtol=4eps()
    @test Λ(MVector(1, 0, 0, 0)) isa MVector{4, Float64}
    @test Λ(MVector(1, 0, 0, 0)) ≈ Float64.(M * v) rtol=4eps()

    # A BigFloat transformation is not narrowed to the element type of the vector
    nB = normalize(big.(n))
    ΛB = Boost(big(0.5), nB)
    MB = boost_matrix(big(0.5), nB)
    @test ΛB([0.3, 1.2, -0.5, 0.7]) isa Vector{BigFloat}
    @test ΛB([0.3, 1.2, -0.5, 0.7]) ≈ MB * big.([0.3, 1.2, -0.5, 0.7]) rtol=100eps(BigFloat)

    @test_throws DimensionMismatch Λ([1.0, 0.0, 0.0])
    @test_throws DimensionMismatch Λ([1.0, 0.0, 0.0, 0.0, 0.0])

    # ForwardDiff through the transformation applied to an ordinary `Vector`
    v4 = [0.3, 1.2, -0.5, 0.7]
    for η ∈ (0.0, 0.4, 2.0)
        d = ForwardDiff.derivative(t -> Boost(t, n)(v4), η)
        @test d ≈ Float64.(boost_matrix_derivative(big(η), big.(n)) * big.(v4)) atol=8eps(cosh(η))
    end
end

@testitem "Lorentz regressions: Boost promotes the element type" tags=[:unit, :fast] begin
    using ForwardDiff, LinearAlgebra

    n = [0.6, -0.48, 0.64]
    @test Boost(1, [0, 0, 1]) isa Lorentz{Float64}
    @test components(Boost(1, [0, 0, 1])) ≈ components(Boost(1.0, [0.0, 0.0, 1.0]))
    @test Boost(0.7f0, n) isa Lorentz{Float64}
    @test Boost(0.7f0, QuatVec(n...)) isa Lorentz{Float64}
    @test Boost(0.7, Float32.(n)) isa Lorentz{Float64}
    @test Boost(0.7f0, Float32.(n)) isa Lorentz{Float32}
    @test Boost(0.7, big.(n)) isa Lorentz{BigFloat}

    # The result is as precise as its element type claims, so a Float32 rapidity does not
    # limit the precision of a Float64 result, and similarly for BigFloat.
    spinor_norm_defect(L) = abs(sum(components(L) .^ 2) - 1)
    @test spinor_norm_defect(Boost(0.7f0, n)) ≤ 4eps()
    @test spinor_norm_defect(Boost(0.7f0, QuatVec(n...))) ≤ 4eps()
    @test components(Boost(0.7f0, n)) ≈ components(Boost(Float64(0.7f0), n)) rtol=4eps()
    nB = normalize(big.(n))
    @test spinor_norm_defect(Boost(0.7, nB)) ≤ 10eps(BigFloat)
    @test components(Boost(0.7, nB)) ≈ components(Boost(big(0.7), nB)) rtol=10eps(BigFloat)
    @test components(Boost(0.7, QuatVec(nB...))) ≈ components(Boost(big(0.7), nB)) rtol=10eps(BigFloat)
    @test Boost([0, 0, 0]) isa Lorentz{Float64}
    @test Boost([1//3, 0, 0]) isa Lorentz{Float64}
    @test components(Boost([1//3, 0, 0])) ≈ components(Boost(atanh(1/3), [1.0, 0, 0]))

    # Differentiation with respect to the direction alone.  The components are
    # (cosh(η/2), im*sinh(η/2)*n), so this weighted sum has gradient sinh(η/2) .* [2, 3, 4].
    weighted(L) = sum(real.(components(L))) + sum(imag.(components(L)) .* (1:4))
    expected = sinh(0.25) .* [2, 3, 4]
    @test ForwardDiff.gradient(nn -> weighted(Boost(0.5, nn)), n) ≈ expected rtol=4eps()
    @test ForwardDiff.gradient(nn -> weighted(Boost(0.5, quatvec(nn...))), n) ≈ expected rtol=4eps()
end

@testitem "Lorentz regressions: Boost from a velocity" tags=[:unit, :fast] setup=[LorentzReference] begin
    using ForwardDiff, LinearAlgebra
    using .LorentzReference: boost_from_velocity

    function boost_components(v)
        c = components(Boost(v))
        [real(c[1]); imag.(c[2:4])]
    end

    dir = [0.6, -0.48, 0.64]
    setprecision(BigFloat, 256) do
        # Values, including velocities so small that β² underflows
        for s ∈ (0.0, 1e-300, 1e-170, 1e-16, 1e-8, 1e-3, 0.3, 0.9, 0.99)
            v = s .* dir
            b = boost_components(v)
            r = boost_from_velocity(big.(v))
            @test all(isfinite, b)
            @test maximum(abs.(b .- r) ./ abs.(r .+ (r .== 0))) ≤ 10eps() / (1 - s^2)
        end

        # Second derivatives near zero velocity, compared with central finite differences
        # of the BigFloat reference
        W = [0.3, -1.1, 0.7, 0.45]
        F(x) = W ⋅ boost_components(x)
        Fref(x) = W ⋅ boost_from_velocity(x)
        h = big"1e-25"
        function hessian_reference(x)
            e(i) = BigFloat[j == i ? h : 0 for j ∈ 1:3]
            [(Fref(x+e(i)+e(j)) - Fref(x+e(i)-e(j)) - Fref(x-e(i)+e(j)) + Fref(x-e(i)-e(j))) / 4h^2
             for i ∈ 1:3, j ∈ 1:3]
        end
        for s ∈ (0.0, 1e-12, 1e-8, 1e-5, 0.3)
            x = s .* dir
            # Offset the reference point slightly from zero to avoid its special case there
            Href = hessian_reference(big.(x) .+ (s == 0 ? big"1e-40" : big(0)))
            @test maximum(abs.(ForwardDiff.hessian(F, x) .- Href)) ≤ 2eps()
        end
    end

    # The QuatVec and AbstractVector forms agree with the rapidity form
    for v ∈ ([0.5, 0.0, 0.0], [0.2, 0.2, 0.2], [0.0, -0.3, 0.4])
        β = norm(v)
        @test components(Boost(QuatVec(v...))) == components(Boost(v))
        @test components(Boost(v)) ≈ components(Boost(atanh(β), v ./ β)) rtol=4eps()
    end
    @test Boost(QuatVec(0.0, 0.0, 0.0)) == one(Lorentz{Float64})
end

@testitem "Lorentz regressions: Boost input validation" tags=[:unit, :fast] begin
    @test_throws DomainError Boost([0.0, 0.0, 1.0])
    @test_throws DomainError Boost([0.0, 0.0, 1.5])
    @test_throws DomainError Boost(QuatVec(0.8, 0.6, 0.0))
    @test_throws DimensionMismatch Boost([0.1])
    @test_throws DimensionMismatch Boost([0.1, 0.2, 0.3, 0.4])
    @test_throws DimensionMismatch Boost(0.5, [0.0, 1.0])
    @test_throws DomainError Boost([NaN, 0.0, 0.0])
end

@testitem "Lorentz regressions: symbolic Boost" tags=[:unit] begin
    using Symbolics

    # The domain check of the velocity forms must not require a `Bool` from a symbolic
    # comparison.
    Symbolics.@variables a b c η
    @test Boost(QuatVec(a, b, c)) isa Lorentz{Symbolics.Num}
    @test Boost([a, b, c]) isa Lorentz{Symbolics.Num}
    @test Boost(η, [0.6, -0.48, 0.64]) isa Lorentz{Symbolics.Num}
    @test Boost(0.5, [a, b, c]) isa Lorentz{Symbolics.Num}

    # Substituting numbers into the symbolic result reproduces the numerical boost.
    v = [0.3, -0.2, 0.4]
    L = Boost([a, b, c])
    numeric(z) = Symbolics.build_function(z, a, b, c; expression=Val(false))(v...)
    w = components(L)
    @test [complex(numeric(real(z)), numeric(imag(z))) for z ∈ w] ≈ components(Boost(v)) rtol=4eps()
end

@testitem "Lorentz regressions: Zygote through the action" tags=[:unit] setup=[LorentzReference] begin
    using Zygote
    using .LorentzReference: boost_matrix, boost_matrix_derivative

    n = [0.6, -0.48, 0.64]
    v4 = [0.3, 1.2, -0.5, 0.7]
    weights = [1.0, -2.0, 3.0, 0.5]
    for η ∈ (0.0, 0.4, 2.0)
        M = boost_matrix(big(η), big.(n))
        dM = boost_matrix_derivative(big(η), big.(n))
        # Gradient with respect to the vector, for a `Vector` and for a view
        g = Zygote.gradient(v -> weights' * Boost(η, n)(v), v4)[1]
        @test g ≈ Float64.(M' * weights) rtol=8eps(cosh(η))
        gview = Zygote.gradient(v -> weights' * Boost(η, n)(view(v, 2:5)), [0.0; v4])[1]
        @test gview ≈ [0.0; Float64.(M' * weights)] rtol=8eps(cosh(η))
        # Derivative with respect to the rapidity
        d = Zygote.gradient(t -> weights' * Boost(t, n)(v4), η)[1]
        @test d ≈ Float64(weights' * (dM * big.(v4))) rtol=8eps(cosh(η))
    end
end

@testitem "Lorentz regressions: constructors" tags=[:unit, :fast] begin
    for w ∈ (1.0, -2.0, 1.0im, 3, 0.5f0)
        L = Lorentz(w)
        @test L isa Lorentz
        @test components(L) ≈ components(Lorentz(w, 0, 0, 0))
    end
    @test components(Lorentz(1.0)) == components(one(Lorentz{Float64}))
    @test components(Lorentz(-2.0)) == components(-one(Lorentz{Float64}))
    @test Lorentz(Float64) === Lorentz{Float64}
    @test Lorentz(Float32) === Lorentz{Float32}
    @test Lorentz(QuaternionF64) === Lorentz{Float64}
    @test components(Lorentz([2.0])) == components(one(Lorentz{Float64}))
    @test components(Lorentz([1.0, 0.0, 0.0])) == components(Lorentz(0.0, 1.0, 0.0, 0.0))
    @test components(Lorentz([1.0, 2, 3, 4])) ≈ components(Lorentz(rotor(1.0, 2, 3, 4)))
    @test_throws DimensionMismatch Lorentz([1.0, 2.0])
    # Normalization uses the complex spinor norm, so the bivector 𝐭𝐱 becomes the rotation 𝐢
    @test components(Lorentz(QuatVec(0, 1.0im, 0, 0))) ≈ components(Lorentz(imx))
    # ga_components is unaffected by the move to `components`
    @test ga_components(Boost(0.7, [0.0, 0.0, 1.0])) ≈ [cosh(0.35), 0, 0, 0, 0, 0, sinh(0.35), 0]
end

@testitem "Lorentz regressions: ℂ helpers" tags=[:unit, :fast] begin
    import Quaternionic: ℂconj, ℂreal, ℂimag, ℂreim, RB
    using LinearAlgebra

    q = Quaternion(1.0+2im, 3.0-1im, 0.5im, 2.0+0im)
    @test ℂconj(q) isa Quaternion{ComplexF64}
    @test components(ℂconj(q)) == conj.(components(q))
    @test ℂreim(q) == (ℂreal(q), ℂimag(q))
    @test ℂreal(q) == Quaternion(1.0, 3.0, 0.0, 2.0)
    @test ℂimag(q) == Quaternion(2.0, -1.0, 0.5, 0.0)

    # For Λ = R*B, ℂconj(Λ) = R*inv(B), and abs(ℂreal(Λ)) = cosh(η/2)
    R = rotor(0.4, 1.1, -0.7, 0.3)
    B = Boost(0.7, normalize([0.3, -0.5, 0.8]))
    Λ = Lorentz(R) * B
    @test ℂconj(Λ) isa Lorentz{Float64}
    @test components(ℂconj(Λ)) ≈ components(Lorentz(R) * inv(B)) atol=4eps()
    @test abs(ℂreal(Λ)) ≈ cosh(0.35) rtol=4eps()
end

@testitem "Lorentz regressions: RB and BR at large rapidity" tags=[:unit, :fast] setup=[LorentzReference] begin
    import Quaternionic: RB, BR, Rv
    using LinearAlgebra
    using .LorentzReference: action_matrix

    Rx = Rotor(cos(0.3), sin(0.3), 0, 0)
    n̂ = normalize([0.3, -0.5, 0.8])
    relerr(A, B) = maximum(abs.(components(A) .- components(B))) / maximum(abs.(components(B)))
    for η ∈ (2.0, 10.0, 20.0, 30.0)
        Bη = Boost(η, n̂)
        # Build the products without renormalization, so that they are accurate
        ΛRB = Lorentz{Float64}(Quaternion(Lorentz(Rx)) * Quaternion(Bη))
        ΛBR = Lorentz{Float64}(Quaternion(Bη) * Quaternion(Lorentz(Rx)))

        R, B = RB(ΛRB)
        @test B isa Lorentz{Float64}
        @test relerr(B, Bη) ≤ 4eps()
        @test min(relerr(R, Rx), relerr(-R, Rx)) ≤ 4eps()
        # In Λ = R * B, the boost acts first and then the rotation
        v = [0.3, 1.2, -0.5, 0.7]
        @test ΛRB(v) ≈ Lorentz(R)(B(v)) rtol=16eps()

        B2, R2 = BR(ΛBR)
        @test B2 isa Lorentz{Float64}
        @test relerr(B2, Bη) ≤ 4eps()
        @test min(relerr(R2, Rx), relerr(-R2, Rx)) ≤ 4eps()
        @test ΛBR(v) ≈ B2(Lorentz(R2)(v)) rtol=16eps()
    end
end

@testitem "Lorentz regressions: KAN return types" tags=[:unit, :fast] setup=[LorentzReference] begin
    import Quaternionic: KAN
    using LinearAlgebra
    using .LorentzReference: boost_matrix, action_matrix

    for T ∈ (Float32, Float64, BigFloat)
        Λ = Lorentz(rotor(T(0.4), T(1.1), T(-0.7), T(0.3))) * Boost(T(0.6), normalize(T[0.2, -0.3, 0.25]))
        Rₖ, Rₐ, Rₙ = KAN(Λ)
        @test Rₖ isa Rotor{T}
        @test Rₐ isa Lorentz{T}
        @test Rₙ isa Lorentz{T}
        # Rₐ is a boost along 𝐳 that can be used directly as a Lorentz transformation
        c = ga_components(Rₐ)
        φ = 2asinh(c[7])
        @test c[[2, 3, 4, 5, 6, 8]] ≈ zeros(T, 6) atol=10eps(T)
        @test action_matrix(Rₐ, T) ≈ boost_matrix(φ, T[0, 0, 1]) rtol=10eps(T)
        @test Lorentz(Rₖ) * Rₐ * Rₙ isa Lorentz{T}
        @test components(Lorentz(Rₖ) * Rₐ * Rₙ) ≈ components(Λ) rtol=20eps(T)
    end
end

@testitem "Lorentz regressions: composition at large rapidity" tags=[:unit, :fast] setup=[LorentzReference] begin
    using LinearAlgebra
    using .LorentzReference: boost_matrix, rotation_matrix, action_matrix

    # The action of a composition is compared with the product of explicit BigFloat
    # matrices.  The relative error is measured against the largest matrix entry, which is
    # of order cosh(η).  The spinor norm of a rounded rotor with large rapidity is
    # ill-conditioned, but the action is not, so it must remain accurate in this sense.
    relerr(L, M) = Float64(maximum(abs.(L .- M)) / maximum(abs.(M)))
    setprecision(BigFloat, 256) do
        n̂ = normalize([0.3, -0.5, 0.8])
        x̂, ŷ, ẑ = [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]
        for (η₁, n₁, η₂, n₂) ∈ (
            (10.0, n̂, 10.0, n̂),  # total rapidity 20
            (19.0, n̂, 19.0, n̂),  # total rapidity 38
            (30.0, n̂, 30.0, n̂),  # total rapidity 60
            (19.0, ẑ, 19.0, ẑ),
            (10.0, x̂, 10.0, ŷ),
        )
            Λ = Boost(η₁, n₁) * Boost(η₂, n₂)
            M = boost_matrix(big(η₁), big.(n₁)) * boost_matrix(big(η₂), big.(n₂))
            @test relerr(action_matrix(Λ, Float64), M) ≤ 20eps()
        end

        # A rotation composed with a large boost
        R = rotor(0.4, 1.1, -0.7, 0.3)
        Λ = Lorentz(R) * Boost(38.0, ẑ)
        M = rotation_matrix(big.(components(R))...) * boost_matrix(big(38.0), big.(ẑ))
        @test relerr(action_matrix(Λ, Float64), M) ≤ 20eps()

        # Float32, at a total rapidity of 16
        n̂32 = normalize(Float32[0.3, -0.5, 0.8])
        Λ = Boost(8f0, n̂32) * Boost(8f0, n̂32)
        M = boost_matrix(big(8.0), big.(n̂32)) * boost_matrix(big(8.0), big.(n̂32))
        @test relerr(action_matrix(Λ, Float32), M) ≤ 20eps(Float32)
    end
end
