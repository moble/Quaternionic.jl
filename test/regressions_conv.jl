# Regression tests for conversions, distances, and alignment.

@testitem "distance2 near-identity series matches BigFloat" tags=[:validation, :fast] begin
    using DoubleFloats: Double64
    import Quaternionic: absvec

    # Rotors separated from the identity by a rotation through θ are a distance θ/2 away.
    # The reference is computed in BigFloat from the very same rotor, so that the rounding
    # of its components does not enter.  The grid straddles the threshold at which
    # `distance2` switches from its series to the closed form, and a tolerance of a few
    # ulps detects an error in any coefficient of the series, or in the threshold.
    θs = exp10.(range(-10, log10(0.2), length=400))
    @testset "$T" for T ∈ (Float32, Float64, Double64)
        n̂ = normalize(quatvec(T(1), T(-2), T(3)))
        worst = maximum(θs) do θ
            R = exp(T(θ) * n̂ / 2)
            d = distance2(R, one(R))
            Rb = Rotor{BigFloat}(R)
            ref = atan(absvec(Rb), Rb[1])^2
            Float64(abs(d - ref) / (eps(T) * ref))
        end
        @test worst ≤ 4
        # The same for `distance`, which is the square root
        worst_sqrt = maximum(θs) do θ
            R = exp(T(θ) * n̂ / 2)
            d = distance(one(R), R)
            Rb = Rotor{BigFloat}(R)
            ref = atan(absvec(Rb), Rb[1])
            Float64(abs(d - ref) / (eps(T) * ref))
        end
        @test worst_sqrt ≤ 4
        # Away from the identity, the distance is half the angle, up to the rounding
        # incurred by forming the two rotors
        R₀ = rotor(T(0.4), T(-0.2), T(0.7), T(0.5))
        @test all(θs) do θ
            abs(distance(exp(T(θ) * n̂ / 2) * R₀, R₀) - T(θ) / 2) ≤ 10eps(T)
        end
    end
end

@testitem "distance2 on symbolic and Lorentz rotors" tags=[:unit, :fast] begin
    import Symbolics

    # For a symbolic rotor, the comparison that selects the series does not evaluate to a
    # `Bool`, so the closed form is used, and it must still evaluate symbolically.
    Symbolics.@variables w x y z a b c d
    Q₁ = Rotor{Symbolics.Num}(w, x, y, z)
    Q₂ = Rotor{Symbolics.Num}(a, b, c, d)
    expr = distance2(Q₁, Q₂)
    @test expr isa Symbolics.Num
    R₁ = rotor(0.4, -0.2, 0.7, 0.5)
    R₂ = rotor(-0.1, 0.3, 0.2, 0.9)
    f = Symbolics.build_function(expr, [w, x, y, z, a, b, c, d]; expression=Val(false))
    @test f([components(R₁)...; components(R₂)...]) ≈ distance2(R₁, R₂) rtol=1e-14

    # The Lorentz group has no positive-definite bi-invariant metric
    L = rotor(1.0 + 0.2im, 0.3, 0.1im, 0.2)
    @test_throws ArgumentError distance2(L, L)
    @test_throws ArgumentError distance(L, rotor(1.0, 0, 0, 0))
end

@testitem "distance2 gradients and the series branch" tags=[:unit, :fast] begin
    using ForwardDiff

    # Dual numbers must take the series branch at the identity, where the closed form has
    # a NaN derivative, and agree with plain floats elsewhere
    R₀ = rotor(0.4, -0.2, 0.7, 0.5)
    g = ForwardDiff.gradient(u -> distance2(rotor(u...), R₀), collect(components(R₀)))
    @test g ≈ zeros(4) atol=1e-14
    @test !any(isnan, g)
    # A point where the scalar part of the quotient is zero must not divide by zero
    v = ForwardDiff.gradient(u -> distance2(rotor(u...), rotor(1.0, 0, 0, 0)), [0.0, 1.0, 0.0, 0.0])
    @test !any(isnan, v)
end

@testitem "distance2 Hessians under ReverseDiff take the series" tags=[:unit, :fast] begin
    using ForwardDiff, ReverseDiff

    # `ReverseDiff.TrackedReal` is not an `AbstractFloat`, but it must still take the
    # series at the identity.  The closed form has a singular second derivative there, so
    # the ReverseDiff Hessian of the closed form was wrong by O(1), with no NaN.
    R₀ = rotor(0.4, -0.2, 0.7, 0.5)
    f(u) = distance2(rotor(u...), R₀)
    for u ∈ (collect(components(R₀)), collect(components(rotor(0.3, 0.5, -0.1, 0.8))))
        Hr = ReverseDiff.hessian(f, u)
        Hf = ForwardDiff.hessian(f, u)
        @test !any(isnan, Hr)
        @test Hr ≈ Hf atol=1e-14
    end
    # The Hessian at the minimum, compared with the explicit formula.  At the identity,
    # `distance2(rotor(u...), R₀)` is the squared length of the vector part of the
    # normalized `u/R₀`, to second order, so its Hessian is twice the projector onto the
    # tangent space spanned by `imx * R₀`, `imy * R₀`, and `imz * R₀`.
    u₀ = collect(components(R₀))
    P = sum(e -> (v = collect(components(e * R₀)); v * v'), (imx, imy, imz))
    @test ReverseDiff.hessian(f, u₀) ≈ 2P atol=1e-14
end

@testitem "to_spherical_coordinates is accurate near the poles" tags=[:validation, :fast] begin
    using ForwardDiff

    # The polar angle was computed with acos, which lost half the digits near θ = 0
    @testset "$T" for T ∈ (Float16, Float32, Float64)
        θs = T === Float16 ? T[1e-3, 1e-2, 0.1] : T[1e-12, 1e-8, 1e-6, 1e-3, 0.1]
        for θ ∈ θs, ϕ ∈ T[0.7, -2.5]
            R = from_spherical_coordinates(θ, ϕ)
            Rb = Rotor{BigFloat}(R)
            θref = 2atan(hypot(Rb[2], Rb[3]), hypot(Rb[1], Rb[4]))
            @test abs(to_spherical_coordinates(R)[1] - θref) ≤ 2eps(T) * θref
        end
    end
    # The near-pole rotor found by the property test of round trips
    R = rotor(-0.7074310290501881, -0.007599356757375066, -0.0013817012286136486, -0.7067401784358812)
    θ, ϕ = to_spherical_coordinates(R)
    S = from_spherical_coordinates(θ, ϕ)
    @test maximum(abs, to_rotation_matrix(S)[:, 3] - to_rotation_matrix(R)[:, 3]) ≤ 10eps(Float64)
    # The gradient is finite and matches the exact value at a small polar angle
    s = 1e-9
    F(u) = to_spherical_coordinates(rotor(u...))[1]
    u = [1.0, s, -2s, 0.5s]
    @test F(u) ≈ 2atan(√5 * s, hypot(1.0, 0.5s)) rtol=4eps()
    g = ForwardDiff.gradient(F, u)
    @test all(isfinite, g)
    @test g[2:3] ≈ [2, -4] / √5 rtol=1e-6
    # The normalization scales out, even where the sum of squares would underflow
    @test to_spherical_coordinates(1e-170 * rotor(0.4, -0.2, 0.7, 0.5)) ≈
        to_spherical_coordinates(rotor(0.4, -0.2, 0.7, 0.5))
end

@testitem "to_euler_phases handles integers" tags=[:validation, :fast] begin
    q = quaternion(0.4, -0.2, 0.7, 0.5)
    z = to_euler_phases(q)
    # Integer input used to throw an `InexactError`
    @test to_euler_phases(quaternion(1, 2, 3, 4)) == to_euler_phases(quaternion(1.0, 2, 3, 4))
    @test to_euler_phases(rotor(1, 0, 0, 0)) == to_euler_phases(rotor(1.0, 0, 0, 0))
    # The normalization scales out
    for s ∈ (1e-100, 1e-10, 1e10, 1e100)
        @test to_euler_phases(s * q) ≈ z atol=4eps()
    end
    # Agreement with a BigFloat evaluation
    for _ ∈ 1:100
        R = randn(RotorF64)
        @test to_euler_phases(R) ≈ to_euler_phases(Rotor{BigFloat}(R)) atol=4eps()
    end
end

@testitem "from_euler_phases has no NaN gradients under ReverseDiff" tags=[:unit, :fast] begin
    using ForwardDiff, ReverseDiff

    # The sign test used `abs` of a complex number that is often exactly zero, whose
    # derivative ReverseDiff recorded as NaN
    wsum(R) = sum((1:4) .* components(R))
    f(x) = wsum(from_euler_phases(cis(x[1]), cis(x[2]), cis(x[3])))
    for x ∈ ([0.3, 1.1, -0.4], [0.3, 0.5, -0.1], [1.0, 0.5, 2.0])
        gr = ReverseDiff.gradient(f, x)
        gf = ForwardDiff.gradient(f, x)
        @test !any(isnan, gr)
        @test gr ≈ gf atol=1e-14
    end
end

@testitem "conversions keep their values after leaving the literal macros" tags=[:unit, :fast] begin
    using StaticArrays

    # Explicit formulas, independent of the package's implementation
    q = quaternion(0.3, 0.5, -0.1, 0.8)
    n = abs2(q)
    w, x, y, z = components(q)
    ℛ = [
        w^2+x^2-y^2-z^2  2(x*y-w*z)  2(x*z+w*y);
        2(x*y+w*z)  w^2-x^2+y^2-z^2  2(y*z-w*x);
        2(x*z-w*y)  2(y*z+w*x)  w^2-x^2-y^2+z^2
    ] / n
    M = to_rotation_matrix(q)
    @test M isa SMatrix{3,3,Float64}
    @test M ≈ ℛ atol=4eps()
    v = randn(QuatVecF64)
    @test M * vec(v) ≈ vec(rotor(q)(v)) atol=8eps()
    @test to_euler_angles(q) isa SVector{3,Float64}
    @test to_euler_phases(q) isa SVector{3,ComplexF64}
    @test to_spherical_coordinates(q) isa SVector{2,Float64}
end

@testitem "to_euler_angles round trips at gimbal lock" tags=[:unit, :fast] begin
    for α ∈ (-2.0, 0.3, 1.7), γ ∈ (-0.4, 0.9), β ∈ (0.0, Float64(π))
        R = from_euler_angles(α, β, γ)
        R′ = from_euler_angles(to_euler_angles(R))
        @test distance(R, R′) ≤ 10eps()
    end
end

@testitem "conversions reject complex input with an ArgumentError" tags=[:unit, :fast] begin
    using LinearAlgebra

    L = rotor(1.0 + 0.1im, 0.2, 0.3im, 0.1)
    @test_throws ArgumentError to_euler_angles(L)
    @test_throws ArgumentError to_euler_phases(L)
    @test_throws ArgumentError to_spherical_coordinates(L)
    @test_throws ArgumentError from_rotation_matrix(ComplexF64[1 0 0; 0 1 0; 0 0 1])
    a⃗ = [quatvec(1.0 + 0im, 0, 0), quatvec(0, 1.0 + 0im, 0)]
    b⃗ = [quatvec(0, 1.0 + 0im, 0), quatvec(1.0 + 0im, 0, 0)]
    @test_throws ArgumentError align(a⃗, b⃗)
    # `to_rotation_matrix` of a Lorentz rotor is its complex-orthogonal action on bivectors
    Mx = to_rotation_matrix(L)
    @test transpose(Mx) * Mx ≈ I atol=1e-14
end

@testitem "from_rotation_matrix and align: types and sign convention" tags=[:validation, :fast] begin
    using DoubleFloats: Double64
    using ForwardDiff, ReverseDiff

    # The result has the same element type as the input, including Float16
    for T ∈ (Float16, Float32, Float64, Double64, BigFloat)
        R = randn(Rotor{T})
        @test from_rotation_matrix(to_rotation_matrix(R)) isa Rotor{T}
        a⃗ = randn(QuatVec{T}, 5)
        @test align(a⃗, conj(R).(a⃗)) isa Rotor{T}
    end

    # The scalar part of the result is positive, whatever the element type, so that the
    # LAPACK and generic backends agree
    for _ ∈ 1:100
        R = randn(RotorF64)
        M = to_rotation_matrix(R)
        R₁ = from_rotation_matrix(M)
        R₂ = from_rotation_matrix(Double64.(M))
        R₃ = from_rotation_matrix(BigFloat.(M))
        @test R₁[1] > 0 && R₂[1] > 0 && R₃[1] > 0
        @test abs(R₁ - (R[1] < 0 ? -R : R)) ≤ 50eps()
        @test abs(R₁ - RotorF64(R₂)) ≤ 50eps()
        @test abs(R₁ - RotorF64(R₃)) ≤ 50eps()
        a⃗ = randn(QuatVecF64, 6)
        b⃗ = randn(QuatVecF64, 6)
        A₁ = align(a⃗, b⃗)
        A₂ = align(QuatVec{Double64}.(a⃗), QuatVec{Double64}.(b⃗))
        @test A₁[1] > 0 && A₂[1] > 0
        @test abs(A₁ - RotorF64(A₂)) ≤ 1e-12
    end

    # Rotations through π have zero scalar part, and the first nonzero component is positive
    for R ∈ (rotor(imx), rotor(imy), rotor(imz), rotor(0, 1, 1, 0), rotor(0, -1, 1, 0), rotor(0, 0, -1, 1))
        R′ = from_rotation_matrix(to_rotation_matrix(R))
        @test abs(R′[1]) ≤ 10eps()
        i = findfirst(c -> abs(c) > 10eps(), components(R′))
        @test R′[i] > 0
        @test min(abs(R′ - R), abs(R′ + R)) ≤ 10eps()
    end

    # A slightly non-orthogonal matrix gives the nearby rotor
    R = randn(RotorF64)
    M = to_rotation_matrix(R) + 1e-9 * randn(3, 3)
    R′ = from_rotation_matrix(M)
    @test min(abs(R′ - R), abs(R′ + R)) ≤ 1e-8

    # Because the sign is fixed, ReverseDiff and ForwardDiff differentiate the same branch
    wsum(R) = sum((1:4) .* components(R))
    f(x) = wsum(from_rotation_matrix(to_rotation_matrix(rotor(x...))))
    for x ∈ (normalize([0.3, -0.7, 0.5, 1.1]), normalize([-0.2, 0.4, 0.1, -0.9]))
        @test ReverseDiff.gradient(f, x) ≈ ForwardDiff.gradient(f, x) atol=1e-12
    end
end

@testitem "align on QuatVecs: docstring matrix, direction, and allocations" tags=[:validation, :fast] begin
    using LinearAlgebra

    a⃗ = randn(QuatVecF64, 6)
    b⃗ = randn(QuatVecF64, 6)
    R = align(a⃗, b⃗)
    # The matrix given in the docstring, whose dominant eigenvector is (w, x, y, z)
    S = sum(vec(a⃗[i]) * vec(b⃗[i])' for i in eachindex(a⃗, b⃗))
    s⃗ = [S[2,3]-S[3,2], S[3,1]-S[1,3], S[1,2]-S[2,1]]
    M = [tr(S) -s⃗'; -s⃗ (S + S' - tr(S) * I)]
    v = eigen(Symmetric(M)).vectors[:, end]
    @test min(abs(rotor(v...) - R), abs(rotor(v...) + R)) ≤ 1e-14

    # The result rotates b⃗ onto a⃗
    R₀ = randn(RotorF64)
    b⃗′ = conj(R₀).(a⃗)
    R′ = align(a⃗, b⃗′)
    @test maximum(abs, R′.(b⃗′) - a⃗) ≤ 1e-14

    # Forming S no longer allocates for each point
    function allocations(N)
        a⃗, b⃗, w = randn(QuatVecF64, N), randn(QuatVecF64, N), rand(N)
        align(a⃗, b⃗, w)
        align(a⃗, b⃗)
        (@allocated(align(a⃗, b⃗, w)), @allocated(align(a⃗, b⃗)))
    end
    small, large = allocations(10), allocations(1000)
    @test large[1] ≤ small[1] + 256
    @test large[2] ≤ small[2] + 256
end

@testitem "to_float_array of a single quaternion" tags=[:unit, :fast] begin
    f = to_float_array(quaternion(1, 2, 3, 4))
    @test f isa Vector{Float64}
    @test f == [1.0, 2.0, 3.0, 4.0]
    @test to_float_array(imx) == [0.0, 1.0, 0.0, 0.0]
    @test to_float_array(rotor(1.0f0, 0, 0, 0)) isa Vector{Float32}
end

@testitem "Zygote differentiates the conversions to static arrays" tags=[:unit, :fast] begin
    using ForwardDiff, Zygote

    # Zygote cannot differentiate the `@SVector [...]` and `@SMatrix [...]` literal forms,
    # which these functions used to build their results
    u = [0.3, 0.5, -0.1, 0.8]
    fs = (
        m -> to_rotation_matrix(quaternion(m...))[1, 2],
        m -> to_rotation_matrix(quaternion(m...))[3, 1],
        m -> to_euler_angles(quaternion(m...))[1],
        m -> to_euler_angles(quaternion(m...))[2],
        m -> real(to_euler_phases(quaternion(m...))[1]),
        m -> imag(to_euler_phases(quaternion(m...))[2]),
        m -> to_spherical_coordinates(quaternion(m...))[1],
        m -> to_spherical_coordinates(quaternion(m...))[2],
    )
    for f ∈ fs
        @test Zygote.gradient(f, u)[1] ≈ ForwardDiff.gradient(f, u) atol=1e-14
    end
end
