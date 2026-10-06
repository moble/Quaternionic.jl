# Smoke tests for the oldest Julia versions allowed by the compat bounds.
#
# The workspace test environment (test/Project.toml) cannot be instantiated on Julia
# versions before 1.9, so these checks run in a fresh environment instead.  The CI job that
# uses this file runs, from the root of the package,
#
#     julia -e 'using Pkg; Pkg.develop(path="."); Pkg.add(["ChainRulesCore", "ForwardDiff", "StaticArrays"])'
#     julia test/compat16_smoke.jl
#
# The file is not included by test/runtests.jl, and it contains no test items.

using Test
using Quaternionic
import Quaternionic: squad!
using ChainRulesCore, ForwardDiff, StaticArrays, LinearAlgebra, Random

@testset "Julia $VERSION smoke tests" begin
    Random.seed!(1)
    q = quaternion(1.0, 2.0, 3.0, 4.0)
    p = quaternion(-0.5, 0.25, 1.5, -2.0)
    R = rotor(1, 2, 3, 4)
    v = quatvec(1.0, -2.0, 0.5)
    Rs = [rotor(1, 0.1i, 0.2, -0.1i) for i ∈ 1:5]
    t = collect(1.0:5.0)

    @testset "Algebra" begin
        @test q * inv(q) ≈ 1
        @test q / q == 1
        @test (q / p) * p ≈ q
        @test quaternion(2e200, 0, 0, 0) / quaternion(1e200, 1e200, 0, 0) ≈ quaternion(1.0, -1.0, 0, 0)
        @test typeof(R * 2.0) <: Quaternion
        @test typeof(R + 1) <: Quaternion
        @test typeof(R * (1 + 2im)) <: Quaternion
        @test iszero(quatvec(quaternion(5, 1, 2, 3))[1])
        @test iszero(QuatVec(5.0, 1.0, 2.0, 3.0)[1])
        @test iszero(convert(QuatVecF64, quaternion(5.0, 1, 2, 3))[1])
        @test (v ⋅ v) ≈ abs2vec(v)
        @test v × v == quatvec(0, 0, 0)
        @test imx^2 == -1
        @test typeof(imx^3) == typeof(imx^2)
        @test zero(QuatVecF64) == 0
        @test !isequal(quaternion(5, 1, 2, 3), quatvec(1, 2, 3))
    end

    @testset "Math functions" begin
        @test exp(log(q)) ≈ q
        @test sqrt(q)^2 ≈ q
        @test sqrt(p)^2 ≈ p
        @test sqrt(quaternion(1e300, 1e300, 0, 0))^2 ≈ quaternion(1e300, 1e300, 0, 0)
        @test sqrt(quaternion(3e-310, 1e-310, 0, 0))^2 ≈ quaternion(3e-310, 1e-310, 0, 0)
        @test log(quaternion(1e300, 1e300, 0, 0)) ≈ quaternion(log(sqrt(2) * 1e300), π / 4, 0, 0)
        @test R^0.5 * R^0.5 ≈ R
        @test R^3 ≈ R * R * R
        @test R^0.5 isa Rotor
        @test exp(big(1) * imx) ≈ cos(big(1)) + sin(big(1)) * imx
    end

    @testset "Conversions" begin
        @test from_rotation_matrix(to_rotation_matrix(R)) ≈ R
        @test from_euler_angles(to_euler_angles(R)) ≈ R
        @test from_euler_phases(to_euler_phases(R)) ≈ R
        @test from_spherical_coordinates(to_spherical_coordinates(rotor(1, 0.2, -0.1, 0.3))) isa Rotor
        @test from_float_array(to_float_array([q, p])) == [q, p]
        @test sprint(show, quaternion(1, 2, 3, 4)) == "1 + 2𝐢 + 3𝐣 + 4𝐤"
        @test sprint(show, MIME("text/plain"), R) isa String
    end

    @testset "Distance, interpolation, and alignment" begin
        @test distance(R, R) == 0
        @test distance2(rotor(1, 0, 0, 0), rotor(cos(0.1), sin(0.1), 0, 0)) ≈ 0.01
        @test slerp(R, -R, 0.0) ≈ R
        @test all(squad(Rs, t, t) .≈ Rs)
        out = similar(Rs, 3)
        squad!(out, Rs, t, [1.5, 2.5, 3.5])
        @test out ≈ squad(Rs, t, [1.5, 2.5, 3.5])
        @test unflip([R, -R, R]) == [R, R, R]
        A = [R, -R, R]
        @test unflip!(A) === A
        s, ∂s1, ∂s2, ∂sτ = slerp∂slerp(R, rotor(2, 1, 0, 1), 0.3)
        @test s ≈ slerp(R, rotor(2, 1, 0, 1), 0.3)
        L, dL = log∂log(R)
        @test L ≈ log(R)
        E, dE = exp∂exp(quatvec(0.1, 0.2, 0.3))
        @test E ≈ exp(quatvec(0.1, 0.2, 0.3))
        @test align([quatvec(1.0, 0, 0), quatvec(0, 1.0, 0)], [quatvec(0, 1.0, 0), quatvec(-1.0, 0, 0)]) isa Rotor
    end

    @testset "ForwardDiff" begin
        @test ForwardDiff.derivative(τ -> components(exp(τ * imx))[2], 0.0) ≈ 1
        @test ForwardDiff.derivative(τ -> components(log(rotor(cos(τ), sin(τ), 0, 0)))[2], 0.0) ≈ 1
        f(τ) = components(sqrt(quaternion(1.0 + τ, 2, 3, 4)))[1]
        @test ForwardDiff.derivative(f, 0.1) ≈ (f(0.1 + 1e-6) - f(0.1 - 1e-6)) / 2e-6 rtol=1e-8
        g(x) = components(log(Quaternion(x...)))
        x = [0.3, -0.5, 0.7, 0.2]
        @test 1e160 * ForwardDiff.jacobian(g, 1e160 * x) ≈ ForwardDiff.jacobian(g, x)
    end

    @testset "ChainRulesCore" begin
        # Before Julia 1.9, Requires loads the extension when ChainRulesCore is loaded.
        y, pullback = rrule(abs2, q)
        @test y == abs2(q)
        @test unthunk(pullback(1.0)[2]) ≈ 2q
        y, pullback = rrule(distance2, R, rotor(2, 1, 0, 1))
        @test y ≈ distance2(R, rotor(2, 1, 0, 1))
        @test pullback(1.0)[2] isa Union{Rotor, Quaternion, AbstractThunk, Tangent}
    end

    @testset "Lorentz and random" begin
        B = Boost(quatvec(0.3, 0.0, 0.4))
        @test B([1.0, 0, 0, 0])[1] ≈ 1 / sqrt(1 - 0.25)
        @test Quaternionic.RB(B * Lorentz(R))[1] isa Rotor
        @test Lorentz(2.0) isa Rotor
        @test randn(QuaternionF64) isa QuaternionF64
        @test randn(RotorF32) isa RotorF32
        @test iszero(randn(QuatVecF64)[1])
    end
end
