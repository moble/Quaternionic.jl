# Regression tests for the math functions in src/math.jl: `abs`, `norm`, `inv`, `log`,
# `exp`, `sqrt`, `angle`, and `^`.

@testmodule MathRegressions begin
    using Quaternionic
    import ForwardDiff

    """
        rotorpower(R, s)

    Return the components of the unit rotor `R` raised to the real power `s`, computed in
    `BigFloat` directly from the axis-angle form of `R`.  This is independent of the
    package's own `log` and `^`.
    """
    function rotorpower(R, s)
        c = big.(components(R))
        c = c / sqrt(sum(abs2, c))
        vn = sqrt(c[2]^2 + c[3]^2 + c[4]^2)
        θ = atan(vn, c[1])
        n̂ = iszero(vn) ? [big(0), big(0), big(1)] : collect(c[2:4]) / vn
        [cos(big(s)*θ); sin(big(s)*θ) * n̂]
    end

    """
        maxerror(q, ref)

    Return the largest absolute difference between the components of `q` and `ref`.
    """
    maxerror(q, ref) = Float64(maximum(abs, big.(collect(components(q))) - big.(collect(ref))))

    """
        nthderivative(f, t, n)

    Return the `n`th derivative of `f` at `t`, computed with nested ForwardDiff.
    """
    function nthderivative(f, t, n)
        n == 0 && return f(t)
        ForwardDiff.derivative(τ -> nthderivative(f, τ, n-1), t)
    end
end


@testitem "math: complex log chooses its branch by q[1]/abs(q)" tags=[:unit, :fast] begin
    using Random
    # The one-argument `atan` returns the angle whose cosine has a nonnegative real part,
    # so the branch depends on the sign of real(q[1]/abs(q)), not of real(q[1]).  The two
    # differ for about 17% of these inputs.
    rng = Random.Xoshiro(2)
    ok(a, b) = maximum(abs, components(a) - components(b)) ≤ 1e-10 * max(1, maximum(abs, components(b)))
    nfail = count(1:2000) do _
        q = Quaternion{ComplexF64}(randn(rng, ComplexF64, 4)...)
        !ok(exp(log(q)), q)
    end
    @test nfail == 0
    nfail = count(1:2000) do _
        q = Quaternion{ComplexF64}(randn(rng, ComplexF64, 4)...)
        R = Rotor{ComplexF64}((components(q) / abs(q))...)
        !ok(exp(log(R)), R)
    end
    @test nfail == 0
    # A specific failing input from the review
    q = Quaternion{ComplexF64}(0.2565-0.4387im, 0.0510-0.0888im, 0.7445+0.5311im, -0.2963-0.5571im)
    @test ok(exp(log(q)), q)
end


@testitem "math: complex log with zero scalar part" tags=[:unit, :fast] begin
    # The general formula would divide by q[1] == 0 and return NaN.
    for T ∈ (Float64, Float32, BigFloat)
        qi = Quaternion{Complex{T}}(0, 1, 0, 0)
        @test log(qi) ≈ Quaternion{Complex{T}}(0, T(π)/2, 0, 0)
        ri = Rotor{Complex{T}}(0, 1, 0, 0)
        @test log(ri) ≈ QuatVec{Complex{T}}(0, T(π)/2, 0, 0)
        @test angle(ri) ≈ T(π)
        @test !isnan(angle(Lorentz(rotor(imx))))
        @test exp(0.5 * log(Lorentz(rotor(imx)))) ≈ Lorentz(rotor(1, 1, 0, 0))
        @test qi^T(0.5) ≈ sqrt(qi)
    end
    # A pure-vector complex quaternion, raised to a non-integer power
    v = quatvec(1+0.2im, 0.3-0.5im, 0.7im)
    @test v^2.0 ≈ v^2
    q = Quaternion{ComplexF64}(0, 0.3+0.4im, 1-0.2im, 0.5im)
    @test exp(log(q)) ≈ q
end


@testitem "math: complex log is differentiable at zero scalar part" tags=[:unit, :fast] begin
    using ForwardDiff
    # A branch that returned the angle π/2 whenever q[1] == 0 gave zero derivatives with
    # respect to q[1] there.  Compare with central differences in BigFloat.
    base = [0.0+0im, 0.6+0.2im, 0.3+0im, -0.5im]
    logq(c) = collect(components(log(Quaternion{eltype(c)}(c...))))
    function logR(c)
        n = sqrt(sum(z -> z*z, c))  # Normalize by the complex spinor norm
        collect(components(log(Rotor{eltype(c)}((c ./ n)...))))
    end
    e₁ = [1, 0, 0, 0]
    for f ∈ (logq, logR)
        d = ForwardDiff.derivative(t -> f(base .+ t .* e₁), 0.0)
        h = big"1e-30"
        dref = setprecision(BigFloat, 256) do
            b = Complex{BigFloat}.(base)
            (f(b .+ h .* e₁) - f(b .- h .* e₁)) / (2h)
        end
        @test maximum(abs, d - dref) ≤ 1e-14
        @test maximum(abs, d) > 0.5
    end
    q = Quaternion{ComplexF64}(base...)
    @test exp(log(q)) ≈ q
end


@testitem "math: norm and isapprox for complex components" tags=[:unit, :fast] begin
    import LinearAlgebra: norm
    qa = Quaternion{ComplexF64}(1, 2im, 3, 4)
    @test norm(qa) isa Float64
    @test norm(qa) ≈ sqrt(30)
    @test norm(Quaternion{ComplexF64}(1, 1im, 0, 0)) ≈ sqrt(2)
    @test abs2(Quaternion{ComplexF64}(1, 1im, 0, 0)) == 0  # The spinor norm is unchanged
    @test qa ≈ qa + 1e-12
    @test qa ≉ qa + 1e-3
    L = Boost(0.3, [0, 0, 1.0]) * Lorentz(rotor(1, 2, 3, 4))
    @test norm(L) ≈ norm(components(L))
    @test L ≈ L + 1e-12
    @test norm(quaternion(1.0, 2, 4, 10)) == 11.0
end


@testitem "math: sqrt of a QuatVec is a Quaternion" tags=[:unit, :fast] begin
    for T ∈ (Float64, Float32, Float16, BigFloat)
        v = quatvec(T(1), T(2), T(3))
        s = sqrt(v)
        @test s isa Quaternion{T}
        @test s * s ≈ v rtol=10eps(T)
        @test s ≈ sqrt(quaternion(v))
        @test !iszero(s[1])
    end
    vc = QuatVec{ComplexF64}(1, 2im, 3)
    @test sqrt(vc) isa Quaternion{ComplexF64}
    @test sqrt(vc) * sqrt(vc) ≈ vc
end


@testitem "math: sqrt does not overflow or underflow" tags=[:unit, :fast] begin
    # Compare with Base's complex `sqrt`, which handles all of these inputs.
    for (a, b) ∈ (
        (-1.0, 1e-200), (-1.0, 1e-160), (1e-310, 0.0), (1e308, 1e308), (6e307, 6e307),
        (-1e-200, 1e-200), (1e-310, 1e-310), (-1e300, 1e-300), (floatmax(), floatmax()),
        (1e-320, -1e-321), (3.0, 4.0), (-3.0, 4.0), (-1e308, 1e-300), (1e308, 1e-20),
        (-1e308, 1e-20), (1e-310, 1e-320), (-1e-310, 1e-320), (-4e35, 1e-310),
    )
        r = sqrt(quaternion(a, b, 0, 0))
        c = sqrt(complex(a, b))
        @test r[1] ≈ real(c) rtol=4eps()
        @test r[2] ≈ imag(c) rtol=4eps()
        @test iszero(r[3]) && iszero(r[4])
    end
    for (a, b) ∈ ((1f38, 1f38), (-1f0, 1f-30), (1f-40, 0f0), (-floatmax(Float32), 1f0))
        r = sqrt(quaternion(a, b, 0, 0))
        c = sqrt(complex(a, b))
        @test r[1] ≈ real(c) rtol=4eps(Float32)
        @test r[2] ≈ imag(c) rtol=4eps(Float32)
    end
    @test sqrt(Rotor{Float64}(-1.0, 1e-200, 0, 0)) ≈ rotor(5e-201, 1.0, 0, 0)
    # The direction of a vector part that is tiny compared to a negative scalar part
    for q ∈ (quaternion(-1e308, 1e-10, 2e-10, 0), quaternion(-1e308, 1e-300, -3e-301, 2e-300),
        quaternion(-1.0, 1e-310, 3e-311, 0), quaternion(-3f35, 1f-42, 2f-43, 0))
        T = eltype(components(q))
        r = sqrt(q)
        v = big.(collect(vec(q)))
        @test vec(r) / absvec(r) ≈ v / sqrt(sum(abs2, v)) rtol=4eps(T)
        @test absvec(r) ≈ sqrt(-q[1]) rtol=4eps(T)
    end
end


@testitem "math: log is accurate for tiny, huge, and nearly unit inputs" tags=[:unit, :fast] begin
    using Random
    # Compare with Base's complex `log` along one quaternionic axis.
    rng = Random.Xoshiro(7)
    for _ ∈ 1:2000
        z = complex(randn(rng), randn(rng)) * exp10(rand(rng, -300:300))
        rand(rng) < 0.3 && (z = complex(real(z), imag(z) * exp10(-rand(rng, 1:20))))
        l = log(quaternion(real(z), imag(z), 0, 0))
        lz = log(z)
        @test l[1] ≈ real(lz) atol=4eps()*max(1, abs(real(lz)))
        @test l[2] ≈ imag(lz) atol=4eps()*max(1, abs(imag(lz)))
    end
    # Near |q| = 1, the scalar part must keep its relative accuracy.
    @test log(quaternion(1.0, 1e-8, 0, 0))[1] ≈ real(log(1.0 + 1e-8im)) rtol=4eps()
    @test log(quaternion(1.0 - 1e-9, 3e-5, 0, 0))[1] ≈ real(log(1.0 - 1e-9 + 3e-5im)) rtol=4eps()
    # Subnormal and nearly negative-real inputs gave Inf or NaN.
    @test log(quaternion(-1.0, 1e-320, 0, 0)) ≈ quaternion(0, π, 0, 0)
    @test log(Rotor{Float64}(-1.0, 1e-320, 0, 0)) ≈ quatvec(π, 0, 0)
    @test log(quaternion(-1e10, 1e-310, 0, 0)) ≈ quaternion(log(1e10), π, 0, 0)
    @test log(quaternion(1e-170, 0, 0, 0)) ≈ quaternion(log(1e-170), 0, 0, 0) rtol=2eps()
    @test log(quaternion(1f-23, 0, 0, 0)) ≈ quaternion(log(1f-23), 0, 0, 0)
    @test all(isfinite, components(log(quaternion(1e-310, 1e-310, 0, 0))))
end

@testitem "math: derivatives of log at huge and tiny norms" tags=[:unit, :fast] begin
    import ForwardDiff
    using LinearAlgebra: norm
    # Since log(sq) = log(s) + log(q) for a real s > 0, the Jacobian of `log` at sq is that
    # at q divided by s.  The derivatives of `x / y` and `atan(y, x)` overflowed or underflowed
    # beyond about 1e154 or below about 1e-146, giving wrong values or NaN.
    f(x) = collect(components(log(Quaternion(x...))))
    for q ∈ (Quaternion(0.3, -0.5, 0.7, 0.2), Quaternion(-0.3, -0.5, 0.7, 0.2),
             Quaternion(0.9, 1e-9, 0, 0), Quaternion(-0.9, 0, 0, 0))
        x = collect(components(q))
        J = ForwardDiff.jacobian(f, x)
        for s ∈ (1e100, 1e160, 1e300, 1e-140, 1e-160, 1e-300)
            @test f(s * x) ≈ f(x) + [log(s), 0, 0, 0] rtol=4eps()
            Js = ForwardDiff.jacobian(f, s * x)
            @test all(isfinite, Js)
            @test norm(s * Js - J) ≤ 8eps() * norm(J)
        end
    end
    # The rescaling also applies to other float types.
    @test log(Quaternion{Float32}(1f30, 2f30, 0, 0)) ≈ Quaternion{Float32}(log(5f0) / 2 + log(1f30), atan(2f0), 0, 0)
    @test log(Quaternion{BigFloat}(1e200, 1, 2, 3))[1] ≈ log(big(1e200)) rtol=4eps(BigFloat)
end


@testitem "math: integer and rational components" tags=[:unit, :fast] begin
    @test log(quaternion(0, 0, 0, 0)) == quaternion(-Inf, 0, 0, 0)
    @test log(quaternion(-2, 0, 0, 0)) ≈ quaternion(log(2), 0, 0, π)
    @test log(quaternion(-2//1, 0, 0, 0)) isa QuaternionF64
    @test sqrt(quaternion(2, 0, 0, 0)) ≈ quaternion(√2, 0, 0, 0)
    @test sqrt(quaternion(4, 0, 0, 0)) isa QuaternionF64
    @test sqrt(quaternion(1, 2, 3, 4))^2 ≈ quaternion(1, 2, 3, 4)
    @test sqrt(Rotor{Int}(0, 1, 0, 0)) ≈ rotor(1, 1, 0, 0)
    @test angle(quaternion(-1, 0, 0, 0)) ≈ 2π
    @test quaternion(-4, 0, 0, 0)^0.5 ≈ quaternion(0, 0, 0, 2) atol=4eps()
    @test (@inferred log(quaternion(1//2, 1//3, 0, 0))) isa QuaternionF64
    @test (@inferred log(Rotor{Rational{Int}}(1//2, 1//3, 0, 0))) isa QuatVecF64
    @test (@inferred exp(quatvec(0//1, 0//1, 0//1))) isa RotorF64
    @test (@inferred exp(quatvec(1//2, 0//1, 0//1))) isa RotorF64
end


@testitem "math: inv does not overflow or underflow" tags=[:unit, :fast] begin
    @test inv(quaternion(1e200)) ≈ quaternion(1e-200)
    @test inv(quaternion(1e-200)) ≈ quaternion(1e200)
    @test inv(quaternion(1e200, 1e200, 0, 0)) ≈ quaternion(0.5e-200, -0.5e-200, 0, 0)
    @test inv(quatvec(1e-200, 0, 0)) ≈ quatvec(-1e200, 0, 0)
    @test inv(quatvec(1e-200, 0, 0)) isa QuatVecF64
    @test inv(quaternion(1f30, 0, 0, 0)) ≈ quaternion(1f-30)
    @test inv(quaternion(Float16(1e-3), 0, 0, 0)) ≈ quaternion(Float16(1e3))
    # A largest component below floatmin, whose rescaling factor 2⁻ᵉ is not representable
    @test inv(quaternion(1e-308)) ≈ quaternion(1/big(1e-308))
    @test inv(quaternion(1e-308, 1e-308, 0, 0)) ≈ quaternion(0.5/big(1e-308), -0.5/big(1e-308), 0, 0)
    @test inv(quatvec(1e-308, 0, 0)) ≈ quatvec(-1/big(1e-308), 0, 0)
    @test inv(quaternion(Float16(3e-5))) ≈ quaternion(1/Float16(3e-5))  # A subnormal input
    @test inv(quaternion(1f-41))[1] == Inf  # The true inverse overflows.
    @test !any(isnan, components(inv(quaternion(1f-39, 1f-39, 0, 0))))
    @test inv(quaternion(1e-320)) == quaternion(Inf, -0.0, -0.0, -0.0)
    q = quaternion(1.0, 2.0, 3.0, 4.0)
    @test inv(q) == conj(q) / abs2(q)
    @test all(isnan, components(inv(quaternion(0.0))))  # 0/0, as before
end


@testitem "math: exp of a Rotor" tags=[:unit, :fast] begin
    R = rotor(1.0, 2.0, 3.0, 4.0)
    @test exp(R) isa QuaternionF64
    @test exp(R) ≈ exp(quaternion(R))
    @test ℯ^R ≈ exp(quaternion(R))
end


@testitem "math: angle is accurate, with accurate gradients, near -1" tags=[:unit, :fast] setup=[MathRegressions] begin
    using ForwardDiff
    # angle = 2atan(|v⃗|, w), whose gradient is computed analytically here.
    function ∇angle(x)
        w, v = x[1], x[2:4]
        a = sqrt(sum(abs2, v))
        n = w^2 + a^2
        [-2a/n; 2w .* v ./ (a*n)]
    end
    for x ∈ ([-1.0, 1e-17, 2e-17, 0], [-1.0, 1e-8, 2e-8, 0], [0.3, 0.5, -0.2, 0.1])
        @test ForwardDiff.gradient(x -> angle(Rotor(x...)), x) ≈ ∇angle(x) rtol=10eps()
        @test ForwardDiff.gradient(x -> angle(quaternion(x...)), x) ≈ ∇angle(x) rtol=10eps()
    end
    # The kinks at ±1 give zero derivatives rather than NaN
    @test ForwardDiff.gradient(x -> angle(Rotor(x...)), [-1.0, 0, 0, 0]) == zeros(4)
    @test ForwardDiff.gradient(x -> angle(Rotor(x...)), [1.0, 0, 0, 0]) == zeros(4)
    @test angle(rotor(-1.0)) == 2π
    @test angle(quaternion(0.0)) == 0
    @test angle(exp(6.2 * imz / 2)) ≈ 6.2
end


@testitem "math: Rotor^s is accurate across the series threshold" tags=[:unit, :fast] setup=[MathRegressions] begin
    using .MathRegressions: rotorpower, maxerror
    n̂ = [1.0, -2.0, 3.0] / sqrt(14)
    for T ∈ (Float64, Float32, Float16), s ∈ (0.3, 2.5, 10.0, -0.7)
        worst = maximum([0; exp10.(range(-9, 0.4, length=200))]) do θ
            R = Rotor{T}(cos(T(θ)/2), (sin(T(θ)/2) .* T.(n̂))...)
            maxerror(R^T(s), rotorpower(R, T(s)))
        end
        @test worst ≤ 2eps(T) * max(1, abs(s))
    end
    # For high precision, the series in y = |f 𝐯|² must stop at a smaller threshold.  Here
    # y = (sθ/2)² lies on either side of it.
    setprecision(BigFloat, 256) do
        y₀ = Quaternionic.trigseriestolerance(BigFloat)
        n̂ = [big(1), -2, 3] / sqrt(big(14))
        worst = maximum(Iterators.product(big.((1e-6, 1e-5, 3e-5)), big.((0.5, 0.99, 1.01, 3.0)))) do (θ, f)
            s = 2sqrt(f*y₀)/θ
            R = Rotor{BigFloat}(cos(θ/2), (sin(θ/2) .* n̂)...)
            maxerror(R^s, rotorpower(R, s)) / max(1, Float64(s))
        end
        @test worst ≤ 2eps(BigFloat)
        R = Rotor{BigFloat}(1, big"1e-9", 0, 0)
        @test maxerror(R^big"12247.5", rotorpower(R, big"12247.5")) ≤ 2eps(BigFloat) * 12247.5
    end
end


@testitem "math: exp is accurate across the series threshold" tags=[:unit, :fast] begin
    # The series for cos and sinc stop after the y⁵ term, so their threshold must be smaller
    # for high precision.
    for p ∈ (53, 113, 256, 1000)
        setprecision(BigFloat, p) do
            y₀ = Quaternionic.trigseriestolerance(BigFloat)
            worst = maximum(big.((0.5, 0.99, 1.01, 2.0, 100.0))) do f
                t = sqrt(f * y₀ / 6)
                a = sqrt(big(6)) * t  # The norm of the vector part
                ref = [cos(a); sin(a) / a .* [t, 2t, -t]]
                e₁ = maximum(abs, collect(components(exp(quatvec(t, 2t, -t)))) - ref)
                e₂ = maximum(abs, collect(components(exp(quaternion(0, t, 2t, -t)))) - ref)
                max(e₁, e₂)
            end
            @test worst ≤ 2eps(BigFloat)
        end
    end
    for T ∈ (Float64, Float32, Float16)
        @test Quaternionic.trigseriestolerance(T) == Quaternionic.seriestolerance(T)
    end
end


@testitem "math: Rotor^s near -1 and at the identity" tags=[:unit, :fast] setup=[MathRegressions] begin
    using ForwardDiff
    using .MathRegressions: rotorpower, maxerror
    # A tiny vector part must keep its direction.
    R = Rotor{Float64}(-1.0, 1e-17, 0, 0)
    @test maxerror(R^0.5, rotorpower(R, 0.5)) ≤ 2eps()
    @test (R^0.5)[2] ≈ 1
    R32 = Rotor{Float32}(-1f0, 1f-8, 0, 0)
    @test maxerror(R32^0.5f0, rotorpower(R32, 0.5f0)) ≤ 2eps(Float32)
    @test all(isfinite, components(Rotor{Float64}(-1.0, 1e-310, 0, 0)^0.5))
    # Exactly at -1, the result follows the π𝐤 convention of `log`.
    @test (-one(RotorF64))^0.5 ≈ rotor(0, 0, 0, 1) atol=eps()
    # A dual-number exponent must work at -1, and slerp of antipodal rotors.
    @test ForwardDiff.derivative(s -> components((-one(RotorF64))^s), 0.3) ≈
        π .* [-sinpi(0.3), 0, 0, cospi(0.3)]
    @test ForwardDiff.derivative(τ -> components(slerp(one(RotorF64), -one(RotorF64), τ)), 0.3) ≈
        π .* [-sinpi(0.3), 0, 0, cospi(0.3)]
    # The result type is the promoted type in every branch.
    for R ∈ (Rotor{Float32}(1, 0, 0, 0), Rotor{Float32}(-1, 0, 0, 0), rotor(1f0, 2, 3, 4))
        @test (@inferred R^0.3) isa RotorF64
    end
    @test (@inferred rotor(Float16(1), 2, 3, 4)^0.5f0) isa RotorF32
    @test (@inferred Rotor{Float64}(-1, 0, 0, 0)^big(0.3)) isa Rotor{BigFloat}
    @test (@inferred slerp(rotor(1f0, 2, 3, 4), rotor(2f0, 1, 3, 4), 0.3)) isa RotorF64
end


@testitem "math: Rotor^s with complex exponents and Lorentz rotors" tags=[:unit, :fast] begin
    # A real rotor raised to a complex power near the identity
    R0 = rotor(1, 1e-3, 0, 0)
    s = 0.5 + 0.1im
    @test R0^s ≈ exp(s * log(R0))
    @test one(RotorF64)^(0.5im) ≈ one(Rotor{ComplexF64})
    @test exp(0.05imx/2)^s ≈ exp(s * log(exp(0.05imx/2)))
    # Lorentz rotors
    L = Boost(0.3, [0, 0, 1.0]) * Lorentz(rotor(1, 2, 3, 4))
    M = Lorentz(rotor(1, 0, 0, 0.1))
    @test L^0.3 ≈ exp(0.3 * log(quaternion(L)))
    @test (L^0.5)^2 ≈ L
    @test L^0.5 isa Lorentz{Float64}
    @test Boost(0.3, [0, 0, 1.0])^0.5 ≈ Boost(0.15, [0, 0, 1.0])
    @test one(L)^0.5 ≈ one(L)
    @test Lorentz(rotor(imx))^0.5 ≈ Lorentz(rotor(1, 1, 0, 0))
    @test quaternion(slerp(L, M, 0.3)) ≈ (quaternion(M) / quaternion(L))^0.3 * quaternion(L)
    @test squad([L, M, L*M, M*L], [0.0, 1, 2, 3], [0.5, 1.5]) isa Vector
    Ls = Boost(1e-6, [0, 0, 1.0]) * Lorentz(exp(1e-6imx))
    @test Ls^2.5 ≈ Ls * Ls * sqrt(Ls)
    # A quaternion exponent gives the quaternion exp(s log(q))
    R = rotor(1.0, 2, 3, 4)
    p = quaternion(0.5, 0.1, 0, 0)
    @test R^p ≈ quaternion(R)^p
end


@testitem "math: rational exponents" tags=[:unit, :fast] setup=[MathRegressions] begin
    using .MathRegressions: rotorpower, maxerror
    # These were method ambiguities with Base's `^(::Number, ::Rational)`.
    R = rotor(big(1.0), 2, 3, 4)
    @test maxerror(R^(1//3), rotorpower(R, big(1)/3)) ≤ 10eps(BigFloat)
    q = quaternion(big(1.0), 2, 3, 4)
    @test maxerror((q^(1//2))^2, components(q)) ≤ 100eps(BigFloat)
    v = quatvec(big(1.0), 2, 3)
    @test maxerror((v^(1//2))^2, components(v)) ≤ 100eps(BigFloat)
    @test rotor(1.0, 2, 3, 4)^(1//2) ≈ rotor(1.0, 2, 3, 4)^0.5
    @test (@fastmath rotor(1.0, 2, 3, 4)^(1//2)) ≈ rotor(1.0, 2, 3, 4)^0.5
    @test quaternion(1, 2, 3, 4)^(1//2) ≈ sqrt(quaternion(1, 2, 3, 4))
    q1, q2 = rotor(1.0, 2, 3, 4), rotor(2.0, 1, 3, 4)
    @test slerp(q1, q2, 1//3) ≈ slerp(q1, q2, 1/3)
end


@testitem "math: integer powers match Base.power_by_squaring" tags=[:unit, :fast] begin
    # The powers are computed by a hand-written loop (Enzyme's reverse mode aborts on
    # `Base.power_by_squaring`), which must give identical results and types.
    qs = (quatvec(1.0, 2, 3), quaternion(1.0, 2, 3, 4), rotor(1.0, 2, 3, 4), quaternion(1, 2, 3, 4),
        quatvec(1, 2, 3), quaternion(1f0, 2, 3, 4), rotor(1f0, 2, 3, 4), Rotor{Int}(1, 0, 0, 0),
        quaternion(0.3, -1.1, 0.7, 0.2), Lorentz(rotor(1, 2, 3, 4)), imx,
        quaternion(true, false, true, false))
    mismatches = count(Iterators.product(qs, 0:20)) do (q, n)
        a = q^n
        b = Base.power_by_squaring(q, n)
        !(typeof(a) === typeof(b) && components(a) == components(b))
    end
    @test mismatches == 0
    @test imx^0 isa Quaternion{Int} && imx^1 isa Quaternion{Int} && imx^3 isa Quaternion{Int}
    mismatches = count(Iterators.product(qs, 1:20)) do (q, n)
        eltype(components(q)) <: Integer && return false
        a = q^(-n)
        b = inv(Base.power_by_squaring(q, n))
        !(typeof(a) === typeof(b) && components(a) == components(b))
    end
    @test mismatches == 0
end


@testitem "math: reverse mode at rotors with zero scalar part" tags=[:unit, :fast] setup=[MathRegressions] begin
    using ForwardDiff, ReverseDiff, Random
    # `x = abs2vec(q)/q[1]^2` used to be computed before testing q[1] > 0, so reverse
    # mode recorded a division by zero and returned NaN gradients.
    W = randn(Random.Xoshiro(1234), 4)
    wsum(q) = sum(W .* components(q))
    for x ∈ ([0.0, 0.3, -0.4, 1.2], [0.0, 1.0, 0.0, 0.0], [1e-300, 0.3, -0.4, 1.2])
        for s ∈ (1.3, 0.5, -0.7, 50.0)
            f(x) = wsum(Rotor(x...)^s)
            g = ReverseDiff.gradient(f, x)
            @test all(isfinite, g)
            @test g ≈ ForwardDiff.gradient(f, x)
        end
    end
    fs(x) = wsum(slerp(one(RotorF64), Rotor(x...), 0.37))
    x = [0.0, 0.6, 0.8, 0.0]
    @test ReverseDiff.gradient(fs, x) ≈ ForwardDiff.gradient(fs, x)
    @test all(isfinite, ReverseDiff.gradient(fs, x))
end


@testitem "math: nested derivatives of exp and log near the identity" tags=[:unit] setup=[MathRegressions] begin
    using .MathRegressions: nthderivative
    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    v = vec(n̂)
    E(t) = collect(components(exp(t * n̂ / 2)))
    Eref(t) = [cos(t/2); sin(t/2) .* v]
    Q(t) = collect(components(exp(quaternion(0.5, (t .* v)...))))
    Qref(t) = exp(0.5) .* [cos(t); sin(t) .* v]
    Λ(t) = collect(components(log(Rotor(cos(t), (sin(t) .* v)...))))
    Λref(t) = [0; t .* v]
    for t ∈ (1e-8, 1e-5, 1e-3), n ∈ 1:4
        @test maximum(abs, nthderivative(E, t, n) - nthderivative(Eref, t, n)) ≤ 1e-14
        @test maximum(abs, nthderivative(Q, t, n) - nthderivative(Qref, t, n)) ≤ 1e-14
        @test maximum(abs, nthderivative(Λ, t, n) - nthderivative(Λref, t, n)) ≤ 1e-14
    end
end


@testitem "math: Taylor-series coefficients" tags=[:unit, :fast] begin
    import Quaternionic: cosseries, sincseries, atanseries, logseries
    # Each series is compared with its exact function at points where the truncation error
    # is far below an error in any one coefficient.
    setprecision(BigFloat, 256) do
        for y ∈ big.((1e-1, 1e-2))
            @test abs(cosseries(y) - cos(sqrt(y))) ≤ y^6 / factorial(big(12))
            @test abs(sincseries(y) - sin(sqrt(y)) / sqrt(y)) ≤ y^6 / factorial(big(13))
        end
        # At the threshold, the truncation error is below ε.
        y = Quaternionic.trigseriestolerance(BigFloat)
        @test abs(cosseries(y) - cos(sqrt(y))) ≤ eps(BigFloat)
        @test abs(sincseries(y) - sin(sqrt(y)) / sqrt(y)) ≤ eps(BigFloat)
        x = Quaternionic.seriestolerance(BigFloat)
        @test abs(atanseries(x) - atan(sqrt(x)) / sqrt(x)) ≤ eps(BigFloat)
        for x ∈ big.((1e-2, 1e-3))
            @test abs(atanseries(x) - atan(sqrt(x)) / sqrt(x)) ≤ x^8 / 17
            @test abs(logseries(x) - atan(sqrt(x)) * sqrt(1 + x) / sqrt(x)) ≤ x^8 / 16
        end
    end
end
