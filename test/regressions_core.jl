# Regression tests for the core types: construction, promotion, conversion, equality,
# hashing, display, and the `_sincu` helper.

@testitem "core: QuatVec equality with quaternions and numbers" tags=[:unit, :fast] begin
    v = quatvec(1.0, 2.0, 3.0)
    q0 = quaternion(0.0, 1.0, 2.0, 3.0)
    q5 = quaternion(5.0, 1.0, 2.0, 3.0)

    # The scalar part of the other operand must be zero
    @test v == q0
    @test q0 == v
    @test v != q5
    @test q5 != v
    @test !isequal(v, q5)
    @test !isequal(q5, v)
    @test v != rotor(1.0, 1.0, 2.0, 3.0)
    @test isequal(v, q0)
    @test isequal(q0, v)

    # `==` is transitive through a QuatVec
    @test !(q5 == v && v == q0)

    # A QuatVec equals a number exactly when both are zero
    @test zero(QuatVecF64) == 0
    @test 0 == zero(QuatVecF64)
    @test 𝐢 - 𝐢 == 0
    @test 0.0 == 𝐢 - 𝐢
    @test count(==(0), [𝐢 - 𝐢, quatvec(0.0, 0.0, 0.0), 𝐢]) == 2
    @test quatvec(1.0, 0.0, 0.0) != 0
    @test zero(QuatVecF64) != 1
    @test 1 != zero(QuatVecF64)
end

@testitem "core: isequal agrees with hash across types" tags=[:unit, :fast] begin
    pairs = [
        (quatvec(1.0, 2.0, 3.0), quaternion(0.0, 1.0, 2.0, 3.0)),
        (quatvec(1.0, 2.0, 3.0), quaternion(-0.0, 1.0, 2.0, 3.0)),
        (quatvec(1.0, 2.0, 3.0), quaternion(5.0, 1.0, 2.0, 3.0)),
        (quatvec(1, 2, 3), quaternion(0.0, 1.0, 2.0, 3.0)),
        (zero(QuatVecF64), 0.0),
        (zero(QuatVecF64), -0.0),
        (zero(QuatVecF64), 0),
        (quaternion(1.0), 1.0),
        (quaternion(1.0), 1),
        (quaternion(1.0, -0.0, 0.0, 0.0), 1.0),
        (quaternion(0.0), -0.0),
        (quaternion(NaN), NaN),
        (quaternion(2.5), 5//2),
        (quaternion(1.0 + 2.0im), 1.0 + 2.0im),
    ]
    for (a, b) in pairs
        @test isequal(a, b) == isequal(b, a)
        if isequal(a, b)
            @test hash(a) == hash(b)
        end
    end
    @test isequal(quatvec(1.0, 2.0, 3.0), quaternion(0.0, 1.0, 2.0, 3.0))
    @test !isequal(quatvec(1.0, 2.0, 3.0), quaternion(-0.0, 1.0, 2.0, 3.0))
    @test isequal(zero(QuatVecF64), 0.0)
    @test !isequal(zero(QuatVecF64), -0.0)
    @test isequal(quaternion(1.0), 1.0)
    @test !isequal(quaternion(1.0, -0.0, 0.0, 0.0), 1.0)
    @test !isequal(quaternion(0.0), -0.0)
    @test isequal(quaternion(NaN), NaN)

    # Set and Dict lookups merge values that are `isequal`
    v = quatvec(1.0, 2.0, 3.0)
    q0 = quaternion(0.0, 1.0, 2.0, 3.0)
    @test length(Set([v, q0, quatvec(1, 2, 3)])) == 1
    @test haskey(Dict(v => 1), q0)
    @test !haskey(Dict(v => 1), quaternion(5.0, 1.0, 2.0, 3.0))
end

@testitem "core: QuatVec never stores a scalar part" tags=[:unit, :fast] begin
    for T in (Float64, Float32, Int)
        q = quaternion(T(1), T(2), T(3), T(4))
        @test components(convert(QuatVec{T}, q))[1] == 0
        @test components(QuatVec{T}(1, 2, 3, 4))[1] == 0
        @test iszero(QuatVec{T}(1))
        w = QuatVec{T}[]
        push!(w, q)
        push!(w, T(2))
        @test all(x -> iszero(components(x)[1]), w)
        @test iszero(w[2])
        @test w[1] == quatvec(q)
    end

    # The multiplicative identity is not a pure vector, so these return Quaternions
    @test one(QuatVecF64) === one(QuaternionF64)
    @test oneunit(QuatVecF64) === one(QuaternionF64)
    @test oneunit(quatvec(1.0, 2.0, 3.0)) === one(QuaternionF64)
    # `ones` keeps the requested element type, so the scalar part of the identity is
    # discarded, and no hidden scalar part is stored
    @test ones(QuatVecF64, 2) isa Vector{QuatVecF64}
    @test ones(QuatVecF64, 2, 3) isa Matrix{QuatVecF64}
    @test ones(QuatVecF64, (2, 3)) isa Matrix{QuatVecF64}
    @test ones(QuatVecF64) isa Array{QuatVecF64,0}
    @test all(x -> iszero(components(x)), ones(QuatVecF64, 3))
    @test zeros(QuatVecF64, 2) isa Vector{QuatVecF64}
    @test all(iszero, zeros(QuatVecF64, 2))

    # The UnionAll types follow the same conventions
    @test isone(one(QuatVec))
    @test isone(prod(QuatVec[]))
    @test iszero(zero(Rotor))
    @test iszero(sum(Rotor[]))

    # A product of two equal-looking vectors does not depend on a hidden scalar
    x = convert(QuatVec{Float64}, quaternion(1.0, 2.0, 3.0, 4.0))
    y = quatvec(1.0, 2.0, 3.0, 4.0)
    @test x * x == y * y
    @test hash(x) == hash(y)
end

@testitem "core: Rotor promotion and conversion" tags=[:unit, :fast] begin
    R = rotor(1, 2, 3, 4)

    # Mixing a Rotor with a plain number gives a Quaternion
    @test promote_type(RotorF64, Float64) === QuaternionF64
    @test promote_type(RotorF32, Float64) === QuaternionF64
    @test promote_type(RotorF64, Int) === QuaternionF64
    @test promote_type(RotorF64, ComplexF64) === Quaternion{ComplexF64}
    @test promote_type(RotorF64, QuaternionF32) === QuaternionF64
    @test promote_type(RotorF64, RotorF32) === RotorF64
    @test typeof(promote(R, 2.0)) === Tuple{QuaternionF64, QuaternionF64}
    @test promote(R, 2.0)[2] == quaternion(2.0)
    @test [R, 2.0] isa Vector{QuaternionF64}
    @test [R, 1im] isa Vector{Quaternion{ComplexF64}}

    # Conversion to Rotor{T} normalizes, like `rotor` and `convert(Rotor, x)`
    @test convert(Rotor, 3.0) == rotor(1.0)
    @test convert(RotorF64, 3.0) == rotor(1.0)
    @test convert(RotorF64, -3.0) == rotor(-1.0, 0, 0, 0)
    @test convert(RotorF64, quaternion(1, 2, 3, 4)) ≈ R
    @test convert(RotorF64, R) === R
    @test convert(RotorF32, R) === RotorF32(R)
    @test convert(Rotor{BigFloat}, quaternion(1, 2, 3, 4)) ≈ rotor(big.((1, 2, 3, 4))...) atol=10eps(BigFloat)
    w = RotorF64[]
    push!(w, 2.0)
    push!(w, quaternion(1, 2, 3, 4))
    push!(w, quatvec(0, 0, 3))
    @test w[1] == rotor(1.0)
    @test w[2] ≈ R
    @test w[3] == rotor(0, 0, 0, 1)
    @test all(x -> abs2(quaternion(x)) ≈ 1, w)
    w[1] = quaternion(0, 5, 0, 0)
    @test w[1] == rotor(0, 1, 0, 0)

    # The zero quaternion is not a unit quaternion, so `zero` gives a Quaternion.  But `zeros`
    # keeps the element type, filling the array with unnormalized zero rotors (rather than
    # the NaN that normalizing them would give), so that it can be used for preallocation.
    @test zero(RotorF64) === zero(QuaternionF64)
    @test zeros(RotorF64, 2) isa Vector{RotorF64}
    @test zeros(RotorF64, 2, 3) isa Matrix{RotorF64}
    @test zeros(RotorF64, (2,)) isa Vector{RotorF64}
    @test zeros(RotorF64, Int32(2)) isa Vector{RotorF64}
    @test zeros(RotorF64) isa Array{RotorF64,0}
    @test zeros(Rotor{Float32}, 2) isa Vector{Rotor{Float32}}
    @test all(x -> iszero(components(x)), zeros(RotorF64, 3))
    let R⃗ = zeros(RotorF64, 4)
        for i ∈ eachindex(R⃗)
            R⃗[i] = rotor(1.0, i, 0.5i, -0.25i)
        end
        @test R⃗ isa Vector{RotorF64}
        @test all(x -> abs2(quaternion(x)) ≈ 1, R⃗)
        @test squad(R⃗, 1.0:4.0, 2.5) isa RotorF64
    end

    # Rounding a Rotor gives a Quaternion
    @test round(rotor(1.0, 1.1, 0.2, 0.3)) === quaternion(1.0, 1.0, 0.0, 0.0)
    @test round(quaternion(1.2, 2.7, 0.2, 0.3)) === quaternion(1.0, 3.0, 0.0, 0.0)
end

@testitem "core: single-argument rotor" tags=[:unit, :fast] begin
    @test rotor(2) === rotor(2, 0, 0, 0)
    @test rotor(2) isa RotorF64
    @test rotor(true) isa RotorF64
    @test rotor(-2.0) == rotor(-1.0, 0, 0, 0)
    @test all(isnan, components(rotor(0)))
    @test all(isnan, components(rotor(-0.0)))
    @test isnan(rotor(NaN))
    @test rotor(2.0 + 1im) ≈ Rotor{ComplexF64}(1, 0, 0, 0)
    @test Base.return_types(rotor, (Vector{Int},)) == [RotorF64]
    @test rotor([2]) === rotor(2)
end

@testitem "core: indexing" tags=[:unit, :fast] begin
    q = quaternion(1.0, 2.0, 3.0, 4.0)
    @test q[end] == 4.0
    @test q[begin] == 1.0
    @test firstindex(q) == 1
    @test lastindex(q) == 4
    @test q[begin:end] == [1.0, 2.0, 3.0, 4.0]
    @test q[2:end] == [2.0, 3.0, 4.0]
    @test q[:] == [1.0, 2.0, 3.0, 4.0]
    @test q[CartesianIndex()] === q
    v = quatvec(2.0, 3.0, 4.0)
    @test v[end] == 4.0
    @test v[:] == [0.0, 2.0, 3.0, 4.0]
end

@testitem "core: conversion to real types" tags=[:unit, :fast] begin
    @test Float64(quaternion(2.0)) === 2.0
    @test Float32(quaternion(2.0)) === 2.0f0
    @test Int(quaternion(3)) === 3
    @test Bool(quaternion(true)) === true
    @test convert(Float64, quaternion(2.0)) === 2.0
    @test round(Int, quaternion(2.2)) === 2
    @test Float64(zero(QuatVecF64)) === 0.0
    x = zeros(2)
    x[1] = quaternion(5.0)
    @test x[1] === 5.0
    @test_throws InexactError Float64(quaternion(1.0, 1.0, 0.0, 0.0))
    @test_throws InexactError Float64(quatvec(1.0, 0.0, 0.0))
    @test_throws InexactError Int(quaternion(2.5))
end

@testitem "core: isapprox with complex components" tags=[:unit, :fast] begin
    a1 = quaternion(1.0 + 1im, 2, 3, 4)
    a2 = a1 + 1e-12
    @test a1 ≈ a2
    @test a1 ≈ a2 atol=1e-8
    @test a1 ≉ a1 + 1e-3
    @test a1 ≉ a2 rtol=0 atol=1e-14
    @test quaternion(1.0 + 0im) ≈ 1.0
    @test 1.0 ≈ quaternion(1.0 + 0im)
    @test quaternion(1.0 + 1e-12im) ≈ 1.0
    @test quaternion(1.0 + 1im) ≉ 1.0
    @test Rotor{ComplexF64}(1, 0, 0, 0) ≈ quaternion(1.0)
    @test quaternion(1.0) ≈ Rotor{ComplexF64}(1, 0, 0, 0)
    @test quatvec(1.0 + 1im, 2, 3) ≈ quatvec(1.0 + 1im, 2, 3 + 1e-12)
    @test a1 ≈ a1 + 1e-6 atol=1e-5
end

@testitem "core: plain-text display" tags=[:unit, :fast] begin
    str(x) = sprint(show, x)
    # The sign of a complex component stays inside its parentheses
    @test str(quaternion(0, -1 + 2im, 0, 0)) == "0 + 0im + (-1 + 2im)𝐢 + (0 + 0im)𝐣 + (0 + 0im)𝐤"
    @test str(quaternion(0.0, -1.0 + 2.0im, -3.0 - 4.0im, 1.0 - 2.0im)) ==
        "0.0 + 0.0im + (-1.0 + 2.0im)𝐢 + (-3.0 - 4.0im)𝐣 + (1.0 - 2.0im)𝐤"
    @test str(quatvec(-1 - 2im, 0, 0)) == " + (-1 - 2im)𝐢 + (0 + 0im)𝐣 + (0 + 0im)𝐤"
    # A single negative real may keep its sign outside the parentheses
    @test str(quaternion(1.0, -3.0e-9, 2.0, -4.0)) == "1.0 - (3.0e-9)𝐢 + 2.0𝐣 - 4.0𝐤"
    @test str(quaternion(1//2, -1//3, 0//1, 0//1)) == "1//2 - (1//3)𝐢 + (0//1)𝐣 + (0//1)𝐤"
    # Non-finite and Bool components are followed by `*`, so that the output parses
    @test str(quaternion(1.0, NaN, -Inf, Inf)) == "1.0 + NaN*𝐢 - Inf*𝐣 + Inf*𝐤"
    @test str(quaternion(true, false, true, false)) == "true + false*𝐢 + true*𝐣 + false*𝐤"
    for q in (
        quaternion(1.0, 0.0, 2.0, NaN),
        quaternion(true, false, true, false),
        quaternion(1.5, -Inf, 2.0, -3.0e-9),
        quaternion(1.0 - 2.0im, -1.0 + 2.0im, 0.5im, -3.0),
    )
        p = eval(Meta.parse(repr(q)))
        @test isequal(components(p), components(q))
    end
    # Vector components are printed with the IOContext, like the scalar part
    @test sprint(show, quaternion(1/3, 2/3, 1.0, 2.0); context=:compact => true) ==
        "0.333333 + 0.666667𝐢 + 1.0𝐣 + 2.0𝐤"
end

@testitem "core: Lorentz compositions at large rapidity" tags=[:unit, :fast] begin
    relerr(a, b) = maximum(abs, components(a) .- components(b)) / maximum(abs, components(b))
    ẑ = [0.0, 0.0, 1.0]
    for η in (20.0, 30.0, 40.0, 60.0)
        exact = setprecision(BigFloat, 256) do
            Boost(big(η), big.(ẑ))
        end
        B = Boost(η / 2, ẑ)
        @test B * B isa Rotor{ComplexF64}
        @test relerr(B * B, exact) < 1e-14
        @test relerr(Boost(η / 4, ẑ) * Boost(3η / 4, ẑ), exact) < 1e-14
        @test relerr(B / Boost(-η / 2, ẑ), exact) < 1e-14
        @test all(isfinite, components(B * B))

        # Products of noncollinear boosts with a rotation
        x̂, ŷ = [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]
        R = rotor(1.0, 2.0, 3.0, 4.0)
        P = Lorentz(R) * Boost(η / 3, ŷ) * Boost(η / 3, x̂) * Boost(η / 3, ẑ)
        Pb = setprecision(BigFloat, 256) do
            Lorentz(rotor(big.((1, 2, 3, 4))...)) * Boost(big(η) / 3, big.(ŷ)) *
                Boost(big(η) / 3, big.(x̂)) * Boost(big(η) / 3, big.(ẑ))
        end
        @test relerr(P, Pb) < 1e-14
    end

    # Real rotors are still renormalized
    R1 = Rotor{Float64}(2.0, 0.0, 0.0, 0.0)  # Deliberately not normalized
    @test R1 * R1 == rotor(1.0)
end

@testitem "core: long chains of Lorentz products at moderate rapidity" tags=[:unit, :fast] begin
    using Random
    using LinearAlgebra: normalize
    # Each factor is a rotation seen in a frame boosted by rapidity 0.7.  Its spinor norm is
    # 1 only to rounding error, so a chain of products drifts unless products of complex
    # rotors are renormalized while that is accurate.  The reference multiplies the same
    # factors, each normalized exactly, in BigFloat.
    rng = Xoshiro(1234)
    n̂ = normalize(randn(rng, 3))
    P, maxerr, maxdefect = setprecision(BigFloat, 128) do
        P = Lorentz(rotor(1.0))
        Pb = Lorentz(rotor(big(1.0)))
        maxerr = maxdefect = 0.0
        for i ∈ 1:3000
            F = Boost(0.7, n̂) * Lorentz(rotor(randn(rng, 4))) * Boost(-0.7, n̂)
            P = P * F
            Pb = Pb * rotor(Complex{BigFloat}.(components(F)))
            maxerr = max(maxerr, Float64(maximum(abs, components(P) .- components(Pb))))
            maxdefect = max(maxdefect, abs(sum(z -> z * z, components(P)) - 1))
        end
        P, maxerr, maxdefect
    end
    @test P isa Rotor{ComplexF64}
    @test maxerr < 3e-14
    @test maxdefect < 1e-14
end

@testitem "core: _sincu derivatives near zero" tags=[:unit, :fast] begin
    using ForwardDiff
    import Quaternionic: _sincu

    # The `k`th derivative of sin(x)/x, from its Taylor series, evaluated in BigFloat
    function reference(x, k)
        setprecision(BigFloat, 1024) do
            b = big(x)
            sum(0:120) do n
                p = 2n
                p < k ? zero(BigFloat) :
                    (-1)^n / factorial(big(2n + 1)) * factorial(big(p)) / factorial(big(p - k)) * b^(p - k)
            end
        end
    end
    nth_derivative(f, n, t) =
        n == 0 ? f(t) : ForwardDiff.derivative(τ -> nth_derivative(f, n - 1, τ), t)

    # The series is used below 3 for Float64 and Float32, and the closed form above it, so
    # points on both sides of the threshold are included, where the errors are largest.  The
    # measured errors there are at most about 3 ulps for k ≤ 4 and 10 ulps for k = 5.
    for T in (Float64, Float32)
        points = vcat(
            T[3.5, 3.1, 3.02, 2.98, 2.5, 2.0, 1.5, 1.0, 1e-1, 1e-2, 1e-4, 1e-6, 1e-8, 1e-12, 1e-16, 1e-20, -1e-8, -2.9],
            [prevfloat(T(3), 2), prevfloat(T(3)), T(3), nextfloat(T(3)), nextfloat(T(3), 2)],
        )
        bad = count((x, k) for x in points for k in 0:5) do (x, k)
            d = nth_derivative(_sincu, k, x)
            r = reference(x, k)
            # The `k`th derivative near zero is of order 1/(k+1)
            abs(d - r) > (k ≤ 4 ? 8 : 32) * eps(T) * max(abs(r), 1 / (k + 1))
        end
        @test bad == 0
    end

    # BigFloat Duals also take the series, below a precision-dependent threshold, which is
    # about 0.0285 for 256 bits
    setprecision(BigFloat, 256) do
        bad = count(
            (x0, k) for x0 in (1e-1, 0.03, 0.0285, 0.028, 1e-3, 1e-5, 1e-10, 1e-20) for k in 0:4
        ) do (x0, k)
            x = BigFloat(x0)
            d = nth_derivative(_sincu, k, x)
            r = reference(x, k)
            abs(d - r) > 1e-66 * max(abs(r), 1 / (k + 1))
        end
        @test bad == 0
    end
end

@testitem "core: Hessians of exp near a zero vector part" tags=[:unit, :fast] begin
    using ForwardDiff
    W = [0.3, -0.7, 1.1, 0.5]
    f(x) = sum(W .* components(exp(quaternion(x...))))
    g(x) = sum(W[2:4] .* vec(exp(quatvec(x...))))

    # The references are independent of the package: the explicit formula
    # exp(a + v⃗) = eᵃ (cos|v⃗| + v⃗ sin|v⃗|/|v⃗|), evaluated in 512-bit BigFloat, where sin(r)/r
    # suffers no cancellation, and differentiated by central finite differences with a step
    # of 2⁻¹³⁰, whose truncation and rounding errors are both below 1e-70.
    function expref(a, v)
        r = sqrt(sum(abs2, v))
        vcat(exp(a) * cos(r), exp(a) * sin(r) / r .* v)
    end
    fref(x) = sum(W .* expref(x[1], x[2:4]))
    gref(x) = sum(W[2:4] .* expref(zero(eltype(x)), x)[2:4])
    function fd_hessian(F, p)
        setprecision(BigFloat, 512) do
            x = big.(p)
            h = big(2)^-130
            n = length(x)
            e(i) = [j == i ? h : zero(h) for j in 1:n]
            [(F(x + e(i) + e(j)) - F(x + e(i) - e(j)) - F(x - e(i) + e(j)) + F(x - e(i) - e(j))) / (4h^2)
                for i in 1:n, j in 1:n]
        end
    end

    for s in (1e-4, 1e-8, 1e-12, 1e-20)
        p = [0.7, s, -2s, 0.5s]
        H = ForwardDiff.hessian(f, p)
        Hb = fd_hessian(fref, p)
        @test maximum(abs, H - Hb) < 4e-15 * maximum(abs, Hb)
        pv = [s, -2s, 0.5s]
        Hv = ForwardDiff.hessian(g, pv)
        Hvb = fd_hessian(gref, pv)
        @test maximum(abs, Hv - Hvb) < 1e-14 * max(maximum(abs, Hvb), s)
    end
end

@testitem "core: random quaternion statistics" tags=[:unit] begin
    using Random
    mean(f, x) = sum(f, x) / length(x)
    mean(x) = mean(identity, x)
    var(x) = (m = mean(x); sum(y -> (y - m)^2, x) / (length(x) - 1))
    rng = MersenneTwister(20240917)
    N = 200_000
    q = randn(rng, QuaternionF64, N)
    for i in 1:4
        cᵢ = [x[i] for x in q]
        @test abs(mean(cᵢ)) < 5 * sqrt(1 / 4 / N)
        @test abs(var(cᵢ) - 1 / 4) < 0.01
    end
    @test abs(mean(abs2, q) - 1) < 0.01

    v = randn(rng, QuatVecF64, N)
    @test all(x -> iszero(x[1]), v)
    for i in 2:4
        @test abs(var([x[i] for x in v]) - 1 / 3) < 0.01
    end

    # The rotation angle of a uniformly random rotation has E[cos θ] = -1/2, where
    # cos θ = 2w² - 1 for a unit rotor with scalar part w
    R = randn(rng, RotorF64, N)
    @test abs(mean(r -> 2r[1]^2 - 1, R) + 1 / 2) < 0.01
end

@testitem "core: ReverseDiff through rotor and complex norms" tags=[:unit, :fast] begin
    using ReverseDiff, ForwardDiff

    # `_hypot` squares complex components by multiplication, which ReverseDiff can handle
    f(x) = real(abs(quaternion(complex(x[1], x[2]), complex(x[3], 0.1), 0.3, 0.2)))
    x₀ = [0.7, 0.2, -0.4]
    @test ReverseDiff.gradient(f, x₀) ≈ ForwardDiff.gradient(f, x₀) rtol=1e-14

    # `rotor` normalizes without a broadcast, so recorded tapes can be replayed
    W = [0.3, -0.2, 0.5, 0.7]
    g(x) = sum(W .* components(rotor(x[1], x[2], x[3], x[4])))
    h(x) = sum(W .* components(rotor(x[1], x[2], x[3], x[4]) * rotor(x[4], x[3], x[2], x[1])))
    q₀ = [0.4, -0.3, 0.8, 0.1]
    for F in (g, h)
        tape = ReverseDiff.GradientTape(F, q₀)
        @test ReverseDiff.gradient!(similar(q₀), tape, q₀) ≈ ForwardDiff.gradient(F, q₀) rtol=1e-14
        compiled = ReverseDiff.compile(tape)
        @test ReverseDiff.gradient!(similar(q₀), compiled, q₀) ≈ ForwardDiff.gradient(F, q₀) rtol=1e-14
    end
end

@testitem "core: reverse-mode gradients with respect to Rotor arguments" tags=[:unit, :fast] begin
    using ForwardDiff, Zygote
    # `convert(Rotor{T}, x)` normalizes, so the projection of a cotangent onto a `Rotor`
    # argument must not use `convert`; a cotangent is not a unit quaternion.  Each gradient
    # is compared with ForwardDiff with respect to the raw components.
    q = quaternion(0.3, -0.2, 0.5, 0.7)
    f(R) = (R * q)[1]
    R = rotor(1.0, 2.0, 3.0, 4.0)
    gz = Zygote.gradient(f, R)[1]
    gf = ForwardDiff.gradient(c -> f(Rotor{eltype(c)}(c)), Vector(components(R)))
    @test components(gz) ≈ gf rtol=1e-14

    fv(Rs) = sum(x -> (x * q)[2], Rs)
    Rs = [rotor(1.0, 2.0, 3.0, 4.0), rotor(0.5, -1.0, 0.2, 0.3)]
    gzv = Zygote.gradient(fv, Rs)[1]
    gfv = ForwardDiff.gradient(
        c -> fv([Rotor{eltype(c)}(c[1:4]), Rotor{eltype(c)}(c[5:8])]),
        vcat(Vector(components(Rs[1])), Vector(components(Rs[2])))
    )
    @test vcat(Vector(components(gzv[1])), Vector(components(gzv[2]))) ≈ gfv rtol=1e-14
end
