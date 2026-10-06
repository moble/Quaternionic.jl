# Regression tests for fixes to src/algebra.jl and src/Quaternionic.jl.

@testitem "Division: mixed types are not short-circuited" begin
    # `q / p` once returned exactly 1 whenever `p == q`, and `==` between a `QuatVec` and
    # another quaternion compares only the vector parts.  So a quaternion divided by a
    # `QuatVec` with the same vector part returned 1, whatever its scalar part.
    q = quaternion(5.0, 1.0, 2.0, 3.0)
    v = quatvec(1.0, 2.0, 3.0)
    @test q / v ≈ q * inv(v)
    @test v / q ≈ v * inv(q)
    @test q / v ≉ one(q)
    @test v / q ≉ one(q)

    R = rotor(1.0, 2.0, 3.0, 4.0)
    w = quatvec(R[2], R[3], R[4])
    @test R / w ≈ R * inv(w)
    @test w / R ≈ w * conj(R)
    @test v / v ≈ one(QuaternionF64)
    @test typeof(v / v) === QuaternionF64
end

@testitem "Division: q/q is exactly one for real floats" begin
    using Random
    rng = Random.MersenneTwister(1234)
    for T in (Float64, Float32, Float16, BigFloat)
        # The scales keep `abs2(q)` well inside the range of each type.
        scales = T === Float16 ? (-1:1) : T === Float32 ? (-8:8) : (-60:60)
        qs = [randn(rng, Quaternion{T}) * T(10)^rand(rng, scales) for _ in 1:2000]
        @test all(qs) do q
            r = q / q
            typeof(r) === typeof(q) && components(r) == components(one(q))
        end
    end
    Rs = [randn(rng, RotorF64) for _ in 1:2000]
    @test all(R -> (r = R / R; r isa RotorF64 && components(r) == components(one(R))), Rs)
    vs = [randn(rng, QuatVecF64) for _ in 1:2000]
    @test all(v -> components(v / v) == components(one(QuaternionF64)), vs)
end

@testitem "Division: agrees with a BigFloat reference" begin
    using Random
    rng = Random.MersenneTwister(7)
    @test all(1:1000) do _
        q = randn(rng, QuaternionF64)
        p = randn(rng, QuaternionF64)
        ref = big(q) * inv(big(p))
        maximum(abs, components(q / p - ref)) ≤ 4eps(Float64) * abs(ref)
    end
end

@testitem "Division: the result type does not depend on the values" begin
    a32 = quaternion(1f0, 2f0, 3f0, 4f0)
    @test typeof(a32 / quaternion(1, 2, 3, 4)) === typeof(a32 / quaternion(1, 2, 3, 5))
    @test (@inferred a32 / quaternion(1, 2, 3, 4)) isa QuaternionF32

    ar = quaternion(1//2, 1//3, 0//1, 1//1)
    @test (@inferred ar / ar) isa Quaternion{Rational{Int}}
    @test ar / ar == one(ar)
    @test typeof(ar / ar) === typeof(ar / quaternion(1//1, 0//1, 0//1, 0//1))

    @test (@inferred Rotor{Float32}(1, 0, 0, 0) / Rotor{Int}(1, 0, 0, 0)) isa RotorF32
    @test (@inferred QuatVec{Float32}(0, 1, 2, 3) / QuatVec{Int}(0, 1, 2, 3)) isa QuaternionF32
end

@testitem "Division: reverse-mode derivative at q == p" begin
    # A backend whose `==` compares only primal values used to take the shortcut branch
    # and return a zero derivative at q == p.
    using ForwardDiff, ReverseDiff
    W = [0.3, -1.1, 0.7, 0.45]
    wsum(q) = sum(W .* components(q))
    P0 = quaternion(0.3, 1.2, -0.5, 0.7)
    R0 = rotor(0.3, -0.7, 0.2, 0.5)
    f1(x) = wsum(quaternion(x...) / P0)
    f2(x) = wsum(P0 / quaternion(x...))
    f3(x) = wsum(slerp(R0, rotor(x...), 0.37))
    for (f, x) in ((f1, collect(components(P0))), (f2, collect(components(P0))),
                   (f3, collect(components(R0))))
        ref = ForwardDiff.gradient(f, x)
        @test !iszero(ref)
        @test ReverseDiff.gradient(f, x) ≈ ref atol=10eps()
    end
end

@testitem "Division: FastDifferentiation variables" begin
    # The shortcut evaluated `p == q` in a boolean context, which threw a `TypeError` for
    # FastDifferentiation, whose `==` returns a `Node`.
    import FastDifferentiation
    FastDifferentiation.@variables u1 u2 u3 u4
    Qn = quaternion(u1, u2, u3, u4)
    p = quaternion(1.0, 2.0, 3.0, 4.0)
    @test (Qn / p) isa Quaternion{FastDifferentiation.Node}
    @test (Qn / Qn) isa Quaternion{FastDifferentiation.Node}
    x = [0.3, 0.1, 0.2, 0.4]
    f = FastDifferentiation.make_function(collect(components(Qn / p)), [u1, u2, u3, u4])
    @test f(x) ≈ collect(components(quaternion(x...) / p))
end

@testitem "Cross product extends LinearAlgebra.cross" begin
    using LinearAlgebra
    @test Quaternionic.:× === LinearAlgebra.:×
    @test Quaternionic.:× === LinearAlgebra.cross
    a = quatvec(1.0, 2.0, 3.0)
    b = quatvec(-2.0, 0.5, 1.0)
    # With both `LinearAlgebra` and `Quaternionic` loaded, `×` is not ambiguous.
    @test a × b ≈ quatvec(cross([1.0, 2.0, 3.0], [-2.0, 0.5, 1.0])...)
    @test cross(a, b) == a × b
    @test [1, 0, 0] × [0, 1, 0] == [0, 0, 1]
end

@testitem "Normalized cross product: type stability and null vectors" begin
    for T in (Int, Rational{Int}, BigInt, Float32, Float64)
        a = QuatVec{T}(0, 1, 0, 0)
        @test typeof(a ×̂ QuatVec{T}(0, 2, 0, 0)) === typeof(a ×̂ QuatVec{T}(0, 0, 1, 0))
        @test iszero(a ×̂ QuatVec{T}(0, 2, 0, 0))
        @test a ×̂ QuatVec{T}(0, 0, 3, 0) == quatvec(0, 0, 1)
    end
    @test (@inferred quatvec(1, 0, 0) ×̂ quatvec(2, 0, 0)) isa QuatVecF64
    @test (@inferred quatvec(1//1, 0, 0) ×̂ quatvec(2//1, 0, 0)) isa QuatVecF64

    # A nonzero null vector has zero spinor norm, so it is returned unchanged.
    a = quatvec(1.0, 0, 0)
    b = quatvec(0, 1.0, 1.0im)
    @test a ×̂ b == a × b
    @test !iszero(a ×̂ b)
end

@testitem "Division: divisors whose abs2 overflows or underflows" begin
    # `abs2(p)` overflows or underflows when `abs(p)` is far from 1, and division by it once
    # gave NaN, even for `q / q`.  Now both arguments are rescaled by a power of two first.
    using Random
    rng = Random.MersenneTwister(2468)
    relerr(r, ref) = Float64(abs(big(r) - ref) / abs(ref))
    for T in (Float64, Float32)
        E = exponent(floatmax(T))
        # Each pair of exponents makes `abs2(p)` overflow or underflow, while the quotients
        # stay well inside the range of `T`.  The last pair makes the components of `p`
        # subnormal.
        for (eq, ep) in ((3E÷4, 3E÷4), (-3E÷4, -3E÷4), (3E÷4, E÷2 + 8), (-E÷2 - 8, -3E÷4),
                         (-E - 4, -E - 4))
            @test all(1:200) do _
                q = randn(rng, Quaternion{T}) * ldexp(one(T), eq)
                p = randn(rng, Quaternion{T}) * ldexp(one(T), ep)
                s = 3 * ldexp(one(T), eq)
                components(p / p) == components(one(p)) &&
                    components(quatvec(p) / quatvec(p)) == components(one(p)) &&
                    relerr(q / p, big(q) / big(p)) < 8eps(T) &&
                    relerr(s / p, big(s) / big(p)) < 8eps(T) &&
                    relerr(quatvec(q) / quatvec(p), big(quatvec(q)) / big(quatvec(p))) < 8eps(T)
            end
        end
    end
    q = quaternion(1e200, 2e200, 3e200, 4e200)
    t = quaternion(1e-170, 2e-170, 3e-170, 4e-170)
    @test q / q == one(q)
    @test t / t == one(t)
    @test quaternion(2e200, 0, 0, 0) / quaternion(1e200, 1e200, 0, 0) ≈ quaternion(1, -1, 0, 0)
    @test 1 / quaternion(1e-200) ≈ quaternion(1e200)
end

@testitem "Division: dual numbers are rescaled too" begin
    using ForwardDiff
    W = [0.3, -1.1, 0.7, 0.45]
    for scale in (1e200, 1e-170)
        P = quaternion(3.0, -1.0, 2.0, 5.0) * scale
        x0 = [0.2, 1.3, -0.4, 0.8] * scale
        for f in (x -> sum(W .* components(P / quaternion(x...))),
                  x -> sum(W .* components(quaternion(x...) / P)))
            g = ForwardDiff.gradient(f, x0)
            ref = map(1:4) do i
                h = big(1e-30) * abs(big(x0[i]))
                xp = big.(x0); xp[i] += h
                xm = big.(x0); xm[i] -= h
                (f(xp) - f(xm)) / 2h
            end
            @test all(isfinite, g)
            @test maximum(abs, g - ref) ≤ 4eps() * maximum(abs, ref)
        end
    end
end
