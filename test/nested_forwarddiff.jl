# Higher-order derivatives computed by nesting `ForwardDiff.derivative`, evaluated exactly
# at the points where `iszerovalue` switches a function onto its Taylor-series branch.  See
# issue #113: with ForwardDiff ≥ 1.0, `iszero(::Dual)` also checks the partials, so
# `iszerovalue` must strip *every* level of a nested dual before testing for zero.

@testmodule NestedForwardDiff begin
    using ForwardDiff

    """
        nth_derivative(f, n, t)

    Return the `n`th derivative of `f` at `t`, computed by nesting `n` calls to
    `ForwardDiff.derivative`.  The result has the same shape as `f(t)`.
    """
    nth_derivative(f, n, t) =
        n == 0 ? f(t) : ForwardDiff.derivative(τ -> nth_derivative(f, n - 1, τ), t)

    """
        compare_derivatives(f, g; orders=1:4, t=0.0, atol=1e-12)

    Return `true` if the derivatives of orders `orders` of `f` and `g` at `t` all agree to
    within `atol`.  Here, `f` is the function under test, and `g` is a reference expression
    without branches, so that its derivatives are reliable.  Both must return a vector.
    """
    function compare_derivatives(f, g; orders=1:4, t=0.0, atol=1e-12)
        all(orders) do n
            fₙ = nth_derivative(f, n, t)
            gₙ = nth_derivative(g, n, t)
            ok = isapprox(fₙ, gₙ; atol, rtol=0)
            ok || @info "Mismatch in derivative" n fₙ gₙ
            ok
        end
    end
end

@testitem "nested ForwardDiff: value and iszerovalue" tags=[:unit, :fast] begin
    using ForwardDiff
    import Quaternionic: value, iszerovalue
    using ForwardDiff: Dual

    # A nested dual whose value is zero, but with nonzero partials at both levels
    x = Dual{:outer}(Dual{:inner}(0.0, 1.0), Dual{:inner}(1.0, 0.0))
    @test value(x) === 0.0
    @test iszerovalue(x)
    @test iszerovalue(x * imx)
    @test iszerovalue(vec(x * imx))

    # A triply nested dual
    y = Dual{:outermost}(x, x)
    @test value(y) === 0.0
    @test iszerovalue(y)

    # A nested dual whose value is nonzero
    z = Dual{:outer}(Dual{:inner}(2.0, 0.0), Dual{:inner}(0.0, 0.0))
    @test value(z) === 2.0
    @test !iszerovalue(z)
    @test !iszerovalue(z * imx)
end

@testitem "nested ForwardDiff: exp" tags=[:validation, :fast] setup=[NestedForwardDiff] begin
    using .NestedForwardDiff: nth_derivative, compare_derivatives

    # The examples from issue #113
    @test nth_derivative(t -> exp(t * imx / 2)[1], 2, 0.0) ≈ -1/4
    @test nth_derivative(t -> exp(Quaternion(0.0, t, 0.0, 0.0))[1], 2, 0.0) ≈ -1

    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    nx, ny, nz = vec(n̂)
    for s ∈ (0.0, 0.5)
        # exp(::Quaternion)
        @test compare_derivatives(
            t -> components(exp(Quaternion(s, t*nx, t*ny, t*nz))),
            t -> exp(s) * [cos(t), sin(t)*nx, sin(t)*ny, sin(t)*nz],
        )
    end
    # exp(::QuatVec)
    @test compare_derivatives(
        t -> components(exp(t * n̂ / 2)),
        t -> [cos(t/2), sin(t/2)*nx, sin(t/2)*ny, sin(t/2)*nz],
    )
end

@testitem "nested ForwardDiff: log" tags=[:validation, :fast] setup=[NestedForwardDiff] begin
    using .NestedForwardDiff: compare_derivatives

    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    nx, ny, nz = vec(n̂)
    for s ∈ (1.0, 2.0)
        # log(::Quaternion)
        @test compare_derivatives(
            t -> components(log(Quaternion(s, t*nx, t*ny, t*nz))),
            t -> [log(s^2 + t^2)/2, atan(t/s)*nx, atan(t/s)*ny, atan(t/s)*nz],
        )
    end
    # log(::Rotor)
    @test compare_derivatives(
        t -> components(log(Rotor(cos(t), sin(t)*nx, sin(t)*ny, sin(t)*nz))),
        t -> [zero(t), t*nx, t*ny, t*nz],
    )
end

@testitem "nested ForwardDiff: sincu" tags=[:unit, :fast] setup=[NestedForwardDiff] begin
    using .NestedForwardDiff: compare_derivatives
    import Quaternionic: _sincu

    # sin(t)/t = 1 - t²/6 + t⁴/120 - ⋯
    @test compare_derivatives(t -> [_sincu(t)], t -> [1 - t^2/6 + t^4/120])
end

@testitem "nested ForwardDiff: Boost" tags=[:validation, :fast] setup=[NestedForwardDiff] begin
    using .NestedForwardDiff: compare_derivatives

    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    nx, ny, nz = vec(n̂)
    # Orders beyond 4 would exceed the accuracy of the Taylor series in `Boost`
    @test compare_derivatives(
        t -> [f(c) for c ∈ components(Boost(t * n̂)) for f ∈ (real, imag)],
        t -> let ch = cosh(atanh(t)/2), sh = sinh(atanh(t)/2)
            [ch, zero(t), zero(t), sh*nx, zero(t), sh*ny, zero(t), sh*nz]
        end,
    )
end

@testitem "nested ForwardDiff: iszerovalue on complex" tags=[:unit, :fast] begin
    using ForwardDiff
    import Quaternionic: iszerovalue
    using ForwardDiff: Dual

    # Duals inside a `Complex` must be stripped just like bare duals
    x = Dual{:outer}(Dual{:inner}(0.0, 1.0), Dual{:inner}(1.0, 0.0))
    z = Dual{:outer}(Dual{:inner}(2.0, 0.0), Dual{:inner}(0.0, 0.0))
    @test iszerovalue(complex(x, x))
    @test !iszerovalue(complex(x, z))
    @test !iszerovalue(complex(z, x))
    @test iszerovalue(vec(Quaternion(complex(z), complex(x), complex(x), complex(x))))
    @test iszerovalue(complex(ForwardDiff.Dual(0.0, 1.0)))
end

@testitem "nested ForwardDiff: complex exp, log, sqrt" tags=[:validation, :fast] setup=[NestedForwardDiff] begin
    using .NestedForwardDiff: compare_derivatives

    # Real and imaginary parts of every component, as one real vector
    flat(q) = [f(c) for c ∈ components(q) for f ∈ (real, imag)]
    # The quaternion `s + c t n̂`, in which the vector part is complex when `c` is
    qc(s, c, t, n) = Quaternion(complex(s), c*t*n[1], c*t*n[2], c*t*n[3])
    # Reference results written as `A + B n̂`, with the same flattening as `flat`
    ref(A, B, n) = flat(Quaternion(A, B*n[1], B*n[2], B*n[3]))

    # Base's elementary functions of a `Complex` branch on the values of the real and
    # imaginary parts, and some of those branches drop higher-order derivatives of duals
    # whose values are exactly zero; `atanh`, for example, returns `x` itself when `x == 0`.
    # So the references are built from power series, which need only `+` and `*`.  Forty
    # terms are far more than enough for the arguments used here, which are at most 0.3.
    N = 40
    pcos(w) = evalpoly(w^2, ntuple(k -> Float64((-1)^(k-1) / factorial(big(2k-2))), N))
    psin(w) = w * evalpoly(w^2, ntuple(k -> Float64((-1)^(k-1) / factorial(big(2k-1))), N))
    patan(w) = w * evalpoly(w^2, ntuple(k -> (-1)^(k-1) / (2k-1), N))
    plog1p(u) = evalpoly(u, (0.0, ntuple(k -> (-1)^(k-1) / k, N)...))
    psqrt1p(u) = evalpoly(u, ntuple(k -> prod(j -> (1/2 - j) / (j + 1), 0:k-2; init=1.0), N))

    n = vec(normalize(quatvec(1.0, -2.0, 3.0)))
    # Evaluating at t = 0 tests the branches for a vanishing vector part, and t = 0.3 tests
    # the generic branches.  The closed forms below are valid for `real(s) > 0`.
    for s ∈ (2.0 + 0.5im, 0.5 - 1.0im), c ∈ (1.0 + 0.0im, 0.8 + 0.6im), t₀ ∈ (0.0, 0.3)
        # exp(s + c t n̂) = exp(s) [cos(ct) + n̂ sin(ct)]
        @test compare_derivatives(
            t -> flat(exp(qc(s, c, t, n))),
            t -> ref(exp(s) * pcos(c*t), exp(s) * psin(c*t), n);
            orders=0:4, t=t₀, atol=1e-10,
        )
        # log(s + c t n̂) = log(s² + c²t²)/2 + n̂ atan(ct/s)
        @test compare_derivatives(
            t -> flat(log(qc(s, c, t, n))),
            t -> ref(log(s) + plog1p((c*t/s)^2) / 2, patan(c*t/s), n);
            orders=0:4, t=t₀, atol=1e-10,
        )
        # sqrt(s + c t n̂) = a + n̂ ct/(2a), where a = √((s + √(s² + c²t²))/2)
        @test compare_derivatives(
            t -> flat(sqrt(qc(s, c, t, n))),
            t -> let a = sqrt(s) * psqrt1p((psqrt1p((c*t/s)^2) - 1) / 2)
                ref(a, c*t / (2a), n)
            end;
            orders=0:4, t=t₀, atol=1e-10,
        )
    end
    for c ∈ (1.0 + 0.0im, 0.8 + 0.6im), t₀ ∈ (0.0, 0.3)
        # exp(::QuatVec{<:Complex}) and log(::Rotor{<:Complex})
        @test compare_derivatives(
            t -> flat(exp(quatvec(c*t*n[1], c*t*n[2], c*t*n[3]))),
            t -> ref(pcos(c*t), psin(c*t), n);
            orders=0:4, t=t₀, atol=1e-10,
        )
        @test compare_derivatives(
            t -> flat(log(Rotor(pcos(c*t), psin(c*t)*n[1], psin(c*t)*n[2], psin(c*t)*n[3]))),
            t -> ref(zero(c*t), c*t, n);
            orders=0:4, t=t₀, atol=1e-10,
        )
    end
end

@testitem "complex sqrt with a vanishing vector part" tags=[:unit, :fast] begin
    # Where the vector part vanishes and `real(s) ≤ 0`, `sqrt` returns a complex scalar
    for s ∈ (-4.0 + 0.0im, 3.0im, -3.0im, 0.0im, -1.0 + 2.0im)
        q = Quaternion(s, 0.0im, 0.0im, 0.0im)
        r = sqrt(q)
        # `≈` on a complex quaternion would compare its complex-valued `abs`, so compare
        # components instead
        @test components(r) ≈ [sqrt(s), 0, 0, 0]
        @test components(r * r) ≈ components(q)
    end
    # Where `real(s) > 0`, the result squares back to `q` and has zero vector part
    for s ∈ (4.0 + 0.0im, 2.0 + 0.5im, 0.5 - 1.0im)
        q = Quaternion(s, 0.0im, 0.0im, 0.0im)
        r = sqrt(q)
        @test components(r) ≈ [sqrt(s), 0, 0, 0]
        @test iszero(vec(r))
    end
end
