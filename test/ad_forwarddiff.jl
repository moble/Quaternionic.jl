# Tests of `QuaternionicForwardDiffExt`, which extracts values and derivatives from
# quaternion-valued results of ForwardDiff and of DifferentiationInterface's
# `AutoForwardDiff`, and of the derivatives that ForwardDiff computes natively through the
# package's source.  Every reference value is a central finite difference evaluated with
# 256-bit `BigFloat`s, or a function whose derivatives are known to be correct.

@testmodule ForwardDiffTestUtils begin
    using Quaternionic
    using ForwardDiff

    """
        flat(y)

    Return the real numbers stored in `y` as a vector, with complex numbers split into their
    real and imaginary parts.  `y` may be a real or complex number, a quaternion, or an
    array of these.
    """
    flat(y::Real) = [float(y)]
    flat(y::Complex) = [float(real(y)), float(imag(y))]
    flat(y::AbstractQuaternion) = reduce(vcat, map(flat, Tuple(components(y))))
    flat(y::AbstractArray) = reduce(vcat, map(flat, y))

    const weights = (1.0, -0.7, 0.4, 1.3, 0.9, -1.1, 0.6, -0.3)

    """
        sc(y)

    Return a fixed weighted sum of the real numbers in `flat(y)`, which turns a quaternion
    result (with real or complex components) into a scalar loss whose derivatives involve
    every component.
    """
    function sc(y)
        c = flat(y)
        sum(weights[i] * c[i] for i ∈ eachindex(c))
    end

    """
        bigfd(f, x, n=1)

    Return the `n`th derivative of `f` at the real number `x`, for `n` from 1 to 3, as a
    vector of `BigFloat`s that is laid out as `flat(f(x))`.  It is computed by a central
    finite difference with 256-bit precision, whose truncation and rounding errors are both
    far below the precision of `Float64`.
    """
    function bigfd(f, x::Real, n::Integer=1)
        setprecision(BigFloat, 256) do
            X = BigFloat(x)
            F(t) = flat(f(t))
            if n == 1
                h = big"1e-25"
                (F(X + h) - F(X - h)) / 2h
            elseif n == 2
                h = big"1e-18"
                (F(X + h) - 2F(X) + F(X - h)) / h^2
            elseif n == 3
                h = big"1e-14"
                (F(X + 2h) - 2F(X + h) + 2F(X - h) - F(X - 2h)) / (2h^3)
            else
                throw(ArgumentError("Only derivatives of orders 1 through 3 are implemented"))
            end
        end
    end

    """
        bigfdH(f, x)

    Return the Hessian of the scalar function `f` at the real vector `x`, computed by
    central finite differences with 256-bit precision.
    """
    function bigfdH(f, x::AbstractVector)
        setprecision(BigFloat, 256) do
            X = BigFloat.(x)
            h = big"1e-18"
            n = length(X)
            e(i) = [j == i ? h : zero(h) for j ∈ 1:n]
            [
                (f(X + e(i) + e(j)) - f(X + e(i) - e(j))
                 - f(X - e(i) + e(j)) + f(X - e(i) - e(j))) / 4h^2
                for i ∈ 1:n, j ∈ 1:n
            ]
        end
    end

    """
        relerr(a, b)

    Return the largest absolute difference between the entries of `a` and `b`, relative to
    the largest entry of `b` or to 1, whichever is larger, as a `Float64`.
    """
    relerr(a, b) = Float64(maximum(abs, flat(a) - flat(b)) / max(1, maximum(abs, flat(b))))

    """
        nth_derivative(f, n, t)

    Return the `n`th derivative of `f` at `t`, computed by nesting `n` calls to
    `ForwardDiff.derivative`.
    """
    function nth_derivative(f, n, t)
        n == 0 ? f(t) : ForwardDiff.derivative(τ -> nth_derivative(f, n - 1, τ), t)
    end

    """
        hemisphere_sign(R)

    Return the sign, ±1, of the first nonzero component of `R`, which is the sign by which
    `from_rotation_matrix(to_rotation_matrix(R))` differs from `R` itself.
    """
    function hemisphere_sign(R)
        c = first(filter(!iszero, collect(components(R))))
        c < 0 ? -1 : 1
    end
end

@testitem "ForwardDiff: derivatives of quaternion outputs" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    using ForwardDiff
    using .ForwardDiffTestUtils: bigfd, relerr

    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    R₀ = Rotor{Float64}(0.9, 0.3, -0.4, 0.2)  # This rotor is deliberately not normalized.

    # Each case is a function of a real number, the type of its derivative, and the point.
    cases = [
        (t -> Quaternion(sin(t), t^2, exp(t), 1/t), QuaternionF64, 0.7),
        (t -> exp(Quaternion(0.3t, t^2, -t, 0.5)), QuaternionF64, 0.4),
        # A `Rotor` output has a `Quaternion` derivative.  The non-normalizing constructor
        # keeps this one off the unit sphere.
        (t -> Rotor{typeof(t)}(t*R₀[1], R₀[2]*t^2, R₀[3], R₀[4]+t), QuaternionF64, 0.5),
        (t -> exp(t * n̂ / 2), QuaternionF64, 0.3),
        (t -> rotor(1.0, t, -t^2, 0.5), QuaternionF64, 0.6),
        # A `QuatVec` output has a `QuatVec` derivative.
        (t -> t * quatvec(1.0, 2.0, 3.0), QuatVecF64, 0.3),
        (t -> log(exp(t * n̂ / 2)), QuatVecF64, 0.8),
        (t -> quatvec(sin(t), cos(t), t^3), QuatVecF64, 1.1),
        # Complex components, as in `Lorentz` rotors, have complex derivatives.
        (t -> Boost(t, n̂), Quaternion{ComplexF64}, 0.3),
        (t -> Boost(t, n̂), Quaternion{ComplexF64}, 0.0),
        (t -> Boost(t * n̂), Quaternion{ComplexF64}, 0.4),
        (t -> Boost(t, n̂) * exp(t * imz / 3), Quaternion{ComplexF64}, 0.25),
        (t -> Quaternion(complex(t, t^2), 1, im*t, 0), Quaternion{ComplexF64}, 0.5),
    ]
    for (f, D, t) ∈ cases
        d = ForwardDiff.derivative(f, t)
        @test d isa D
        @test relerr(d, bigfd(f, t)) < 1e-14
        # None of these derivatives is zero; before version 4.4.5, those of results with
        # complex components were.
        @test !iszero(d)
    end

    # Constant outputs have zero derivatives, whose types follow the same rules.
    for (y, D) ∈ (
        (Quaternion(1.0, 2.0, 3.0, 4.0), QuaternionF64),
        (rotor(1.0, 2.0, 3.0, 4.0), QuaternionF64),
        (quatvec(2.0, 3.0, 4.0), QuatVecF64),
        (Boost(0.4, n̂), Quaternion{ComplexF64}),
    )
        d = ForwardDiff.derivative(t -> y, 0.3)
        @test d isa D
        @test iszero(d)
    end
    @test iszero(ForwardDiff.derivative(t -> imx, 0.3))

    # This complex result is built directly from complex arithmetic.
    f = t -> Lorentz{typeof(t)}(exp(im * t), 0, 0, im * sin(t))
    @test ForwardDiff.derivative(f, 0.2) isa Quaternion{ComplexF64}
    @test components(ForwardDiff.derivative(f, 0.2)) ≈ [im * exp(0.2im), 0, 0, im * cos(0.2)]

    # Quaternion arguments are not differentiable by ForwardDiff; inputs must be real.
    @test !ForwardDiff.can_dual(QuaternionF64)
    @test !ForwardDiff.can_dual(RotorF64)
end

@testitem "ForwardDiff: nested derivatives of quaternion outputs" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    using ForwardDiff
    using .ForwardDiffTestUtils: bigfd, relerr, nth_derivative

    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    R₀ = rotor(0.4, -0.2, 0.7, 0.5)
    cases = [
        (t -> exp(Quaternion(0.3t, t^2, -t, 0.5)), QuaternionF64, 0.4),
        (t -> Rotor{typeof(t)}(t^3, sin(t), 1, t), QuaternionF64, 0.5),
        (t -> slerp(R₀, exp(t * n̂ / 2) * R₀, 0.3), QuaternionF64, 0.2),
        (t -> exp(t * n̂ / 2)^2.5, QuaternionF64, 0.7),
        (t -> log(exp(t * n̂ / 2)), QuatVecF64, 0.8),
        (t -> quatvec(sin(t), cos(t), t^3), QuatVecF64, 1.1),
        (t -> Boost(t, n̂), Quaternion{ComplexF64}, 0.3),
        (t -> Boost(t * n̂), Quaternion{ComplexF64}, 0.4),
    ]
    for (f, D, t) ∈ cases, n ∈ 2:3
        d = nth_derivative(f, n, t)
        @test d isa D
        @test relerr(d, bigfd(f, t, n)) < 1e-12
    end

    # The inner derivative of a nested call has dual components of the outer tag, and the
    # outer derivative must extract them without confusing the two tags.  Here, the inner
    # derivative is (1, x, 0, 0), so the outer function is x ↦ x * (1, x, 0, 0).
    g = x -> x * ForwardDiff.derivative(y -> Quaternion(x + y, x * y, 0, 0), 1.0)
    @test ForwardDiff.derivative(g, 1.5) isa QuaternionF64
    @test components(ForwardDiff.derivative(g, 1.5)) == [1.0, 3.0, 0.0, 0.0]
    # An inner function that does not depend on its own argument has zero inner derivative,
    # even though its components are duals of the outer tag.
    h = x -> ForwardDiff.derivative(y -> Quaternion(x, x^2, 0, 0), 1.0)
    @test iszero(h(2.0))
    @test iszero(ForwardDiff.derivative(h, 2.0))
    # The same holds with complex components.
    b = x -> x * ForwardDiff.derivative(y -> Boost(x * y, n̂), 0.0)
    @test relerr(ForwardDiff.derivative(b, 0.7), bigfd(b, 0.7)) < 1e-14
end

@testitem "ForwardDiff: value and partials with a tag" tags=[:ad, :forwarddiff] begin
    using ForwardDiff
    using ForwardDiff: Dual, Tag

    # ForwardDiff orders tags by the time of their creation, so `TI` is the newer one, and
    # it may be the tag of the outer layer of a nested dual whose inner layer has tag `TO`.
    outerfunction = x -> x
    innerfunction = x -> 2x
    TO = typeof(Tag(outerfunction, Float64))
    TI = typeof(Tag(innerfunction, Float64))
    d(v, ps...) = Dual{TO}(v, ps...)

    # A `Rotor` keeps its type, and its components are not normalized.
    R = Rotor{Dual{TO,Float64,2}}(d(2.0, 1.0, 0.0), d(0.0, 0.0, 1.0), d(1.0, 2.0, 3.0), d(0.5, 4.0, 5.0))
    @test ForwardDiff.value(TO, R) isa RotorF64
    @test components(ForwardDiff.value(TO, R)) == [2.0, 0.0, 1.0, 0.5]
    @test ForwardDiff.partials(TO, R, 1) isa QuaternionF64
    @test components(ForwardDiff.partials(TO, R, 1)) == [1.0, 0.0, 2.0, 4.0]
    @test components(ForwardDiff.partials(TO, R, 2)) == [0.0, 1.0, 3.0, 5.0]

    Q = Quaternion(d(1.0, 1.0, 2.0), d(2.0, 3.0, 4.0), d(3.0, 5.0, 6.0), d(4.0, 7.0, 8.0))
    @test ForwardDiff.value(TO, Q) === Quaternion(1.0, 2.0, 3.0, 4.0)
    @test ForwardDiff.partials(TO, Q, 2) === Quaternion(2.0, 4.0, 6.0, 8.0)

    V = quatvec(d(2.0, 3.0, 4.0), d(3.0, 5.0, 6.0), d(4.0, 7.0, 8.0))
    @test ForwardDiff.value(TO, V) === quatvec(2.0, 3.0, 4.0)
    @test ForwardDiff.partials(TO, V, 1) === quatvec(3.0, 5.0, 7.0)

    # Complex components are split into real and imaginary parts.
    CD = Complex{Dual{TO,Float64,2}}
    L = Rotor{CD}(complex(d(1.0, 1.0, 0.0), d(0.0, 2.0, 0.0)), complex(d(0.0, 0.0, 1.0), d(0.5, 0.0, 3.0)),
        zero(CD), zero(CD))
    @test ForwardDiff.value(TO, L) isa Lorentz{Float64}
    @test components(ForwardDiff.value(TO, L)) == [1.0 + 0.0im, 0.0 + 0.5im, 0.0im, 0.0im]
    @test ForwardDiff.partials(TO, L, 1) isa Quaternion{ComplexF64}
    @test ForwardDiff.partials(TO, L, 1) == Quaternion(
        [complex(ForwardDiff.partials(real(c), 1), ForwardDiff.partials(imag(c), 1)) for c ∈ components(L)]
    )

    # Constants have their own values and zero partials.
    @test ForwardDiff.value(TO, Quaternion(1.0, 2.0, 3.0, 4.0)) === Quaternion(1.0, 2.0, 3.0, 4.0)
    @test iszero(ForwardDiff.partials(TO, Quaternion(1.0, 2.0, 3.0, 4.0), 1))
    @test ForwardDiff.value(TO, imx) === imx

    # For nested duals, the outer layer is extracted first.
    x = Dual{TI}(d(1.0, 2.0), d(3.0, 4.0))
    y = Dual{TI}(d(5.0, 6.0), d(7.0, 8.0))
    N = quaternion(x, y, zero(x), zero(x))
    @test ForwardDiff.value(TI, N) === quaternion(d(1.0, 2.0), d(5.0, 6.0), d(0.0, 0.0), d(0.0, 0.0))
    @test ForwardDiff.partials(TI, N, 1) === quaternion(d(3.0, 4.0), d(7.0, 8.0), d(0.0, 0.0), d(0.0, 0.0))
    # Asking for the value or partials of the newer tag from duals of the older tag returns
    # the duals themselves or zero, as ForwardDiff does for scalars.
    @test ForwardDiff.value(TI, Q) === Q
    @test iszero(ForwardDiff.partials(TI, Q, 1))
end

@testitem "ForwardDiff: DifferentiationInterface" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    import DifferentiationInterface as DI
    using .ForwardDiffTestUtils: bigfd, relerr

    backend = DI.AutoForwardDiff()
    n̂ = normalize(quatvec(1.0, -2.0, 3.0))
    R₀ = Rotor{Float64}(0.9, 0.3, -0.4, 0.2)

    # Each case is a function of a real number, the types of its value and derivative, and
    # the point.
    cases = [
        (t -> exp(Quaternion(0.3t, t^2, -t, 0.5)), QuaternionF64, QuaternionF64, 0.4),
        (t -> Rotor{typeof(t)}(t*R₀[1], R₀[2]*t^2, R₀[3], R₀[4]+t), RotorF64, QuaternionF64, 0.5),
        (t -> exp(t * n̂ / 2), RotorF64, QuaternionF64, 0.3),
        (t -> quatvec(sin(t), cos(t), t^3), QuatVecF64, QuatVecF64, 1.1),
        (t -> Boost(t, n̂), Lorentz{Float64}, Quaternion{ComplexF64}, 0.3),
        (t -> Boost(t * n̂), Lorentz{Float64}, Quaternion{ComplexF64}, 0.4),
    ]
    for (f, Y, D, t) ∈ cases
        y, d = DI.value_and_derivative(f, backend, t)
        @test y isa Y
        @test components(y) == components(f(t))
        @test d isa D
        @test relerr(d, bigfd(f, t)) < 1e-14
        @test DI.derivative(f, backend, t) == d

        # ForwardDiff propagates several tangents at once, in one batch.
        tangents = (1.0, 2.0, -3.0)
        ty = DI.pushforward(f, backend, t, tangents)
        @test length(ty) == 3
        @test all(tyᵢ isa D for tyᵢ ∈ ty)
        @test all(relerr(tyᵢ, tᵢ * d) < 1e-14 for (tyᵢ, tᵢ) ∈ zip(ty, tangents))
        y2, ty2 = DI.value_and_pushforward(f, backend, t, tangents)
        @test y2 isa Y && components(y2) == components(y)
        @test all(ty2 .== ty)

        # Second derivatives nest two tags.
        d2 = DI.second_derivative(f, backend, t)
        @test d2 isa D
        @test relerr(d2, bigfd(f, t, 2)) < 1e-12
        y3, d3, d23 = DI.value_derivative_and_second_derivative(f, backend, t)
        @test y3 isa Y && components(y3) == components(y)
        @test d3 isa D && relerr(d3, d) < 1e-15
        @test d23 isa D && relerr(d23, d2) < 1e-15
    end

    # Functions of a vector are pushed forward along several tangents.
    vcases = [
        (x -> exp(quatvec(x[1], x[2]*x[3], x[1])), RotorF64, QuaternionF64),
        (x -> Quaternion(x[1], x[2], x[3], x[1]*x[2]) * Quaternion(x[3], 1, x[2], 0), QuaternionF64, QuaternionF64),
        (x -> log(rotor(1.0, x...)), QuatVecF64, QuatVecF64),
        (x -> Boost(x[1], normalize(quatvec(x[2], 1.0, x[3]))), Lorentz{Float64}, Quaternion{ComplexF64}),
    ]
    x = [0.3, 0.2, 0.1]
    tangents = ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.3, -0.5, 0.7])
    for (f, Y, D) ∈ vcases
        y, ty = DI.value_and_pushforward(f, backend, x, tangents)
        @test y isa Y
        @test components(y) == components(f(x))
        for (tyᵢ, dx) ∈ zip(ty, tangents)
            @test tyᵢ isa D
            @test relerr(tyᵢ, bigfd(s -> f(x + s * dx), 0.0)) < 1e-14
        end
    end
end

@testitem "ForwardDiff: Hessians on both sides of series thresholds" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    using ForwardDiff
    import Quaternionic: seriestolerance, trigseriestolerance, isnearscalar, issmallvalue
    using .ForwardDiffTestUtils: sc, bigfdH, relerr

    # Each point is placed 1 percent inside or outside the threshold at which the source
    # switches between a Taylor series and the closed form, and the ForwardDiff Hessian is
    # compared with a `BigFloat` Hessian.  `BigFloat` has much smaller thresholds, so the
    # reference is always computed from the closed form.
    tol = seriestolerance(Float64)
    trigtol = trigseriestolerance(Float64)
    n = vec(normalize(quatvec(1.0, -2.0, 3.0)))
    R₂ = rotor(0.4, -0.2, 0.7, 0.5)
    w = 1.3
    hessian_error(f, x) = relerr(ForwardDiff.hessian(f, x), bigfdH(f, x))

    for (δ, inside) ∈ ((-1e-2, true), (1e-2, false))
        # `exp` switches on `abs2vec(q) ≤ trigtol`.
        x = [0.3; sqrt(trigtol * (1 + δ)) * n]
        @test issmallvalue(abs2vec(Quaternion(x...)), trigseriestolerance) == inside
        @test hessian_error(x -> sc(exp(Quaternion(x...))), x) < 1e-12
        @test hessian_error(x -> sc(exp(quatvec(x[2], x[3], x[4]))), x) < 1e-12

        # `log` and `Rotor ^ s` switch on `abs2vec(q) ≤ tol q[1]²`, through `isnearscalar`.
        x = [w; w * sqrt(tol * (1 + δ)) * n]
        @test isnearscalar(Quaternion(x...)) == inside
        @test hessian_error(x -> sc(log(Quaternion(x...))), x) < 1e-12
        @test hessian_error(x -> sc(log(rotor(x...))), x) < 1e-12
        for s ∈ (0.3, -0.7, 2.5)
            @test hessian_error(x -> sc(rotor(x...)^s), x) < 1e-12
        end

        # Inside the `isnearscalar` branch, `Rotor ^ s` switches again, on y ≤ trigtol, where
        # y = s² atan(√x)² and x = abs2vec(q) / q[1]².  For |s| > 1, this happens at a
        # smaller x than the first threshold.
        for s ∈ (2.5, -3.0)
            xₛ = tan(sqrt(trigtol * (1 + δ)) / abs(s))^2
            x = [w; w * sqrt(xₛ) * n]
            @test isnearscalar(Quaternion(x...))
            @test hessian_error(x -> sc(rotor(x...)^s), x) < 1e-12
        end

        # `distance2` switches on x ≤ ε^(1/4)/2, with x = abs2vec(q) / q[1]² for q = R₁ / R₂.
        dtol = sqrt(sqrt(eps(Float64))) / 2
        θ = atan(sqrt(dtol * (1 + δ)))
        R₁ = Rotor(cos(θ), sin(θ) * n...) * R₂
        q = R₁ / R₂
        @test (abs2vec(q) ≤ dtol * q[1]^2) == inside
        @test hessian_error(x -> distance2(rotor(x...), R₂), collect(components(R₁))) < 1e-12
        @test hessian_error(x -> distance2(R₂, rotor(x...)), collect(components(R₁))) < 1e-12

        # `slerp(R₂, R₁, τ)` takes the power of `R₁ / R₂`, which switches at `tol`.
        θ = atan(sqrt(tol * (1 + δ)))
        R₁ = Rotor(cos(θ), sin(θ) * n...) * R₂
        @test isnearscalar(R₁ / R₂) == inside
        for τ ∈ (0.3, 0.8)
            @test hessian_error(x -> sc(slerp(R₂, rotor(x...), τ)), collect(components(R₁))) < 1e-12
        end
    end
end

@testitem "ForwardDiff: Hessians of from_rotation_matrix" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    using ForwardDiff
    using Random
    using .ForwardDiffTestUtils: sc, bigfdH, relerr, hemisphere_sign

    # `from_rotation_matrix(to_rotation_matrix(rotor(x)))` is `±rotor(x)`, with the sign
    # chosen by the value of `rotor(x)`, so its Hessian is that of `±sc(rotor(x))`.  At an
    # exact rotation, the Bar-Itzhack matrix has a triply degenerate eigenvalue, which made
    # the derivatives of the generic eigendecomposition wrong at second order, and NaN at the
    # identity, before version 4.4.5.
    rng = Random.MersenneTwister(1234)
    points = [
        [0.6, -0.3, 0.5, 0.2],      # generic
        [1.0, 0.0, 0.0, 0.0],       # identity
        [1.0, 1e-9, -2e-9, 0.5e-9], # near the identity
        [0.0, 0.3, -0.5, 0.8],      # pure vector
        [-1.0, 1e-3, 2e-3, -1e-3],  # near -1
        [-1.0, 0.0, 0.0, 0.0],      # exactly -1
        [2.0, -0.6, 1.0, 0.4],      # not normalized
        [randn(rng, 4) for _ ∈ 1:3]...,
    ]
    for x ∈ points
        s = hemisphere_sign(rotor(x...))
        f = x -> sc(from_rotation_matrix(to_rotation_matrix(rotor(x...))))
        g = x -> s * sc(rotor(x...))
        H = ForwardDiff.hessian(f, x)
        G = ForwardDiff.hessian(g, x)
        @test all(isfinite, H)
        @test relerr(H, G) < 1e-12
        @test relerr(G, bigfdH(g, x)) < 1e-12
    end
end

@testitem "ForwardDiff: third derivatives of from_rotation_matrix" tags=[:ad, :forwarddiff] setup=[ForwardDiffTestUtils] begin
    using ForwardDiff
    using .ForwardDiffTestUtils: bigfd, relerr, nth_derivative, hemisphere_sign

    # Derivatives are compared along a path through each point, through third order.  The
    # eigenvector is refined by Newton steps from the eigenvector of the values; one step
    # would give the correct first and second derivatives, but a third derivative that is
    # wrong by order 1.
    direction = [0.3, -0.5, 0.7, 0.2]
    points = [
        [0.6, -0.3, 0.5, 0.2],
        [1.0, 0.0, 0.0, 0.0],
        [1.0, 1e-9, -2e-9, 0.5e-9],
        [0.0, 0.3, -0.5, 0.8],
        [-1.0, 1e-3, 2e-3, -1e-3],
        [2.0, -0.6, 1.0, 0.4],
    ]
    for x ∈ points
        s = hemisphere_sign(rotor(x...))
        f = t -> from_rotation_matrix(to_rotation_matrix(rotor((x + t * direction)...)))
        g = t -> s * rotor((x + t * direction)...)
        for n ∈ 1:3
            fₙ = nth_derivative(f, n, 0.0)
            @test fₙ isa QuaternionF64
            @test relerr(fₙ, nth_derivative(g, n, 0.0)) < 1e-12
            @test relerr(fₙ, bigfd(g, 0.0, n)) < 1e-12
        end
    end
end
