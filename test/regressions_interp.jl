# Regression tests for the analytic gradients of `log` and `exp`, for `slerp∂slerp`, and
# for `unflip` and `squad`.  The references are either high-precision finite differences in
# `BigFloat` or ForwardDiff applied to the plain functions.

@testmodule InterpReference begin
    using Quaternionic

    # The off-shell extension of `log` that `∂log` differentiates, written out explicitly.
    function log_reference(q::Quaternion{BigFloat})
        w = q[1]
        a = sqrt(q[2]^2 + q[3]^2 + q[4]^2)
        f = iszero(a) ? 1 / w : atan(a, w) / a
        log(abs(q)) + f * Quaternion{BigFloat}(0, q[2], q[3], q[4])
    end

    # The off-shell extension of `exp` that `∂exp` differentiates, written out explicitly.
    function exp_reference(q::Quaternion{BigFloat})
        a = sqrt(q[2]^2 + q[3]^2 + q[4]^2)
        g = iszero(a) ? one(a) : sin(a) / a
        exp(q[1]) * (cos(a) + g * Quaternion{BigFloat}(0, q[2], q[3], q[4]))
    end

    # Central finite differences of `F` at `q` along each component, in 1024-bit BigFloat.
    function gradient_reference(F, q)
        setprecision(BigFloat, 1024) do
            h = big"1e-120"
            Q = Quaternion{BigFloat}(components(q)...)
            basis = [Quaternion{BigFloat}((1:4 .== k)...) for k in 1:4]
            [(F(Q + h*b) - F(Q - h*b)) / (2h) for b in basis]
        end
    end

    maxabs(q::AbstractQuaternion) = maximum(abs, components(q))
    maxerr(a, b) = maximum(maxabs(Quaternion{BigFloat}(components(x)...) - y) for (x, y) in zip(a, b))

    # Machine epsilon of the working precision; BigFloat inputs are made with 256 bits.
    working_eps(::Type{T}) where {T} = eps(T)
    working_eps(::Type{BigFloat}) = big(2.0)^-255

    export log_reference, exp_reference, gradient_reference, maxabs, maxerr, working_eps
end


@testitem "∂log and log∂log are accurate near the identity and the antipode" setup=[InterpReference] tags=[:unit, :validation] begin
    n̂ = normalize(QuatVecF64(0, 0.3, -0.5, 0.8))
    for T in (Float64, Float32, BigFloat)
        ϵ = working_eps(T)
        for size in (1, 1e-2, 1e-4, 1e-8, 1e-12, 0), antipode in (false, true)
            # The exact antipode is a singular point of the logarithm.
            antipode && size == 0 && continue
            Z, (l, ∂l), ∂l′ = setprecision(BigFloat, 256) do
                R = exp(T(size) * QuatVec{T}(n̂))
                Z = antipode ? -R : R
                Z, log∂log(Z), ∂log(Z)
            end
            ref = gradient_reference(log_reference, Z)
            scale = maximum(maxabs, ref)
            @test maxerr(∂l, ref) ≤ 10ϵ * scale
            @test maxerr(∂l′, ref) ≤ 10ϵ * scale
            lref = setprecision(BigFloat, 1024) do
                quatvec(log_reference(Quaternion{BigFloat}(components(Z)...)))
            end
            @test l isa QuatVec
            @test maxerr([quaternion(l)], [quaternion(lref)]) ≤ 10ϵ * max(maxabs(lref), floatmin(T))
        end
    end
    # The identity itself
    @test ∂log(rotor(1)) == QuaternionF64[1, imx, imy, imz]
    l, ∂l = log∂log(rotor(1.0))
    @test iszero(l)
    @test ∂l == QuaternionF64[1, imx, imy, imz]
    # The antipode returns the same value as `log`
    @test log∂log(-rotor(1.0))[1] == log(-rotor(1.0))
end


@testitem "∂exp and exp∂exp are accurate near zero" setup=[InterpReference] tags=[:unit, :validation] begin
    n̂ = normalize(QuatVecF64(0, 0.3, -0.5, 0.8))
    for T in (Float64, Float32, BigFloat)
        ϵ = working_eps(T)
        for size in (3, 1, 1e-1, 1e-2, 1e-4, 1e-8, 1e-12, 0)
            Z, (e, ∂e), ∂e′ = setprecision(BigFloat, 256) do
                Z = T(size) * QuatVec{T}(n̂)
                Z, exp∂exp(Z), ∂exp(Z)
            end
            ref = gradient_reference(exp_reference, Z)
            @test maxerr(∂e, ref) ≤ 10ϵ
            @test maxerr(∂e′, ref) ≤ 10ϵ
            eref = setprecision(BigFloat, 1024) do
                exp_reference(Quaternion{BigFloat}(components(Z)...))
            end
            @test e isa Rotor
            @test maxerr([e], [eref]) ≤ 10ϵ
            # The tiny vector part is kept, rather than being rounded to zero
            @test maxerr([quaternion(quatvec(e))], [quaternion(quatvec(eref))]) ≤ 10ϵ * max(maxabs(quatvec(eref)), floatmin(T))
        end
    end
end


@testitem "Nested ForwardDiff through the analytic gradients at the identity" setup=[InterpReference] tags=[:unit, :validation] begin
    using ForwardDiff
    D = ForwardDiff.derivative
    v = QuatVecF64(0, 0.3, -0.5, 0.8)
    c = QuatVecF64(0, 0.1, 0.2, -0.3)

    # The values carry the partials of the input
    @test D(t -> log∂log(exp(t*v))[1], 0.0) ≈ v atol=2eps()
    @test D(t -> exp∂exp(t*v)[1], 0.0) ≈ v atol=2eps()

    # Derivatives of the gradients along curves through the special points, compared with
    # finite differences of the BigFloat finite-difference gradients
    Zcurve(t) = Rotor{typeof(t)}(1 + t/5, 3t/10, -t/2, 4t/5)  # Not normalized: off shell
    for t₀ in (0.0, 1e-12, 1e-8, 1e-4, 0.5)
        d = D(t -> ∂log(Zcurve(t)), t₀)
        ref = setprecision(BigFloat, 1024) do
            h = big"1e-50"
            (gradient_reference(log_reference, Zcurve(big(t₀)+h)) .- gradient_reference(log_reference, Zcurve(big(t₀)-h))) ./ (2h)
        end
        @test maxerr(d, ref) ≤ 10eps()
    end
    for t₀ in (0.0, 1e-12, 1e-8, 1e-4, 0.5), offset in (0, 1)
        d = D(t -> ∂exp(t*v + offset*c), t₀)
        ref = setprecision(BigFloat, 1024) do
            h = big"1e-50"
            vb, cb = QuatVec{BigFloat}(components(v)...), QuatVec{BigFloat}(components(c)...)
            (gradient_reference(exp_reference, (big(t₀)+h)*vb + offset*cb) .- gradient_reference(exp_reference, (big(t₀)-h)*vb + offset*cb)) ./ (2h)
        end
        @test maxerr(d, ref) ≤ 10eps()
    end
end


@testitem "slerp∂slerp for nearly equal rotors" tags=[:unit, :validation] begin
    using ForwardDiff, Random
    D = ForwardDiff.derivative
    Random.seed!(1234)
    one_ = Quaternion{Bool}(true, false, false, false)
    n̂ = normalize(QuatVecF64(0, 0.3, -0.5, 0.8))
    for q₁ in randn(RotorF64, 3), size in (1, 1e-2, 1e-4, 1e-8, 1e-12, 0), τ in (0.0, 0.3, 1.0)
        q₂ = exp(size * n̂) * q₁
        s, ∂s∂q₁, ∂s∂q₂, ∂s∂τ = slerp∂slerp(q₁, q₂, τ)
        @test s isa Rotor
        @test distance(s, slerp(q₁, q₂, τ)) ≤ 4eps()
        for (k, b) in enumerate((one_, imx, imy, imz))
            @test ∂s∂q₁[k] ≈ D(ϵ -> slerp(q₁ + ϵ*b, q₂, τ), 0.0) atol=10eps()
            @test ∂s∂q₂[k] ≈ D(ϵ -> slerp(q₁, q₂ + ϵ*b, τ), 0.0) atol=10eps()
        end
        @test ∂s∂τ ≈ slerp∂slerp∂τ(q₁, q₂, τ)[2] atol=4eps()
    end
    # The value at τ = 0 carries the partials of τ
    q₁, q₂ = randn(RotorF64, 2)
    @test D(τ -> slerp∂slerp(q₁, q₂, τ)[1], 0.0) ≈ D(τ -> slerp(q₁, q₂, τ), 0.0) atol=10eps()
    @test D(τ -> slerp∂slerp(q₁, q₁, τ)[1], 0.0) ≈ zero(QuaternionF64) atol=10eps()
end


@testitem "squad derivatives agree with ForwardDiff for finely sampled data" tags=[:unit, :validation] begin
    using ForwardDiff, Random
    D = ForwardDiff.derivative
    R, ω⃗, Ṙ = precessing_nutating_example()
    for N in (100, 1000)
        tin = collect(range(0, 1, length=N))
        Rin = rotor.(R.(tin))
        tout = sort([tin[1:7:end]; tin[end]; rand(Random.MersenneTwister(N), 50)])
        Rplain = squad(Rin, tin, tout)
        R₁, Ω⃗₁, Ṙ₁ = squad(Rin, tin, tout; compute_angular_velocity=true, compute_derivative=true)
        R₂, Ṙ₂ = squad(Rin, tin, tout; compute_derivative=true)
        R₃, Ω⃗₃ = squad(Rin, tin, tout; compute_angular_velocity=true)
        ṘFD = [D(t -> squad(Rin, tin, t), t) for t in tout]
        Ω⃗FD = [2quatvec(ṘFD[k] / Rplain[k]) for k in eachindex(tout)]
        # Both derivatives are limited by the conditioning of the interpolant: rounding
        # errors of size eps in the input rotors become relative errors of about eps/(|Ω| dt)
        # in Ṙ, with |Ω| ≈ 6e-3.  Before the accuracy of `log∂log` was fixed, the analytic
        # derivative had relative errors of about 1e-6 here.
        scale = maximum(abs, ṘFD)
        @test maximum(abs, Ṙ₁ .- ṘFD) ≤ 1e-13 * N * scale
        @test maximum(abs, Ω⃗₁ .- Ω⃗FD) ≤ 2e-13 * N * scale
        # The rotors agree with those computed without derivatives
        @test maximum(distance.(Rplain, R₁)) ≤ 4eps()
        @test R₁ == R₂ == R₃
        @test Ṙ₁ == Ṙ₂
        @test Ω⃗₁ == Ω⃗₃
    end
end


@testitem "Nested ForwardDiff through squad derivatives at the knots" tags=[:unit, :validation] begin
    using ForwardDiff
    D = ForwardDiff.derivative
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    ts = collect(1.0:6.0)
    for t in (1.0, 1.0 + 1e-9, 2.0, 2.5, 3.0, 6.0)
        @test D(τ -> squad(qs, ts, τ; compute_derivative=true)[2], t) ≈ D(s -> D(τ -> squad(qs, ts, τ), s), t) atol=20eps()
    end
end


@testitem "ForwardDiff of squad at the first and last knots" tags=[:unit, :validation] begin
    using ForwardDiff
    D = ForwardDiff.derivative
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    ts = collect(1.0:6.0)
    # One-sided BigFloat finite differences
    function reference(t, side)
        setprecision(BigFloat, 512) do
            Rb = [Rotor{BigFloat}(components(q)...) for q in qs]
            tb, h = big.(ts), side * big"1e-40"
            T = big(t)
            (-3squad(Rb, tb, T) + 4squad(Rb, tb, T + h) - squad(Rb, tb, T + 2h)) / (2h)
        end
    end
    for (t, side) in ((1.0, 1), (6.0, -1))
        ref = reference(t, side)
        @test D(τ -> squad(qs, ts, τ), t) ≈ ref atol=10eps()
        @test D(τ -> squad(qs, ts, -τ), -t) ≈ -ref atol=10eps()
        @test squad(qs, ts, t; compute_derivative=true)[2] ≈ ref atol=10eps()
    end
end


@testitem "squad time types and range errors" tags=[:unit, :validation] begin
    using ForwardDiff
    import Quaternionic: squad!
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    ts = collect(1.0:6.0)
    # Integer and rational times
    @test squad(qs, 1:6, 2.5) ≈ squad(qs, ts, 2.5)
    @test squad(qs, collect(1:6), [2, 3]) ≈ squad(qs, ts, [2.0, 3.0])
    @test squad(qs, collect(1//1:6//1), 5//2) ≈ squad(qs, ts, 2.5)
    # Out-of-range times give the range error at either end
    @test_throws ErrorException squad(qs, ts, 0.5)
    @test_throws ErrorException squad(qs, ts, 6.5)
    @test_throws ErrorException squad(qs, ts, [0.5, 1.0])
    @test_throws ErrorException squad(qs, ts, [5.5, 6.5])
    @test_throws ArgumentError squad(qs, ts, [2.0, 6.5]; validate=true)
    @test_throws ArgumentError squad(qs, ts, [3.0, 2.0]; validate=true)
    @test_throws ArgumentError squad(qs, [1.0, 2.0, 2.0, 4.0, 5.0, 6.0], 2.5; validate=true)
    @test_throws ArgumentError squad(qs[1:5], ts, 2.5)
    # A repeated first time is allowed
    @test squad(qs, ts, [1.0, 1.0]) == [qs[1], qs[1]]
    @test squad(qs, ts, 1.0) == qs[1]
    # Derivatives with respect to the knot times
    e₃ = (1:6 .== 3)
    d = ForwardDiff.derivative(h -> squad(qs, ts .+ h .* e₃, 2.5), 0.0)
    ref = setprecision(BigFloat, 512) do
        Rb = [Rotor{BigFloat}(components(q)...) for q in qs]
        h = big"1e-40"
        (squad(Rb, big.(ts) .+ h .* e₃, big(2.5)) - squad(Rb, big.(ts) .- h .* e₃, big(2.5))) / (2h)
    end
    @test d ≈ ref atol=10eps()
    # The documented four-argument method of `squad!`
    Rout = Vector{RotorF64}(undef, 3)
    squad!(Rout, qs, ts, [1.5, 2.5, 3.5])
    @test Rout == squad(qs, ts, [1.5, 2.5, 3.5])
end


@testitem "squad_control_points and squad infer concrete types" tags=[:unit, :fast] begin
    import Quaternionic: squad_control_points
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    ts = collect(1.0:6.0)
    for i in 1:6
        @test (@inferred squad_control_points(qs, ts, i)) isa Tuple{RotorF64, RotorF64}
    end
    # The control point past the end extrapolates linearly
    @test squad_control_points(qs, ts, 6)[2] ≈ qs[6] * conj(qs[5]) * qs[6]
    @test (@inferred squad_control_points(qs, 1:6, 3)) isa Tuple{RotorF64, RotorF64}
    if VERSION ≥ v"1.10"
        sq(R, t, τ) = squad(R, t, τ)
        sqΩṘ(R, t, τ) = squad(R, t, τ; compute_angular_velocity=true, compute_derivative=true)
        @test (@inferred sq(qs, ts, [1.5, 2.5])) isa Vector{RotorF64}
        @test (@inferred sq(qs, ts, 2.5)) isa RotorF64
        @test (@inferred sqΩṘ(qs, ts, 2.5)) isa Tuple{RotorF64, QuatVecF64, QuaternionF64}
    end
end


@testitem "unflip edge cases" tags=[:unit, :fast] begin
    # Empty arrays along the unflipped dimension
    @test unflip(QuaternionF64[]) == QuaternionF64[]
    @test size(unflip(Array{QuaternionF64}(undef, 3, 0); dim=2)) == (3, 0)
    @test size(unflip(Array{QuaternionF64}(undef, 0, 3); dim=1)) == (0, 3)
    @test unflip!(QuaternionF64[]) == QuaternionF64[]
    # `unflip!` returns its argument
    q = [imx, -imx, imx, -imx]
    p = unflip!(q)
    @test p === q
    @test q == [imx, imx, imx, imx]
    # Complex quaternions, using the real part of the inner product
    qc = [
        Quaternion{ComplexF64}(1 + 0.1im, 0.2im, 0, 0.1),
        -Quaternion{ComplexF64}(1 + 0.2im, 0.1im, 0, 0.1),
        Quaternion{ComplexF64}(1, 0, 0.3im, 0),
    ]
    expected = [qc[1], -qc[2], qc[3]]
    @test unflip(qc) == expected
    @test unflip!(copy(qc)) == expected
    @test slerp(qc[1], qc[2], 0.3; unflip=true) ≈ slerp(qc[1], -qc[2], 0.3)
end


@testitem "Analytic-gradient docstring examples run" tags=[:unit, :fast] begin
    ∂log∂w, ∂log∂x, ∂log∂y, ∂log∂z = ∂log(randn(RotorF64))
    @test ∂log∂w isa Quaternion
    (q₁, q₂), τ = randn(RotorF64, 2), rand()
    s, ∂s∂q₁, ∂s∂q₂, ∂s∂τ = slerp∂slerp(q₁, q₂, τ)
    @test length(∂s∂q₁) == length(∂s∂q₂) == 4
end


@testitem "log∂log near the antipode with an underflowing vector norm" setup=[InterpReference] tags=[:unit, :validation] begin
    # For vector parts below about 1.5e-154, `abs2vec` underflows to zero, but the logarithm
    # still points along the vector part, as `log` does.
    for s in (1e-100, 1e-170, 1e-200, 1e-300)
        Z = Rotor{Float64}(-1, 3s, -5s, 8s)
        l, ∂l = log∂log(Z)
        @test maxabs(quaternion(l) - quaternion(log(Z))) ≤ 4eps()
        @test all(isfinite, components(quaternion(l)))
        ref = setprecision(BigFloat, 4096) do
            Q = Quaternion{BigFloat}(components(Z)...)
            h = absvec(Q) * big"1e-60"
            basis = [Quaternion{BigFloat}((1:4 .== k)...) for k in 1:4]
            [(log_reference(Q + h*b) - log_reference(Q - h*b)) / (2h) for b in basis]
        end
        for i in 1:4
            @test maxerr([∂l[i]], [ref[i]]) ≤ 10eps() * maxabs(ref[i])
        end
    end
end


@testitem "squad handles tout in any order" tags=[:unit, :fast] begin
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    ts = collect(1.0:6.0)
    touts = [3.5, 1.5, 5.5, 1.0, 6.0, 2.0, 2.0, 4.25]
    one_at_a_time = [squad(qs, ts, t) for t in touts]
    @test squad(qs, ts, touts) == one_at_a_time
    R, Ω, Ṙ = squad(qs, ts, touts; compute_angular_velocity=true, compute_derivative=true)
    @test R ≈ one_at_a_time
    @test Ω ≈ [squad(qs, ts, t; compute_angular_velocity=true)[2] for t in touts]
    @test Ṙ ≈ [squad(qs, ts, t; compute_derivative=true)[2] for t in touts]
    # Out-of-range times are still caught in any position.
    @test_throws ErrorException squad(qs, ts, [3.5, 0.5])
    @test_throws ErrorException squad(qs, ts, [3.5, 6.5, 2.0])
end


@testitem "unflip! on Rotor arrays negates exactly" tags=[:unit, :fast] begin
    qs = [exp(quatvec(0.3k, sin(k), 0.2cos(2k))) for k in 1:6]
    r = [qs[1], -qs[2], qs[3], -qs[4]]
    p = unflip!(copy(r))
    @test p isa Vector{RotorF64}
    @test p == unflip(r)
    @test p == qs[1:4]
    @test p[2] === -r[2]
    @test p[4] === -r[4]
    # The same holds along the second dimension of a matrix.
    m = permutedims(hcat(r, r))
    pm = unflip!(copy(m); dim=2)
    @test pm isa Matrix{RotorF64}
    @test pm == unflip(m; dim=2)
    @test pm[1, :] == qs[1:4]
end
