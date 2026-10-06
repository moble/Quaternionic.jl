# Tests of how the AD packages handle quaternions at the boundary of the function being
# differentiated: a quaternion-valued final output, or an array of quaternions.  Each case
# must either give the correct derivative or raise an `ArgumentError`; none may silently
# return derivatives of the scalar parts only, or zeros.  Quaternions inside the function
# are tested elsewhere; here they appear only to check that they are unaffected.  The
# reference derivatives come from `ForwardDiff.derivative`, which returns the full
# quaternion derivative of a quaternion-valued function of a real variable.

@testmodule BoundaryCases begin
    using Quaternionic
    const q₀ = quaternion(0.3, -0.5, 0.7, 0.2)
    scalarmultiple(t) = t * q₀
    rotorpath(t) = rotor(1.0, t, 2t, -t)
    quatvecpath(t) = quatvec(t, t^2, -t)
    exppath(t) = exp(t * quatvec(0.3, -0.4, 0.5)) * q₀
    lorentzpath(t) = Lorentz(quaternion(complex(1.0), complex(0.0, t), 0.0, 0.0))
    pair(x) = [x[1] * q₀, x[2] * conj(q₀)]
    # Fifteen inputs, more than ForwardDiff's default chunk size, so that `jacobian` works
    # in chunks
    many(x) = [sum(x) * q₀, x[1] * exp(quatvec(x[2], x[3], x[4]))]
end

@testitem "AD boundaries: ForwardDiff Jacobians of quaternion arrays" tags=[:ad, :forwarddiff] setup=[BoundaryCases] begin
    using ForwardDiff
    using .BoundaryCases: pair, many, q₀
    J = ForwardDiff.jacobian(pair, [1.5, 2.0])
    @test J isa Matrix{QuaternionF64}
    @test J == [q₀ zero(q₀); zero(q₀) conj(q₀)]
    x = collect(range(0.1, 1.5; length=15))
    J = ForwardDiff.jacobian(many, x)
    @test J isa Matrix{QuaternionF64}
    @test size(J) == (2, 15)
    for i ∈ 1:2, j ∈ 1:15
        d = ForwardDiff.derivative(t -> many([k == j ? t : x[k] for k ∈ 1:15])[i], x[j])
        @test J[i, j] ≈ d atol=4eps()
    end
    # A `Rotor` output has `Quaternion` derivatives, and a `QuatVec` output has `QuatVec`
    # derivatives.
    @test ForwardDiff.jacobian(x -> [rotor(x[1], x[2], 0, 0)], [1.0, 2.0]) isa Matrix{QuaternionF64}
    @test ForwardDiff.jacobian(x -> [quatvec(x[1], x[2], 0)], [1.0, 2.0]) isa Matrix{QuatVecF64}
end

@testitem "AD boundaries: ReverseDiff quaternion outputs" tags=[:ad, :reversediff] setup=[BoundaryCases] begin
    using ReverseDiff
    using .BoundaryCases: q₀, pair
    @test_throws ArgumentError ReverseDiff.gradient(x -> x[1] * q₀, [1.5])
    @test_throws ArgumentError ReverseDiff.jacobian(pair, [1.5, 2.0])
    # Quaternions inside the function are unaffected.
    g = ReverseDiff.gradient(x -> abs2(exp(x[1] * q₀)), [0.4])
    @test g[1] ≈ 2 * abs2(exp(0.4q₀)) * real(q₀)
end

@testitem "AD boundaries: Zygote quaternion outputs" tags=[:ad, :zygote] setup=[BoundaryCases] begin
    using Zygote
    using .BoundaryCases: q₀, scalarmultiple, pair
    @test_throws ArgumentError Zygote.gradient(scalarmultiple, 1.5)
    @test_throws ArgumentError Zygote.jacobian(scalarmultiple, 1.5)
    @test_throws ArgumentError Zygote.jacobian(pair, [1.5, 2.0])
    # Quaternions inside the function, and quaternion arguments, are unaffected.
    @test Zygote.gradient(t -> abs2(exp(t * q₀)), 0.4)[1] ≈ 2 * abs2(exp(0.4q₀)) * real(q₀)
    @test Zygote.gradient(abs2, q₀)[1] ≈ 2q₀
    # An explicit pullback gives each component of the derivative.
    _, back = Zygote.pullback(scalarmultiple, 1.5)
    @test [back(e)[1] for e ∈ (quaternion(1.0), 𝐢 + 0.0, 𝐣 + 0.0, 𝐤 + 0.0)] ≈ collect(components(q₀))
end

@testitem "AD boundaries: DifferentiationInterface" tags=[:ad] setup=[BoundaryCases] begin
    import DifferentiationInterface as DI
    using ForwardDiff, Zygote, ReverseDiff
    import Enzyme
    using .BoundaryCases: q₀, scalarmultiple, rotorpath, quatvecpath, exppath, lorentzpath, pair
    paths = (scalarmultiple, rotorpath, quatvecpath, exppath, lorentzpath)
    # Reverse-mode `derivative` of a quaternion-valued function uses one pullback per unit
    # quaternion, so it gives every component.
    for f ∈ paths, backend ∈ (DI.AutoZygote(), DI.AutoEnzyme(mode=Enzyme.Reverse))
        # Enzyme requires the seed of a `QuatVec` output to be a `QuatVec`, but
        # DifferentiationInterface seeds with `oneunit(y)`, which is a `Quaternion`.
        if f === quatvecpath && backend isa DI.AutoEnzyme
            @test_throws Exception DI.derivative(f, backend, 0.3)
            continue
        end
        d = DI.derivative(f, backend, 0.3)
        d₀ = ForwardDiff.derivative(f, 0.3)
        @test typeof(d) == typeof(d₀)
        @test d ≈ d₀ atol=8eps()
    end
    @test DI.derivative(scalarmultiple, DI.AutoForwardDiff(), 0.3) == q₀
    prep = DI.prepare_derivative(rotorpath, DI.AutoZygote(), 0.3)
    @test DI.derivative(rotorpath, prep, DI.AutoZygote(), 0.4) ≈
        ForwardDiff.derivative(rotorpath, 0.4) atol=8eps()
    @test DI.pushforward(scalarmultiple, DI.AutoZygote(), 1.5, (1.0, 2.0)) == (q₀, 2q₀)
    # Forward-mode Jacobians of arrays of quaternions have quaternion entries.
    @test DI.jacobian(pair, DI.AutoForwardDiff(), [1.5, 2.0]) == [q₀ zero(q₀); zero(q₀) conj(q₀)]
    # Reverse-mode Jacobians of arrays of quaternions raise errors.
    for backend ∈ (DI.AutoZygote(), DI.AutoEnzyme(mode=Enzyme.Reverse), DI.AutoReverseDiff())
        @test_throws ArgumentError DI.jacobian(pair, backend, [1.5, 2.0])
    end
end
