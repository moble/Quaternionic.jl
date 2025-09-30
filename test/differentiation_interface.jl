@testitem "Differentiation Interface" begin
    using DifferentiationInterface, DifferentiationInterfaceTest
    import Enzyme, FastDifferentiation,
        FiniteDifferences, ForwardDiff, Mooncake,
        ReverseDiff, Zygote, ChainRules, ChainRulesCore
    using Random
    ChainRulesCore.debug_mode() = true
    Random.seed!(42)

    backends = [
        AutoEnzyme(),
        # AutoFastDifferentiation(),  # has to wait for v0.4 support
        AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(3,1)),
        AutoForwardDiff(),
        AutoMooncake(config=nothing),
        AutoReverseDiff(),
        AutoZygote(),  # Fails with incorrect results when return type is Quaternionic
        AutoChainRules(Zygote.ZygoteRuleConfig()),  # Same as above
    ]
    # backends = [
    #     AutoEnzyme(),
    #     AutoMooncake(config=nothing),
    # ]
    v⃗ = randn(QuatVecF64)
    v = absvec(v⃗)
    v̂ = v⃗ / v

    g(x, y, z) = Vector([x, y, z])

    f∇f = (
        (
            θ->Vector(components(θ[1] * v⃗)),
            θ->Vector(components(v⃗))
        ),
        (
            θ->Vector(components(v⃗ * θ[1])),
            θ->Vector(components(v⃗))
        ),
        (
            θ->Vector(components(θ[1] * v̂)),
            θ->Vector(components(v̂))
        ),
        (
            θ->Vector(components(v̂ * θ[1])),
            θ->Vector(components(v̂))
        ),
        (
            θ->abs2vec(θ[1] * v⃗),
            # θ->θ[1]^2 * abs2vec(v⃗),
            θ->2 * θ[1] * abs2vec(v⃗),
        ),
        (
            θ->absvec(θ[1] * v⃗),
            # θ->θ[1] * absvec(v⃗),
            θ->absvec(v⃗),
        ),
        (
            θ->vec(θ[1] * v⃗)[1],
            # θ->θ[1] * vec(v⃗)[1],
            θ->vec(v⃗)[1]
        ),
        (
            θ->Vector(components(exp(0 + θ[1] * v̂))),
            #θ->Vector(components(cos(θ[1]) + sin(θ[1]) * v̂)),
            θ->Vector(components(-sin(θ[1]) + cos(θ[1]) * v̂))
        ),
        (
            θ->Vector(components(exp(v̂ * θ[1]))),
            θ->Vector(components(-sin(θ[1]) + cos(θ[1]) * v̂))
        ),
        (
            θ->Vector(components(exp(θ[1] * v⃗))),
            θ->Vector(components(-sin(v * θ[1]) + cos(v * θ[1]) * v̂) .* v)
        ),
        (
            θ->Vector(components(exp(v⃗ * θ[1]))),
            θ->Vector(components(-sin(v * θ[1]) + cos(v * θ[1]) * v̂) .* v)
        ),
        # (
        #     θ->Vector(components(exp(1.2θ[1] + 3.4v̂))),
        #     #θ->Vector(components(exp(1.2θ[1])*(cos(3.4) + sin(3.4) * v̂)),
        #     θ->Vector(components(1.2exp(1.2θ[1])*(cos(3.4) + sin(3.4) * v̂)))
        # ),
    )

    x = Float64[
        0.0,
        0.1,
        π/2,
        3.01,
        π,
    ]
    scenarios = [
        Scenario{:derivative,:out}(f, θ; res1=∇f(θ))
        for θ in x
        for (i,(f, ∇f)) in enumerate(f∇f)
        if i≠6 || θ ≠ 0.0  # skip absvec at 0 because it's not differentiable there
    ]
    test_differentiation(
        backends,  # the backends you want to compare
        scenarios,  # the scenarios you defined,
        correctness=true,  # compares values against the reference
        type_stability=:none,  # :prepared or :full checks type stability with JET.jl
        detailed=true,  # prints a detailed test set
        rtol=1e-8,
    )
end
