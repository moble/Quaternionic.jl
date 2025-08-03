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
        AutoZygote(),  # Fails with return type Vector{QuaternionF64}, where each element
        #                  # is the correct quaternion derivative with the given component
        #                  # cycled to the first component.
        AutoChainRules(Zygote.ZygoteRuleConfig()),  # Same as above
    ]
    v⃗ = randn(QuatVecF64)
    v = absvec(v⃗)
    v̂ = v⃗ / v

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
            θ->Vector(components(exp(θ[1] * v̂))),
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
        for (f, ∇f) in f∇f
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
