@testitem "Differentiation Interface" begin
    using DifferentiationInterface, DifferentiationInterfaceTest
    import Enzyme, FastDifferentiation,
        FiniteDifferences, ForwardDiff, Mooncake,
        ReverseDiff, Zygote#, ChainRules, ChainRulesCore
    using Random
    Random.seed!(42)
    backends = [
        AutoEnzyme(),
        # AutoFastDifferentiation(),  # has to wait for v0.4 support
        AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(3,1)),
        AutoForwardDiff(),
        AutoMooncake(config=nothing), # Fails on SArray mutation without `Vector` conversion
        AutoReverseDiff(),
        AutoZygote(),  # Fails with return type Vector{QuaternionF64}, where each element
        #                  # is the correct quaternion derivative with the given component
        #                  # cycled to the first component.
        # AutoChainRules(Zygote.ZygoteRuleConfig()),  # Same as above
    ]
    v⃗ = randn(QuatVecF64)
    v̂ = v⃗ / absvec(v⃗)
    f(θ) = Vector(components(exp(θ[1] * v̂)))
    ∇f(θ) = Vector(components(-sin(θ[1]) + cos(θ[1]) * v̂))
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
