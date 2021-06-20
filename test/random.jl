@testset verbose=true "Random" begin
    using Random
    Random.seed!(1234)
    rng = MersenneTwister(1234)

    @testset "$BT" for BT in FloatTypes
        for T in [Quaternion{BT}, Rotor{BT}, QuatVec{BT}]
            @test randn(T) isa T
            @test randn(T, 4) isa Vector{T}
            @test randn(T, 2, 3) isa Matrix{T}
            @test randn(T, 2, 3, 4) isa Array{T}
            # @test randn(T, ()) isa T
            @test randn(T, (4,)) isa Vector{T}
            @test randn(T, (2, 3)) isa Matrix{T}
            @test randn(T, (2, 3, 4)) isa Array{T}

            @test randn(rng, T) isa T
            @test randn(rng, T, 4) isa Vector{T}
            @test randn(rng, T, 2, 3) isa Matrix{T}
            @test randn(rng, T, 2, 3, 4) isa Array{T}
            # @test randn(rng, T, ()) isa T
            @test randn(rng, T, (4,)) isa Vector{T}
            @test randn(rng, T, (2, 3)) isa Matrix{T}
            @test randn(rng, T, (2, 3, 4)) isa Array{T}
        end

        N = (BT === BigFloat ? 10_000 : 1_000_000)
        error = maximum(abs, 1.0 .- abs2.(randn(Rotor{BT}, N)))
        @test error ≈ zero(BT) atol=4eps(BT)
    end
end
