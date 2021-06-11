@testset verbose=true "Random" begin
    using Random
    Random.seed!(1234)
    rng = MersenneTwister(1234)

    @test randn_rotor() isa QuaternionF64
    @test randn_rotor(4) isa Vector{QuaternionF64}
    @test randn_rotor(2, 3) isa Matrix{QuaternionF64}
    @test randn_rotor(2, 3, 4) isa Array{QuaternionF64}
    @test randn_rotor(()) isa QuaternionF64
    @test randn_rotor((4,)) isa Vector{QuaternionF64}
    @test randn_rotor((2, 3)) isa Matrix{QuaternionF64}
    @test randn_rotor((2, 3, 4)) isa Array{QuaternionF64}

    @test randn_rotor(rng) isa QuaternionF64
    @test randn_rotor(rng, 4) isa Vector{QuaternionF64}
    @test randn_rotor(rng, 2, 3) isa Matrix{QuaternionF64}
    @test randn_rotor(rng, 2, 3, 4) isa Array{QuaternionF64}
    @test randn_rotor(rng, ()) isa QuaternionF64
    @test randn_rotor(rng, (4,)) isa Vector{QuaternionF64}
    @test randn_rotor(rng, (2, 3)) isa Matrix{QuaternionF64}
    @test randn_rotor(rng, (2, 3, 4)) isa Array{QuaternionF64}

    for BT in [Float64, Float32, Float16]
        QT = Quaternion{BT}
        for T in [BT, QT]
            @test randn_rotor(T) isa QT
            @test randn_rotor(T, 4) isa Vector{QT}
            @test randn_rotor(T, 2, 3) isa Matrix{QT}
            @test randn_rotor(T, 2, 3, 4) isa Array{QT}
            @test randn_rotor(T, ()) isa QT
            @test randn_rotor(T, (4,)) isa Vector{QT}
            @test randn_rotor(T, (2, 3)) isa Matrix{QT}
            @test randn_rotor(T, (2, 3, 4)) isa Array{QT}

            @test randn_rotor(rng, T) isa QT
            @test randn_rotor(rng, T, 4) isa Vector{QT}
            @test randn_rotor(rng, T, 2, 3) isa Matrix{QT}
            @test randn_rotor(rng, T, 2, 3, 4) isa Array{QT}
            @test randn_rotor(rng, T, ()) isa QT
            @test randn_rotor(rng, T, (4,)) isa Vector{QT}
            @test randn_rotor(rng, T, (2, 3)) isa Matrix{QT}
            @test randn_rotor(rng, T, (2, 3, 4)) isa Array{QT}
        end
    end

    for T in [Float64, Float32, Float16]
        error = maximum(abs, 1.0 .- abs2.(randn_rotor(Quaternion{T}, 1_000_000)))
        @test error ≈ zero(T) atol=4eps(T)
    end
end
