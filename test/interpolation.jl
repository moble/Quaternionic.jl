@testset verbose=true "Interpolation" begin
    @testset verbose=true "Unflip" begin
        Random.seed!(1234)
        function measure_flips(q)
            flips = 0
            for i in 2:size(q, 1)
                if q[i].w*q[i-1].w+q[i].x*q[i-1].x+q[i].y*q[i-1].y+q[i].z*q[i-1].z < 0
                    flips += 1
                end
            end
            flips
        end

        @testset verbose=true "$T" for T in FloatTypes
            q = randn(Quaternion{T}, 100_000);
            flips = measure_flips(q)
            @test flips > 0
            p = unflip(q)
            flips = measure_flips(q)
            @test flips > 0
            flips = measure_flips(p)
            @test flips == 0
            unflip!(q)
            flips = measure_flips(q)
            @test flips == 0

            q = randn(Quaternion{T}, 10_000, 3);
            flips = sum(measure_flips(q[:, i]) for i in 1:3)
            @test flips > 0
            p = unflip(q)
            flips = sum(measure_flips(q[:, i]) for i in 1:3)
            @test flips > 0
            flips = sum(measure_flips(p[:, i]) for i in 1:3)
            @test flips == 0
            unflip!(q)
            flips = sum(measure_flips(q[:, i]) for i in 1:3)
            @test flips == 0

            q = randn(Quaternion{T}, 3, 1_000, 3);
            flips = sum(measure_flips(q[i, :, j]) for i in 1:3 for j in 1:3)
            @test flips > 0
            p = unflip(q, dim=2)
            flips = sum(measure_flips(q[i, :, j]) for i in 1:3 for j in 1:3)
            @test flips > 0
            flips = sum(measure_flips(p[i, :, j]) for i in 1:3 for j in 1:3)
            @test flips == 0
            unflip!(q, dim=2)
            flips = sum(measure_flips(q[i, :, j]) for i in 1:3 for j in 1:3)
            @test flips == 0
        end
    end

    @testset verbose=true "Squad" begin
        ω = 0.1
        tin = collect(LinRange(-10, 10, 201))
        Rin = [exp(ω*ti*imz/2) for ti in tin]
        tout = (tin[1:end-1] + tin[2:end]) / 2
        Rout = squad(Rin, tin, tout, validate=true)
        Ravg = [exp(ω*ti*imz/2) for ti in tout]
        @test maximum(distance.(Rout, Ravg)) < 2eps()
        random_signs = [1; rand([-1, 1], length(Rin)-1)...]  # First one must be 1
        Rin_flipped = [Rotor(r*R) for (r, R) in zip(random_signs, Rin)]
        @test Rout == squad(Rin_flipped, tin, tout, validate=true, unflip=true)
        Rout = squad(Rin, tin, tin, validate=true)
        @test maximum(distance.(Rout, Rin)) < 2eps()
        @test distance(Rotor(1), squad(Rin, tin, 0, validate=true)) == 0
        @test distance(Rin[1], squad(Rin, tin, tin[1], validate=true)) == 0
        @test distance(Rin[end], squad(Rin, tin, tin[end], validate=true)) == 0
        @test squad(Rin, tin, Vector{eltype(tin)}()) == Vector{eltype(Rin)}()
    end
end
