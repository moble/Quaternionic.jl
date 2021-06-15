@testset verbose=true "Interpolation" begin
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

