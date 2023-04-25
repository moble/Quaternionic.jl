@testset verbose=true "Interpolation" begin
    @testset verbose=true "Unflip" begin
        Random.seed!(1234)
        function measure_flips(q)
            flips = 0
            for i in 2:size(q, 1)
                if q[i][1]*q[i-1][1]+q[i][2]*q[i-1][2]+q[i][3]*q[i-1][3]+q[i][4]*q[i-1][4] < 0
                    flips += 1
                end
            end
            flips
        end

        @testset verbose=true "$T" for T in FloatTypes
            N = 98_203
            println("Constructing randn(Quaternion{$T}, $N)")
            flush(stdout)
            q = randn(Quaternion{T}, N);  # 98_303 because of https://github.com/JuliaLang/julia/issues/49501
            println("Finished constructing randn(Quaternion{$T}, $N)")
            flush(stdout)
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

    @testset verbose=true "Slerp" begin
        ϵ = 20eps()
        comp = [0, 1, -1]
        Rs = [
            [Rotor(w, x, y, z) for w in comp for x in comp for y in comp for z in comp][2:end]...;
            randn(RotorF64, 20)...
        ]
        directions = [Rotor(1), Rotor(imx), Rotor(imy), Rotor(imz), Rotor(-imx), Rotor(-imy), Rotor(-imz)]

        # Check extremes
        for Q1 in directions
            @test distance(slerp(Q1, Q1, 0.0), Q1) < ϵ
            @test distance(slerp(Q1, Q1, 1.0), Q1) < ϵ
            @test distance(slerp(Q1, -Q1, 0.0), Q1) < ϵ
            @test distance(slerp(Q1, -Q1, 1.0), Q1) < ϵ
            for Q2 in directions
                @test distance(slerp(Q1, Q2, 0.0), Q1) < ϵ
                @test distance(slerp(Q1, Q2, 1.0), Q2) < ϵ
                @test distance(slerp(Q1, -Q2, 0.0), Q1) < ϵ
                @test distance(slerp(Q1, -Q2, 1.0), -Q2) < ϵ
                @test distance(slerp(Q2, Q1, 0.0), Q2) < ϵ
                @test distance(slerp(Q2, Q1, 1.0), Q1) < ϵ
            end
        end

        # Test simple rotations about each axis
        for Q2 in directions[2:end]
            for t in LinRange(0.0, 1.0, 100)
                @test distance(
                    slerp(Rotor(1), Q2, t),
                    Rotor(cos(π * t / 2) + sin(π * t / 2) * Q2)
                ) < ϵ
            end
        end

        # Test that (slerp of rotated rotors) is (rotated slerp of rotors)
        for R in Rs
            for Q2 in directions[1:end]
                for t in LinRange(0.0, 1.0, 100)
                    @test distance(
                        R * slerp(Rotor(1), Q2, t),
                        slerp(R, R*Q2, t)
                    ) < ϵ
                end
            end
        end

        # Check that unflipping works
        for Q in directions
            for t in LinRange(0.0, 1.0, 10)
                @test distance(slerp(Q, -Q, t, unflip=true), Q) < ϵ
            end
        end

    end

    @testset verbose=true "Squad" begin
        ϵ = 20eps()
        comp = [0, 1, -1]
        Rs = [
            [Rotor(w, x, y, z) for w in comp for x in comp for y in comp for z in comp][2:end]...;
            randn(RotorF64, 20)...
        ]
        directions = [Rotor(1), Rotor(imx), Rotor(imy), Rotor(imz), Rotor(-imx), Rotor(-imy), Rotor(-imz)]

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

        t_in = LinRange(0.0, 1.0, 13)
        t_out = LinRange(0.0, 1.0, 37)
        t_out2 = sort(rand(59))
        # squad interpolated onto the inputs should be the identity
        for R1 in Rs
            for R2 in Rs
                R_in = [slerp(R1, R2, t) for t in t_in]
                @test squad(R_in, t_in, t_in) ≈ R_in atol=ϵ
            end
        end

        # squad should be the same as slerp for linear interpolation
        for R in directions
            R_in = [slerp(Rotor(1), R, t) for t in t_in]
            R_out_squad = squad(R_in, t_in, t_out)
            R_out_slerp = [slerp(Rotor(1), R, t) for t in t_out]
            @test R_out_squad ≈ R_out_slerp atol=ϵ
            R_out_squad = squad(R_in, t_in, t_out2)
            R_out_slerp = [slerp(Rotor(1), R, t) for t in t_out2]
            @test R_out_squad ≈ R_out_slerp atol=ϵ
        end
    end

    @testset verbose=true "Squad derivatives" begin
        R, ω⃗, Ṙ = precessing_nutating_example()
        tin = collect(range(0, 1, length=1000))
        for tout in [tin, tin[498]]
            Rsquad, ω⃗squad, Ṙsquad = squad(
                Rotor.(R.(tin)), tin, tout,
                compute_angular_velocity=true, compute_derivative=true
            )
            ω⃗eval = ω⃗.(tout)
            Ṙeval = Ṙ.(tout)
            @test ω⃗squad ≈ ω⃗eval rtol=1e-6
            @test Ṙsquad ≈ Ṙeval rtol=1e-6
            Rsquad2, ω⃗squad2 = squad(
                Rotor.(R.(tin)), tin, tout,
                compute_angular_velocity=true
            )
            @test Rsquad == Rsquad2
            @test ω⃗squad == ω⃗squad2
            Rsquad3, Ṙsquad3 = squad(
                Rotor.(R.(tin)), tin, tout,
                compute_derivative=true
            )
            @test Rsquad == Rsquad3
            @test Ṙsquad == Ṙsquad3
        end
    end
end
