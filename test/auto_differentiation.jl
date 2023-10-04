@testset verbose=true "Auto Diff" begin
    # Make sure everything makes sense to ChainRulesCore
    test_method_tables()

    @testset "Quaternion $T rrules" for T ∈ (BigFloat, Float64, Float32)
        ϵ = max(√eps(T), 1e-9)  # Because Float64 is used internally by test_rrule
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        sv = @SVector[w,x,y,z]
        v4 = [w,x,y,z]
        v3 = [x,y,z]
        v1 = [w]
        test_rrule(Quaternion{T}, sv; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion, sv; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, sv; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion{T}, v4; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion, v4; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, v4; rtol=ϵ, atol=ϵ)
        #test_rrule(Quaternion{T}, v3)  # This method doesn't exist
        test_rrule(Quaternion, v3; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, v3; rtol=ϵ, atol=ϵ)
        #test_rrule(Quaternion{T}, v1)  # This method doesn't exist
        test_rrule(Quaternion, v1; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, v1; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion{T}, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion{T}, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion{T}, w; rtol=ϵ, atol=ϵ)
        test_rrule(Quaternion, w; rtol=ϵ, atol=ϵ)
        test_rrule(quaternion, w; rtol=ϵ, atol=ϵ)
    end
    @testset "Quaternion $T components" for T ∈ (BigFloat, Float64)
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        f(w,x,y,z) = components(Quaternion(w,x,y,z))
        @test all(
            Zygote.jacobian(f, w, x, y, z)
            .==
            eachrow(Matrix{T}(LinearAlgebra.I, 4, 4))
        )
    end
    @testset "abs2 Quaternion $T" for T ∈ [FloatTypes; SymbolicTypes]
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Quaternion{T}(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion{T}([a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion{T}(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Quaternion{T}(a,b,c,d)),

            (a,b,c,d)->abs2(quaternion(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(quaternion([a,b,c,d])),
            (a,b,c,d)->abs2(quaternion(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(quaternion(a,b,c,d)),

            (a,b,c,d)->abs2(Quaternion(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion([a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Quaternion(a,b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            #@info "A" i f(w, x, y, z) ∇
            @test all(∇ .≈ (2w, 2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Quaternion{T}(b,c,d)),

            (a,b,c,d)->abs2(quaternion(b,c,d)),
            (a,b,c,d)->abs2(quaternion([b,c,d])),

            (a,b,c,d)->abs2(Quaternion(b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test isnothing(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Quaternion{T}(a)),

            (a,b,c,d)->abs2(quaternion(a)),
            (a,b,c,d)->abs2(quaternion([a])),

            (a,b,c,d)->abs2(Quaternion(a)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test all(isnothing, ∇[2:4]) && ∇[1] .≈ 2w
        end
    end

    @testset "Rotor $T rrules" for T ∈ (BigFloat, Float64, Float32)
        ϵ = max(√eps(T), 1e-9)  # Because Float64 is used internally by test_rrule
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        sv = @SVector[w,x,y,z]
        v4 = [w,x,y,z]
        v3 = [x,y,z]
        v1 = [w]
        test_rrule(Rotor{T}, sv; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor, sv; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, sv; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor{T}, v4; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor, v4; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, v4; rtol=ϵ, atol=ϵ)
        #test_rrule(Rotor{T}, v3)  # This method doesn't exist
        test_rrule(Rotor, v3; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, v3; rtol=ϵ, atol=ϵ)
        #test_rrule(Rotor{T}, v1)  # This method doesn't exist
        test_rrule(Rotor, v1; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, v1; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor{T}, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor{T}, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor{T}, w; rtol=ϵ, atol=ϵ)
        test_rrule(Rotor, w; rtol=ϵ, atol=ϵ)
        test_rrule(rotor, w; rtol=ϵ, atol=ϵ)
    end
    @testset "Rotor $T components" for T ∈ (BigFloat, Float64)
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        f(w,x,y,z) = components(Rotor(w,x,y,z))
        n = √(w^2 + x^2 + y^2 + z^2)
        @test all(
            Zygote.jacobian(f, w, x, y, z)
            .≈
            (
                [(1 - w^2/n^2) / n, - w*x/n^3, - w*y/n^3, - w*z/n^3],
                [- x*w/n^3,  (1 - x^2/n^2) / n, - x*y/n^3, - x*z/n^3],
                [- y*w/n^3, - y*x/n^3,  (1 - y^2/n^2) / n, - y*z/n^3],
                [- z*w/n^3, - z*x/n^3, - z*y/n^3,  (1 - z^2/n^2) / n],
            )
        )
    end
    @testset "abs2 Rotor $T" for T ∈ FloatTypes
        # Note that we can define
        #   f(a,b,c,d) = sum(abs2, @SVector[a,b,c,d] / √(a^2 + b^2 + c^2 + d^2))
        #   g(a,b,c,d) = sum(abs2, @SVector[false,b,c,d] / √(b^2 + c^2 + d^2))
        #   h(a,b,c,d) = sum(abs2, @SVector[a,false,false,false] / √(a^2))
        # And evaluate the gradients to give us guidance on what we expect to find
        # when differentiating these functions for `rotor` and `Rotor`:
        #   Zygote.gradient(f, w, x, y, z) = (0.0, 0.0, 0.0, 0.0)
        #   Zygote.gradient(g, w, x, y, z) = (nothing, 0.0, 0.0, 0.0)
        #   Zygote.gradient(h, w, x, y, z) = (0.0, nothing, nothing, nothing)
        # But remember that `Rotor{T}`` is basically the same as Quaternion{T} in
        # this context.
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Rotor{T}(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Rotor{T}([a,b,c,d])),
            (a,b,c,d)->abs2(Rotor{T}(Rotor{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Rotor{T}(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Rotor{T}(a,b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test all(∇ .≈ (2w, 2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(rotor(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(rotor([a,b,c,d])),
            (a,b,c,d)->abs2(rotor(Rotor{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(rotor(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(rotor(a,b,c,d)),

            (a,b,c,d)->abs2(Rotor(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Rotor([a,b,c,d])),
            (a,b,c,d)->abs2(Rotor(Rotor{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Rotor(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(Rotor(a,b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test maximum(abs, ∇) < 10eps(T)
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Rotor{T}(b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test isnothing(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(rotor(b,c,d)),
            (a,b,c,d)->abs2(rotor([b,c,d])),

            (a,b,c,d)->abs2(Rotor(b,c,d)),
            (a,b,c,d)->abs2(Rotor([b,c,d])),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test isnothing(∇[1]) && maximum(abs, ∇[2:4]) < 10eps(T)
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(Rotor{T}(a)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test ∇[1] .≈ 2w && all(isnothing, ∇[2:4])
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(rotor(a)),
            (a,b,c,d)->abs2(rotor([a])),

            (a,b,c,d)->abs2(Rotor(a)),
            (a,b,c,d)->abs2(Rotor([a])),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test abs(∇[1]) < 10eps(T) && all(isnothing, ∇[2:4])
        end
    end

    @testset "QuatVec $T rrules" for T ∈ (BigFloat, Float64, Float32)
        ϵ = max(√eps(T), 1e-9)  # Because Float64 is used internally by test_rrule
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        sv = @SVector[w,x,y,z]
        v4 = [w,x,y,z]
        v3 = [x,y,z]
        v1 = [w]
        test_rrule(QuatVec{T}, sv; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec, sv; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, sv; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec{T}, v4; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec, v4; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, v4; rtol=ϵ, atol=ϵ)
        #test_rrule(QuatVec{T}, v3)  # This method doesn't exist
        test_rrule(QuatVec, v3; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, v3; rtol=ϵ, atol=ϵ)
        #test_rrule(QuatVec{T}, v1)  # This method doesn't exist
        test_rrule(QuatVec, v1; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, v1; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec{T}, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, w, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec{T}, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, x, y, z; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec{T}, w; rtol=ϵ, atol=ϵ)
        test_rrule(QuatVec, w; rtol=ϵ, atol=ϵ)
        test_rrule(quatvec, w; rtol=ϵ, atol=ϵ)
    end
    @testset "QuatVec $T components" for T ∈ (BigFloat, Float64)
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        f(w,x,y,z) = components(QuatVec(w,x,y,z))
        @test all(
            Zygote.jacobian(f, w, x, y, z)
            .==
            (
                T[false, false, false, false],
                T[false, true, false, false],
                T[false, false, true, false],
                T[false, false, false, true],
            )
        )
    end
    @testset "abs2 QuatVec $T" for T ∈ [FloatTypes; SymbolicTypes]
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(QuatVec{T}(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(QuatVec{T}([a,b,c,d])),
            (a,b,c,d)->abs2(QuatVec{T}(Quaternion{T}(@SVector[a,b,c,d]))),
            (a,b,c,d)->abs2(QuatVec{T}(a,b,c,d)),

            (a,b,c,d)->abs2(quatvec(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(quatvec([a,b,c,d])),
            (a,b,c,d)->abs2(quatvec(Quaternion{T}(@SVector[a,b,c,d]))),

            (a,b,c,d)->abs2(QuatVec(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(QuatVec([a,b,c,d])),
            (a,b,c,d)->abs2(QuatVec(Quaternion{T}(@SVector[a,b,c,d]))),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test all(∇ .≈ (0, 2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(quatvec(a,b,c,d)),

            (a,b,c,d)->abs2(QuatVec(a,b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test isnothing(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(QuatVec{T}(b,c,d)),

            (a,b,c,d)->abs2(quatvec(b,c,d)),
            (a,b,c,d)->abs2(quatvec([b,c,d])),

            (a,b,c,d)->abs2(QuatVec(b,c,d)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test isnothing(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(QuatVec{T}(a)),

            (a,b,c,d)->abs2(quatvec([a]))
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test all(isnothing, ∇[2:4]) && ∇[1] .≈ 0
        end
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(quatvec(a)),

            (a,b,c,d)->abs2(QuatVec(a)),
        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test all(isnothing, ∇)
        end
    end

    @testset "exp(QuatVec::$T) rrules" for T ∈ (Float64,)#(BigFloat, Float64, Float32)
        ϵ = max(√eps(T), 1e-9)  # Because Float64 is used internally by test_rrule
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        v⃗ = QuatVec{T}(x,y,z)
        R̄ = Quaternion{T}(x, y, z, w)
        test_rrule(exp, v⃗; rtol=ϵ, atol=ϵ)
        test_rrule(exp, v⃗ ⊢ R̄; rtol=ϵ, atol=ϵ)
    end
    @testset verbose=true "exp $T" for T ∈ FloatTypes
        f(θ) = exp((θ / 2) * 𝐣)(𝐢) ⋅ 𝐢  # should equal cos(θ)
        #f(θ) = (exp((θ / 2) * 𝐣)(𝐢))[2]  # should equal cos(θ)
        @test f(T(0)) ≈ 1
        @test Zygote.gradient(f, T(0))[1] ≈ -sin(T(0)) atol=10eps(T)
        @test f(T(π)/2) ≈ 0 atol=10eps(T)
        @test Zygote.gradient(f, T(π)/2)[1] ≈ -sin(T(π)/2)
        @test f(T(π)) ≈ -1
        @test Zygote.gradient(f, T(π))[1] ≈ -sin(T(π)) atol=10eps(T)
        @test f(3T(π)/2) ≈ 0 atol=10eps(T)
        @test Zygote.gradient(f, 3T(π)/2)[1] ≈ -sin(3T(π)/2)
        @test f(2T(π)) ≈ 1
        @test Zygote.gradient(f, 2T(π))[1] ≈ -sin(2T(π)) atol=10eps(T)
        for θ ∈ LinRange(0, 2T(π), 100)
            @test f(θ) ≈ cos(θ) atol=10eps(T)
            @test Zygote.gradient(f, θ)[1] ≈ -sin(θ) atol=10eps(T)
        end

        # for v̂ ∈ normalize.(randn(QuatVec{T}, 5))
        #     g(t) = exp(t * v̂) # = cos(t) + v̂*sin(t)
        #     for t ∈ LinRange(0, 2T(π), 100)
        #         @test g(t) ≈ cos(t) + v̂*sin(t)
        #         @show Zygote.gradient(g, t)
        #         # @test Zygote.jacobian(g, t)[1] ≈ components(-sin(t) + v̂*cos(t))
        #     end
        # end

        # f(x,y,z) = exp(x*𝐢 + y*𝐣 + z*𝐤)(𝐢)[1:4]
        # @test f(T(0), T(0), T(0)) ≈ T[0, 1, 0, 0]

        # g(x,y,z) = components(exp(x*𝐢 + y*𝐣 + z*𝐤)(𝐢))
        # @info "jacobian starting"
        # @show Zygote.jacobian(g, T(0), T(0), T(0))
    end

    # @testset verbose=true "log" begin
    #     ∂log∂q(q) = [
    #         ForwardDiff.derivative(ϵ->log(q+ϵ), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imx), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imy), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imz), 0)
    #     ]
    # end
end
