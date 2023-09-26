@testset verbose=true "Auto Diff" begin
    # Make sure everything makes sense to ChainRulesCore
    test_method_tables()

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
        n = √(w^2 + x^2 + y^2 + z^2)
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
            @test all(isnothing, ∇)  # Not really sure why it's not the following line...
            # @test abs(∇[1]) < 10eps(T) && all(isnothing, ∇[2:4])
        end
        for (i,f) ∈ enumerate([

        ])
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test abs(∇[1]) < 10eps(T) && all(isnothing, ∇[2:4])
        end
    end

    # @testset "abs2 QuatVec $T" for T ∈ FloatTypes
    #     w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
    #     for (i,f) ∈ enumerate([
    #         (a,b,c,d)->abs2(quatvec(a,b,c,d)),
    #     ])
    #         ∇ = Zygote.gradient(f, w, x, y, z)
    #         @test isnothing(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
    #     end
    #     for (i,f) ∈ enumerate([
    #         (a,b,c,d)->abs2(quatvec([a,b,c,d])),
    #         (a,b,c,d)->abs2(quatvec(@SVector[a,b,c,d])),
    #         (a,b,c,d)->abs2(QuatVec{T}(a,b,c,d)),
    #         (a,b,c,d)->abs2(QuatVec{T}([a,b,c,d])),
    #         (a,b,c,d)->abs2(QuatVec{T}(@SVector[a,b,c,d])),
    #     ])
    #         #@show i T Zygote.gradient(f, w, x, y, z)
    #         @test all(Zygote.gradient(f, w, x, y, z) .≈ (0w, 2x, 2y, 2z))
    #     end
    # end

    # ϵ = 200eps()

    # @testset verbose=true "log" begin
    #     ∂log∂q(q) = [
    #         ForwardDiff.derivative(ϵ->log(q+ϵ), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imx), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imy), 0),
    #         ForwardDiff.derivative(ϵ->log(q+ϵ*imz), 0)
    #     ]
    # end
end
