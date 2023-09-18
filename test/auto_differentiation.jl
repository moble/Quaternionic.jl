@testset verbose=true "Auto Diff" begin
    # Make sure everything makes sense to ChainRulesCore
    #test_method_tables()

    @testset "abs2 Quaternion $T" for T ∈ [Float64]# [FloatTypes; SymbolicTypes]
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(quaternion(a,b,c,d)),
            (a,b,c,d)->abs2(quaternion([a,b,c,d])),
            (a,b,c,d)->abs2(quaternion(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion{T}(a,b,c,d)),
            (a,b,c,d)->abs2(Quaternion{T}([a,b,c,d])),
            (a,b,c,d)->abs2(Quaternion{T}(@SVector[a,b,c,d])),
        ])
            @info "A $i: "
            @info Zygote.gradient(f, w, x, y, z)
            @test all(Zygote.gradient(f, w, x, y, z) .≈ (2w, 2x, 2y, 2z))
            println()
            println()
        end
    end

    @testset "abs2 Rotor $T" for T ∈ [Float64]# FloatTypes
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        for (i,f) ∈ enumerate([
            (a,b,c,d)->abs2(rotor(a,b,c,d)),
            (a,b,c,d)->abs2(rotor([a,b,c,d])),
            #(a,b,c,d)->abs2(rotor(@SVector[a,b,c,d])),
            (a,b,c,d)->abs2(Rotor{T}(a,b,c,d)),
            #(a,b,c,d)->abs2(Rotor{T}([a,b,c,d])),  # Fail
            #(a,b,c,d)->abs2(Rotor{T}(@SVector[a,b,c,d])),  # Fail
        ])
            @info "B $i: "
            @info Zygote.gradient(f, w, x, y, z)
            @test maximum(abs, Zygote.gradient(f, w, x, y, z)) < 10eps(T)
            println()
            println()
        end
    end

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
