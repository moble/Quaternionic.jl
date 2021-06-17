@testset verbose=true "Distance" begin
    @testset "Int" begin
        scalars = Int[0, 1, -1]
        quaternions = [Quaternion(a, b, c, d) for a in scalars for b in scalars for c in scalars for d in scalars]
        d(q₁, q₂) = distance2(q₁, q₂)
        for q₁ in quaternions
            for q₂ in quaternions
                d₁₁ = d(q₁, q₁)
                d₁₂ = d(q₁, q₂)
                d₂₂ = d(q₂, q₂)
                d₂₁ = d(q₂, q₁)
                # Real-valued:
                @test d₁₁ isa Int
                @test d₁₂ isa Int
                @test d₂₁ isa Int
                @test d₂₂ isa Int
                # Symmetry:
                @test d₁₂ == d₂₁
                # Invariance:
                for q₃ in [Quaternion(1), 1imx, 1imy, 1imz]
                    @test d₁₂ == d(q₃*q₁, q₃*q₂)
                    @test d₁₂ == d(q₁*q₃, q₂*q₃)
                end
                # Identity:
                @test d₁₁ == 0
                @test d₂₂ == 0
                # Positive-definiteness
                @test (q₁==q₂) ⊻ (d₁₂>0)
            end
        end
    end
    @testset "$T" for T in [BigFloat, Float64]  # FloatTypes
        scalars = [zero(T), one(T), -one(T)]
        quaternions = [a + b*imn for a in scalars for b in scalars for imn in [imx, imy, imz]]
        quaternions = [quaternions...; randn(Quaternion{T}, 23)...]
        # quaternions = [Quaternion(a, b, c, d) for a in scalars for b in scalars for c in scalars for d in scalars]
        # quaternions = [quaternions...; randn(Quaternion{T}, 19)...]
        @testset "rotor=$rotor" verbose=true for rotor in [true, false]
            @testset "sqrt=$return_sqrt" verbose=true for return_sqrt in [true, false]
                d(q₁, q₂) = return_sqrt ? distance(q₁, q₂) : distance2(q₁, q₂)
                for (i, q₁) in enumerate(quaternions)
                    if rotor
                        if iszero(q₁)
                            continue
                        end
                        q₁ = Rotor(q₁)
                    end
                    for q₂ in quaternions[i+1:end]
                        if rotor
                            if iszero(q₂)
                                continue
                            end
                            q₂ = Rotor(q₂)
                        end
                        absq₁, absq₂ = abs(q₁), abs(q₂)
                        rtol = 10eps(T)
                        atol = 10eps(max(absq₁, absq₂))
                        d₁₁ = d(q₁, q₁)
                        d₁₂ = d(q₁, q₂)
                        d₂₁ = d(q₂, q₁)
                        d₂₂ = d(q₂, q₂)
                        @test d₁₁ isa T
                        @test d₁₂ isa T
                        @test d₂₁ isa T
                        @test d₂₂ isa T
                        @test d₁₂ ≈ d₂₁ rtol=rtol atol=atol
                        for q₃ in quaternions
                            if iszero(q₃)
                                continue
                            end
                            if rotor
                                q₃ = Rotor(q₃)
                            else
                                q₃ = q₃ / abs(q₃)
                            end
                            @test d₁₂ ≈ d(q₃*q₁, q₃*q₂) rtol=rtol atol=atol
                            @test d₁₂ ≈ d(q₁*q₃, q₂*q₃) rtol=rtol atol=atol
                        end
                        @test d₁₁ ≈ zero(T) rtol=rtol atol=atol
                        @test d₂₂ ≈ zero(T) rtol=rtol atol=atol
                        if rotor
                            scalar_multiples = (≈(q₁, q₂, rtol=rtol, atol=atol) || ≈(q₁, -q₂, rtol=rtol, atol=atol))
                            @test scalar_multiples ⊻ (d₁₂ > atol)
                        else
                            @test (q₁==q₂) ⊻ (d₁₂>0)
                        end
                    end
                end
            end
        end
    end
end
