@testset verbose=true "Quaternion" begin
    @test eltype([1.0, imx]) === QuaternionF64
    @test [1.0, imx] == QuaternionF64[
        QuaternionF64(1.0, 0.0, 0.0, 0.0),
        QuaternionF64(0.0, 1.0, 0.0, 0.0)
    ]

    @testset "wrappers(::T)" begin
        for T in FloatTypes
            @test typeof(randn(Rotor{T}) + randn(Rotor{T})) === Quaternion{T}
            @test typeof(randn(Rotor{T}) - randn(Rotor{T})) === Quaternion{T}
            @test typeof(randn(Rotor{T}) * randn(Rotor{T})) === Rotor{T}
            @test typeof(randn(Rotor{T}) / randn(Rotor{T})) === Rotor{T}

            @test typeof(randn(Rotor{T}) * randn(QuatVec{T})) === Quaternion{T}
            @test typeof(randn(QuatVec{T}) * randn(Rotor{T})) === Quaternion{T}
            @test typeof(randn(Rotor{T}) / randn(QuatVec{T})) === Quaternion{T}
            @test typeof(randn(QuatVec{T}) / randn(Rotor{T})) === Quaternion{T}

            @test typeof(randn(Rotor{T}) + randn(QuatVec{T})) === Quaternion{T}
            @test typeof(randn(QuatVec{T}) + randn(Rotor{T})) === Quaternion{T}
            @test typeof(randn(Rotor{T}) - randn(QuatVec{T})) === Quaternion{T}
            @test typeof(randn(QuatVec{T}) - randn(Rotor{T})) === Quaternion{T}

            @test typeof(randn(QuatVec{T}) + randn(QuatVec{T})) === QuatVec{T}
            @test typeof(randn(QuatVec{T}) - randn(QuatVec{T})) === QuatVec{T}
            @test typeof(randn(QuatVec{T}) * randn(QuatVec{T})) === Quaternion{T}
            @test typeof(randn(QuatVec{T}) / randn(QuatVec{T})) === Quaternion{T}
        end
    end

    @testset "$Q{T}" for (Q,q) in ((Quaternion, quaternion), (Rotor,rotor), (QuatVec,quatvec))
        for T in Types
            @test Q(T) === Q{T}
            @test Q(Q{T}) === Q{T}
            @test real(Q{T}) === T
        end

        for TypeGroup in [FloatTypes, IntTypes]
            for i1 in 1:length(TypeGroup)
                T1 = TypeGroup[i1]
                for i2 in i1+1:length(TypeGroup)
                    T2 = TypeGroup[i2]
                    if all(Ti ∈ [BigFloat, Symbolics.Num] for Ti ∈ [T1, T2])
                        continue
                    end
                    @test promote_rule(T1, T2) === T1
                    if Q === QuatVec
                        @test promote_rule(Q{T1}, T2) === Quaternion{T1}
                    else
                        @test promote_rule(Q{T1}, T2) === Q{T1}
                    end
                    @test promote_rule(Q{T1}, Q{T2}) === Q{T1}
                end
            end
        end

        for T1 in Types
            for T2 in Types
                if all(Ti ∈ [BigFloat, Symbolics.Num] for Ti ∈ [T1, T2])
                    continue
                end
                for Q1 in QTypes
                    for Q2 in QTypes
                        Qexpected = Q1===Q2 ? Q1 : Quaternion
                        Texpected = promote_type(T1, T2)
                        expected = Qexpected{Texpected}
                        observed = promote_rule(Q1{T1}, Q2{T2})
                        @test expected === observed
                    end
                end
            end
        end

        for TypeGroup in [FloatTypes, IntTypes]
            @test widen(Q{TypeGroup[1]}) === Q{TypeGroup[1]}
            for i in 2:length(TypeGroup)
                T1 = TypeGroup[i-1]
                T2 = TypeGroup[i]
                @test widen(Q{T2}) === Q{T1}
            end
        end

        for T in FloatTypes
            @test float(Q{T}) === Q{T}
            @test big(Q{T}) === Q{big(T)}
            if Q !== Rotor
                @test float(Q{T}(1, 2, 3, 4)) == Q(float(T)(1), float(T)(2), float(T)(3), float(T)(4))
                @test big(Q{T}(1, 2, 3, 4)) == Q(big(T)(1), big(T)(2), big(T)(3), big(T)(4))
                @test big(Q{T}(1, 2, 3, 4)) == Q{big(T)}(1, 2, 3, 4)
                @test typeof(Q{T}(1, 2, 3, 4) * rand(T)) === Q{T}
                @test typeof(Q{T}(1, 2, 3, 4) / rand(T)) === Q{T}
                @test typeof(rand(T) * Q{T}(1, 2, 3, 4)) === Q{T}
                @test typeof(rand(T) / Q{T}(1, 2, 3, 4)) === Q{T}
            end
            @test float(Q{T}(1, 2, 3, 4)) == Q{T}(1, 2, 3, 4)

            @test typeof(Rotor(Q{T}(1, 2, 3, 4))) === Rotor{T}
            @test typeof(QuatVec(Q{T}(1, 2, 3, 4))) === QuatVec{T}

            @test_throws DimensionMismatch q(T[1,2])
            @test_throws DimensionMismatch Q(T[1,2])
        end

        for T in IntTypes
            @test float(Q{T}) === Q{float(T)}
            if Q !== Rotor
                @test float(Q{T}(1, 2, 3, 4)) == Q(float(T)(1), float(T)(2), float(T)(3), float(T)(4))
            end
            for i in 1:4
                v = zeros(T, 4)
                v[i] = one(T)
                v = SVector(v...)
                @test float(q(v)) == q(float(T).(v))
                @test float(Q{T}(v)) == Q{T}(float(T).(v))
                v = Vector(v)
                @test float(Q{T}(v)) == Q{T}(float(T).(v))
            end
        end

        for T in SymbolicTypes
            let q = Q{T}(1, 2, 3, 4)
                if Q !== Rotor
                    @test typeof(a * q) === Q{T}
                    @test typeof(a / q) === Q{T}
                    @test typeof(q * a) === Q{T}
                    @test typeof(q / a) === Q{T}
                end
            end

            @test typeof(Rotor{T}(a,b,c,d) + Rotor{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(Rotor{T}(a,b,c,d) - Rotor{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(Rotor{T}(a,b,c,d) * Rotor{T}(w,x,y,z)) === Rotor{T}
            @test typeof(Rotor{T}(a,b,c,d) / Rotor{T}(w,x,y,z)) === Rotor{T}

            @test typeof(Rotor{T}(a,b,c,d) * QuatVec{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(QuatVec{T}(a,b,c,d) * Rotor{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(Rotor{T}(a,b,c,d) / QuatVec{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(QuatVec{T}(a,b,c,d) / Rotor{T}(w,x,y,z)) === Quaternion{T}

            @test typeof(Rotor{T}(a,b,c,d) + QuatVec{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(QuatVec{T}(a,b,c,d) + Rotor{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(Rotor{T}(a,b,c,d) - QuatVec{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(QuatVec{T}(a,b,c,d) - Rotor{T}(w,x,y,z)) === Quaternion{T}

            @test typeof(QuatVec{T}(a,b,c,d) + QuatVec{T}(w,x,y,z)) === QuatVec{T}
            @test typeof(QuatVec{T}(a,b,c,d) - QuatVec{T}(w,x,y,z)) === QuatVec{T}
            @test typeof(QuatVec{T}(a,b,c,d) * QuatVec{T}(w,x,y,z)) === Quaternion{T}
            @test typeof(QuatVec{T}(a,b,c,d) / QuatVec{T}(w,x,y,z)) === Quaternion{T}
        end
    end

end
