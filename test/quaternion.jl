@testset verbose=true "Quaternion" begin
    @testset "Quaternion{T}" begin
        for T in Types
            @test Quaternion(T) === Quaternion{T}
            @test Quaternion(Quaternion{T}) === Quaternion{T}
        end

        for TypeGroup in [FloatTypes, IntTypes]
            for i1 in 1:length(TypeGroup)
                T1 = TypeGroup[i1]
                for i2 in i1+1:length(TypeGroup)
                    T2 = TypeGroup[i2]
                    @test promote_rule(T1, T2) === T1
                    @test promote_rule(Quaternion{T1}, T2) === Quaternion{T1}
                    @test promote_rule(Quaternion{T1}, Quaternion{T2}) === Quaternion{T1}
                end
            end
        end

        for TypeGroup in [FloatTypes, IntTypes]
            @test widen(Quaternion{TypeGroup[1]}) === Quaternion{TypeGroup[1]}
            for i in 2:length(TypeGroup)
                T1 = TypeGroup[i-1]
                T2 = TypeGroup[i]
                @test widen(Quaternion{T2}) === Quaternion{T1}
            end
        end

        for T in FloatTypes
            @test float(Quaternion{T}) === Quaternion{T}
            @test float(Quaternion{T}(1, 2, 3, 4)) == Quaternion(float(T)(1), float(T)(2), float(T)(3), float(T)(4))
            @test float(Quaternion{T}(1, 2, 3, 4)) == Quaternion{T}(1, 2, 3, 4)
        end

        for T in IntTypes
            @test float(Quaternion{T}) === Quaternion{float(T)}
            @test float(Quaternion{T}(1, 2, 3, 4)) == Quaternion(float(T)(1), float(T)(2), float(T)(3), float(T)(4))
        end
    end

end
