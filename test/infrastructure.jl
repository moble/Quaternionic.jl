@testset verbose=true "Infastructure" begin
    @test bswap(Quaternion(1)) == Quaternion(2^56)
    @test bswap(1imx) == (2^56)imx
    @test bswap(1imy) == (2^56)imy
    @test bswap(1imz) == (2^56)imz
    @test bswap(bswap(Quaternion(1))) == Quaternion(1)
    @test bswap(bswap(1imx)) == 1imx
    @test bswap(bswap(1imy)) == 1imy
    @test bswap(bswap(1imz)) == 1imz

    @testset "Types" begin
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

    @testset "io" begin
        io = IOBuffer()

        Base.show(io, MIME("text/plain"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤"
        Base.show(io, MIME("text/plain"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "1 + 2𝐢 + 3𝐣 + 4𝐤"
        Base.show(io, MIME("text/plain"), Quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "a - b + b*c𝐢 + {c*(d^-1)}𝐣 + {d + e}𝐤"
        Base.show(io, MIME("text/latex"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1 + 2\\,\\mathbf{i} + 3\\,\\mathbf{j} + 4\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "\$a - b + b c\\,\\mathbf{i} + \\frac{c}{d}\\,\\mathbf{j} + \\left\\{d + e\\right\\}\\,\\mathbf{k}\$"

        for T in PrimitiveTypes
            io = IOBuffer()
            q = Quaternion{T}(1, 2, 3, 4)
            write(io, q)
            seekstart(io)
            p = read(io, typeof(q))
            @test q == p
        end
    end

    @testset "hash" begin
        for T1 in [FloatTypes...; IntTypes...]
            for T2 in [FloatTypes...; IntTypes...]
                q1 = Quaternion{T1}(1, 2, 3, 4)
                q2 = Quaternion{T2}(1, 2, 3, 4)
                @test isequal(q1, q2) && hash(q1)==hash(q2)
            end
        end

        for T1 in FloatTypes
            for T2 in FloatTypes
                q1 = Quaternion(T1(0.0))
                q2 = Quaternion(T2(-0.0))
                @test !isequal(q1, q2) && hash(q1)!=hash(q2)
                q1 = Quaternion(T1(NaN))
                q2 = Quaternion(T2(NaN))
                @test isequal(q1, q2) && hash(q1)==hash(q2)
            end
        end
    end

end
