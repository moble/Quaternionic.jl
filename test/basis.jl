
@testset verbose=true "Basis" begin
    @testset "$T" for T in Types
        # Note that, because `Num` from Symbolics is a weird type, we have to
        # be a little more explicit below than we normally would be.  Also,
        # because of signed zeros in the float types, we have to take the
        # absolute value of the difference before comparing to zero.

        # Define basis elements
        u = Quaternion(one(T), zero(T), zero(T), zero(T))
        i = Quaternion(zero(T), one(T), zero(T), zero(T))
        j = Quaternion(zero(T), zero(T), one(T), zero(T))
        k = Quaternion(zero(T), zero(T), zero(T), one(T))
        basis = [u, i, j, k]

        # Check basis elements
        for (index, element) in enumerate(basis)
            for component in 1:4
                if component == index
                    @test element[component] == one(T)
                else
                    @test element[component] == zero(T)
                end
            end
        end

        # Check "real" part
        @test real(u) == one(T)
        @test real(i) == zero(T)
        @test real(j) == zero(T)
        @test real(k) == zero(T)

        # Standard expressions
        @test abs(u * u - (u)) == zero(T)
        @test abs(i * i - (-u)) == zero(T)
        @test abs(j * j - (-u)) == zero(T)
        @test abs(k * k - (-u)) == zero(T)
        @test abs(i * j * k - (-u)) == zero(T)

        # Full multiplication table
        @test abs(u * u - (u)) == zero(T)
        @test abs(u * i - (i)) == zero(T)
        @test abs(u * j - (j)) == zero(T)
        @test abs(u * k - (k)) == zero(T)
        @test abs(i * u - (i)) == zero(T)
        @test abs(i * i - (-u)) == zero(T)
        @test abs(i * j - (k)) == zero(T)
        @test abs(i * k - (-j)) == zero(T)
        @test abs(j * u - (j)) == zero(T)
        @test abs(j * i - (-k)) == zero(T)
        @test abs(j * j - (-u)) == zero(T)
        @test abs(j * k - (i)) == zero(T)
        @test abs(k * u - (k)) == zero(T)
        @test abs(k * i - (j)) == zero(T)
        @test abs(k * j - (-i)) == zero(T)
        @test abs(k * k - (-u)) == zero(T)
    end

    @testset "float($T)" for T in IntTypes
        @test float(Quaternion{T}) === Quaternion{float(T)}
        @test float(Quaternion{T}(1, 2, 3, 4)) == Quaternion(float(T)(1), float(T)(2), float(T)(3), float(T)(4))
    end

    @testset "show" begin
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
    end
end
