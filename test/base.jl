@testset verbose=true "Base" begin
    @testset "Numbers $T" for T in Types
        # Note that, because `Num` from Symbolics is a weird type, we have to
        # be a little more explicit below than we normally would be.  Also,
        # because of signed zeros in the float types, we have to take the
        # absolute value of the difference before comparing to zero.

        # Define basis elements
        u = Quaternion{T}(1)
        i = Quaternion{T}(𝐢)
        j = Quaternion{T}(𝐣)
        k = Quaternion{T}(𝐤)
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
        @test u.w==one(T) && u.x==zero(T) && u.y==zero(T) && u.z==zero(T)
        @test i.w==zero(T) && i.x==one(T) && i.y==zero(T) && i.z==zero(T)
        @test j.w==zero(T) && j.x==zero(T) && j.y==one(T) && j.z==zero(T)
        @test k.w==zero(T) && k.x==zero(T) && k.y==zero(T) && k.z==one(T)

        # Check equality with constants; note that these are *equal*, but not *identical*
        @test u == one(T) + zero(T)*𝐢 == one(T)
        @test i == 𝐢 == imx
        @test j == 𝐣 == imy
        @test k == 𝐤 == imz

        @test Quaternion{T}(1, 2, 3, 4) == Quaternion{T}(SVector{4, T}(1, 2, 3, 4))
        @test Quaternion{T}(0, 2, 3, 4) == Quaternion{T}(2, 3, 4)
        @test Quaternion{T}(1, 0, 0, 0) == Quaternion{T}(1)
        @test Quaternion(T.([1, 2, 3, 4])...) == Quaternion(SVector{4, T}(1, 2, 3, 4))
        @test Quaternion(T.([0, 2, 3, 4])...) == Quaternion(T(2), T(3), T(4))
        @test Quaternion(T.([1, 0, 0, 0])...) == Quaternion(T(1))
        if !(T<:Integer)
            @test Rotor{T}(1, 2, 3, 4) == Rotor{T}(SVector{4, T}(1, 2, 3, 4)/√T(30))
            @test Rotor{T}(0, 2, 3, 4) == Rotor{T}(2, 3, 4)
            @test Rotor{T}(1, 0, 0, 0) == Rotor{T}(1)
            @test Rotor(T.([1, 2, 3, 4])...) == Rotor(SVector{4, T}(1, 2, 3, 4)/√T(30))
            @test Rotor(T.([0, 2, 3, 4])...) == Rotor(T(2), T(3), T(4))
            @test Rotor(T.([1, 0, 0, 0])...) == Rotor(T(1))
        end
        @test QuatVec{T}(1, 2, 3, 4) == QuatVec{T}(SVector{4, T}(0, 2, 3, 4))
        @test QuatVec{T}(0, 2, 3, 4) == QuatVec{T}(2, 3, 4)
        @test QuatVec{T}(1, 0, 0, 0) == QuatVec{T}(1)
        @test QuatVec(T.([0, 2, 3, 4])...) == QuatVec(SVector{4, T}(0, 2, 3, 4))
        @test QuatVec(T.([0, 2, 3, 4])...) == QuatVec(T(2), T(3), T(4))
        @test QuatVec(T.([0, 0, 0, 0])...) == QuatVec(T(0))

        # Test indexing
        q = Quaternion(T(1), T(2), T(3), T(4))
        @test q[[1, 2]] == [T(1), T(2)]
        @test q[[3, 4]] == [T(3), T(4)]
        @test q[[4, 2]] == [T(4), T(2)]
        @test q[[4, 2, 3]] == [T(4), T(2), T(3)]
        @test q[[4, 2, 3, 1]] == [T(4), T(2), T(3), T(1)]

        # Test copy constructor and self-equality
        @test Quaternion(u) == Quaternion{T}(u) == u
        @test Quaternion(i) == Quaternion{T}(i) == i
        @test Quaternion(j) == Quaternion{T}(j) == j
        @test Quaternion(k) == Quaternion{T}(k) == k
        @test Quaternionic.wrapper(Quaternion(u)) == Quaternion
        @test Quaternionic.wrapper(Rotor(u)) == Rotor
        @test Quaternionic.wrapper(QuatVec(u)) == QuatVec
        @test u == one(T)
        @test one(T) == u
        @test i != one(T)
        @test one(T) != i
        @test j != one(T)
        @test one(T) != j
        @test k != one(T)
        @test one(T) != k
        @test isequal(u, u)
        @test !isequal(u, i)
        @test !isequal(u, j)
        @test !isequal(u, k)
        @test !isequal(i, u)
        @test isequal(i, i)
        @test !isequal(i, j)
        @test !isequal(i, k)
        @test !isequal(j, u)
        @test !isequal(j, i)
        @test isequal(j, j)
        @test !isequal(j, k)
        @test !isequal(k, u)
        @test !isequal(k, i)
        @test !isequal(k, j)
        @test isequal(k, k)

        # Check isapprox
        @test u ≈ one(T)
        @test one(T) ≈ u
        if T !== Num  # Num(1) ≉ Num(2) doesn't work
            @test i ≉ one(T)
            @test one(T) ≉ i
            @test j ≉ one(T)
            @test one(T) ≉ j
            @test k ≉ one(T)
            @test one(T) ≉ k
        end

        # Check non-sensical basis elements are not allowed
        @test_throws DomainError one(QuatVec)
        @test_throws DomainError one(QuatVec{T})
        @test_throws DomainError one(QuatVec{T}(1, 2, 3))
        @test_throws DomainError zero(Rotor)
        @test_throws DomainError zero(Rotor{T})
        @test_throws DomainError zero(Rotor{T}(1))

        # Check "real" part
        @test real(u) == one(T)
        @test real(i) == zero(T)
        @test real(j) == zero(T)
        @test real(k) == zero(T)

        # Check "imag" part
        @test imag(u) == [zero(T), zero(T), zero(T)]
        @test imag(i) == [one(T), zero(T), zero(T)]
        @test imag(j) == [zero(T), one(T), zero(T)]
        @test imag(k) == [zero(T), zero(T), one(T)]

        # Check "isreal"
        @test isreal(u)
        @test !isreal(i)
        @test !isreal(j)
        @test !isreal(k)

        # Check "isinteger"
        if T != Num
            @test isinteger(u)
            @test !isinteger(1.2u)
        end
        @test !isinteger(i)
        @test !isinteger(j)
        @test !isinteger(k)

        if T<:AbstractFloat
            # Check "isnan"
            @test !isnan(u)
            @test !isnan(i)
            @test !isnan(j)
            @test !isnan(k)
            @test isnan(T(NaN) + 0imx)
            @test isnan(T(NaN)imx)
            @test isnan(T(NaN)imy)
            @test isnan(T(NaN)imz)

            # Check "isinf"
            @test !isinf(u)
            @test !isinf(i)
            @test !isinf(j)
            @test !isinf(k)
            @test isinf(T(Inf) + 0imx)
            @test isinf(T(Inf)imx)
            @test isinf(T(Inf)imy)
            @test isinf(T(Inf)imz)
        end

        # Check "isone"
        @test isone(u)
        @test !isone(2u)
        @test !isone(i)
        @test !isone(j)
        @test !isone(k)

        # Check "flipsign"
        @test flipsign(u, 1) == u
        @test flipsign(u, -1) == -u
        @test flipsign(i, 1) == i
        @test flipsign(i, -1) == -i
        @test flipsign(j, 1) == j
        @test flipsign(j, -1) == -j
        @test flipsign(k, 1) == k
        @test flipsign(k, -1) == -k

        if T != Num
            # Check "in"
            @test u ∈ 0:2
            @test u ∉ 2:4
            @test i ∉ 0:2
            @test j ∉ 0:2
            @test k ∉ 0:2
        end

    end

    @testset "bswap" begin
        @test bswap(Quaternion(1)) == Quaternion(2^56)
        @test bswap(1imx) == (2^56)imx
        @test bswap(1imy) == (2^56)imy
        @test bswap(1imz) == (2^56)imz
        @test bswap(bswap(Quaternion(1))) == Quaternion(1)
        @test bswap(bswap(1imx)) == 1imx
        @test bswap(bswap(1imy)) == 1imy
        @test bswap(bswap(1imz)) == 1imz
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

    @testset "io" begin
        io = IOBuffer()

        Base.show(io, MIME("text/plain"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤"
        Base.show(io, MIME("text/plain"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "1 + 2𝐢 + 3𝐣 + 4𝐤"
        Base.show(io, MIME("text/plain"), Quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "a - b + b*c𝐢 + (c*(d^-1))𝐣 + (d + e)𝐤"
        Base.show(io, MIME("text/latex"), Quaternion{Float64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion{Int64}(1, 2, 3, 4))
        @test String(take!(io)) == "\$1 + 2\\,\\mathbf{i} + 3\\,\\mathbf{j} + 4\\,\\mathbf{k}\$"
        Base.show(io, MIME("text/latex"), Quaternion(a-b, b*c, c/d, d+e))
        @test String(take!(io)) == "\$a - b + b c\\,\\mathbf{i} + \\frac{c}{d}\\,\\mathbf{j} + \\left(d + e\\right)\\,\\mathbf{k}\$"

        for T in PrimitiveTypes
            io = IOBuffer()
            q = Quaternion{T}(1, 2, 3, 4)
            write(io, q)
            seekstart(io)
            p = read(io, typeof(q))
            @test q == p
        end
    end

end
