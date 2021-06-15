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
        @test u.w==one(T) && u.x==zero(T) && u.y==zero(T) && u.z==zero(T)
        @test i.w==zero(T) && i.x==one(T) && i.y==zero(T) && i.z==zero(T)
        @test j.w==zero(T) && j.x==zero(T) && j.y==one(T) && j.z==zero(T)
        @test k.w==zero(T) && k.x==zero(T) && k.y==zero(T) && k.z==one(T)
        q = Quaternion(T(1), T(2), T(3), T(4))
        @test q[[1, 2]] == [T(1), T(2)]
        @test q[[3, 4]] == [T(3), T(4)]
        @test q[[4, 2]] == [T(4), T(2)]
        @test q[[4, 2, 3]] == [T(4), T(2), T(3)]
        @test q[[4, 2, 3, 1]] == [T(4), T(2), T(3), T(1)]

        # Check equality with constants; note that these are *equal*, but not *identical*
        @test u == one(T) + zero(T)*𝐢 == one(T)
        @test i == 𝐢 == imx
        @test j == 𝐣 == imy
        @test k == 𝐤 == imz

        # Test copy constructor and self-equality
        @test Quaternion(u) == Quaternion{T}(u) == u == Quaternion{T}(:w)
        @test Quaternion(i) == Quaternion{T}(i) == i == Quaternion{T}(:x)
        @test Quaternion(j) == Quaternion{T}(j) == j == Quaternion{T}(:y)
        @test Quaternion(k) == Quaternion{T}(k) == k == Quaternion{T}(:z)
        if T === Float64
            @test u == Quaternion(:w)
            @test i == Quaternion(:x)
            @test j == Quaternion(:y)
            @test k == Quaternion(:z)
        end
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
end
