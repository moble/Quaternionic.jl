@testset verbose=true "Basis" begin
    @testset "$T" for T in Types
        # Note that, because `Num` from Symbolics is a weird type, we have to
        # be a little more explicit below than we normally would be.  Also,
        # because of signed zeros in the float types, we have to take the
        # absolute value of the difference before comparing to zero.

        # Define basis elements
        u = Quaternion{T}(1)
        i = Quaternion{T}(𝐢)
        j = Quaternion{T}(𝐣)
        k = Quaternion{T}(𝐤)

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
