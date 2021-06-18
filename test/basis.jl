# Just a simple set of tests to ensure that the basis of any given quaternionic
# type actually obeys the standard quaternion rules.

@testset verbose=true "Basis" begin
    @testset "$T" for T in Types
        for Q in [Quaternion, Rotor]
            # Define basis elements
            u = Q{T}(1)
            i = Q{T}(𝐢)
            j = Q{T}(𝐣)
            k = Q{T}(𝐤)

            # Standard expressions
            @test u * u ≈ u atol=eps(T)
            @test i * i ≈ -u atol=eps(T)
            @test j * j ≈ -u atol=eps(T)
            @test k * k ≈ -u atol=eps(T)
            @test i * j * k ≈ -u atol=eps(T)

            # Full multiplication table
            @test u * u ≈ u atol=eps(T)
            @test u * i ≈ i atol=eps(T)
            @test u * j ≈ j atol=eps(T)
            @test u * k ≈ k atol=eps(T)
            @test i * u ≈ i atol=eps(T)
            @test i * i ≈ -u atol=eps(T)
            @test i * j ≈ k atol=eps(T)
            @test i * k ≈ -j atol=eps(T)
            @test j * u ≈ j atol=eps(T)
            @test j * i ≈ -k atol=eps(T)
            @test j * j ≈ -u atol=eps(T)
            @test j * k ≈ i atol=eps(T)
            @test k * u ≈ k atol=eps(T)
            @test k * i ≈ j atol=eps(T)
            @test k * j ≈ -i atol=eps(T)
            @test k * k ≈ -u atol=eps(T)
        end
    end
end
