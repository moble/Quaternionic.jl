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

        for Q in [Quaternion, QuatVec]
            # Define basis elements
            u = Q{T}(1)
            i = Q{T}(𝐢)
            j = Q{T}(𝐣)
            k = Q{T}(𝐤)
            basis = [u, i, j, k]

            # Basic self-addition/subtraction
            @test u + u == 2u
            @test u - u == 0u
            @test i + i == 2i
            @test i - i == 0i
            @test j + j == 2j
            @test j - j == 0j
            @test k + k == 2k
            @test k - k == 0k

            # Full addition/subtraction table
            for i1 in 1:4
                for i2 in 1:4
                    a = zeros(T, 4)
                    s = zeros(T, 4)
                    a[i1] += one(T)
                    a[i2] += one(T)
                    s[i1] += one(T)
                    s[i2] -= one(T)
                    @test basis[i1] + basis[i2] == Q{T}(a...)
                    @test basis[i1] - basis[i2] == Q{T}(s...)
                end
            end

            # Normalization
            for q in basis
                n = normalize(q)
                @test typeof(n) === Q{float(T)}
                if Q === QuatVec && q == u
                    continue
                end
                @test q == n
                @test q == normalize(2n)
            end
        end

        let Q = Rotor
            # Define basis elements
            u = Q{T}(1)
            i = Q{T}(𝐢)
            j = Q{T}(𝐣)
            k = Q{T}(𝐤)
            basis = [u, i, j, k]

            # Normalization
            for q in basis
                n = normalize(q)
                @test typeof(n) === Q{float(T)}
                @test q == n
                @test q == normalize(2n)
                twoq = Q{T}(2*(q.components))  # Won't normalize
                @test q == normalize(twoq)
            end
        end

        let Q = QuatVec
            i = Q{T}(𝐢)
            j = Q{T}(𝐣)
            k = Q{T}(𝐤)

            @test i × i == zero(i)
            @test j × j == zero(j)
            @test k × k == zero(k)
            @test i × j == k
            @test j × k == i
            @test k × i == j
            @test j × i == -k
            @test k × j == -i
            @test i × k == -j

            @test i ×̂ i == zero(i)
            @test j ×̂ j == zero(j)
            @test k ×̂ k == zero(k)
            @test i ×̂ j == k
            @test j ×̂ k == i
            @test k ×̂ i == j
            @test j ×̂ i == -k
            @test k ×̂ j == -i
            @test i ×̂ k == -j

            @test i ×̂ 2j == k
            @test j ×̂ 2k == i
            @test k ×̂ 2i == j
            @test j ×̂ 2i == -k
            @test k ×̂ 2j == -i
            @test i ×̂ 2k == -j

            @test 2i ×̂ j == k
            @test 2j ×̂ k == i
            @test 2k ×̂ i == j
            @test 2j ×̂ i == -k
            @test 2k ×̂ j == -i
            @test 2i ×̂ k == -j
        end

    end
end
