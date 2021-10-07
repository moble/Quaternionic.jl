@testset verbose=true "Alignment" begin
    Random.seed!(1234)
    @testset verbose=true "Align QuatVec{$T}" for T in [Float16, Float32, Float64]
        for N in [1, 2, 3, 4, 5, 10, 20]
            a⃗ = randn(QuatVec{T}, N)
            R = randn(Rotor{T})

            # Test the exact result
            b⃗ = R .* a⃗ .* conj(R)
            R′ = align(a⃗, b⃗)
            if N > 1
                @test distance(R, conj(R′)) < 25eps(T)
            end
            @test maximum(abs, a⃗ - R′ .* b⃗ .* conj(R′)) < 40eps(T)
            @test_throws DimensionMismatch align(a⃗, b⃗[2:end])

            # Uniform weights
            w = 17ones(T, size(a⃗)) / T(3)
            R′′ = align(a⃗, b⃗, w)
            if N > 1
                @test distance(R, conj(R′′)) < 25eps(T)
            end
            @test maximum(abs, a⃗ - R′′ .* b⃗ .* conj(R′′)) < 40eps(T)
            @test_throws DimensionMismatch align(a⃗, b⃗, w[2:end])

            # Perturb b⃗ slightly
            δ = √eps(T)
            b⃗′′′ = [b⃗i + QuatVec((2*(rand(T, 3) .- 1/T(2)) * δ/√T(3))...) for b⃗i in b⃗]
            R′′′ = align(a⃗, b⃗′′′)
            if N > 1
                @test distance(R, conj(R′′′)) < 25δ
            end
            @test maximum(abs, a⃗ - R′′′ .* b⃗′′′ .* conj(R′′′)) < 40δ

            # Change first third, but use weights to ignore
            if N > 3
                N′ = N ÷ 3
                b⃗′ = copy(b⃗)
                w′ = copy(w)
                b⃗′[1:N′] = randn(QuatVec{T}, N′)
                w′[1:N′] .= 0
                R1 = align(a⃗, b⃗, w′)
                R2 = align(a⃗, b⃗′, w′)
                @test distance(R, conj(R1)) < 25eps(T)
                @test distance(R, conj(R2)) < 25eps(T)
                @test maximum(abs, (a⃗ - R1 .* b⃗ .* conj(R1))[N′+1:end]) < 40eps(T)
                @test maximum(abs, (a⃗ - R2 .* b⃗ .* conj(R2))[N′+1:end]) < 40eps(T)
            end
        end
    end

    @testset verbose=true "Align Rotor{$T}" for T in [Float16, Float32, Float64]
        for N in [1, 2, 3, 4, 5, 10, 20]
            A = randn(Rotor{T}, N)
            R = randn(Rotor{T})

            # Test the exact result
            B = R .* A
            R′ = align(A, B)
            if N > 1
                @test distance(R, conj(R′)) < 25eps(T)
            end
            @test maximum(abs, A - R′ .* B) < 40eps(T)
            @test_throws DimensionMismatch align(A, B[2:end])

            # Uniform weights
            w = 17ones(T, size(A)) / T(3)
            R′′ = align(A, B, w)
            if N > 1
                @test distance(R, conj(R′′)) < 25eps(T)
            end
            @test maximum(abs, A - R′′ .* B) < 40eps(T)
            @test_throws DimensionMismatch align(A, B, w[2:end])

            # Perturb B slightly
            δ = √eps(T)
            B′′′ = [Rotor(Bi + Quaternion((2*(rand(T, 4) .- 1/T(2)) * δ/√T(3))...)) for Bi in B]
            R′′′ = align(A, B′′′)
            if N > 1
                @test distance(R, conj(R′′′)) < 25δ
            end
            @test maximum(abs, A - R′′′ .* B′′′) < 40δ

            # Change first third, but use weights to ignore
            if N > 3
                N′ = N ÷ 3
                B′ = copy(B)
                w′ = copy(w)
                B′[1:N′] = randn(Rotor{T}, N′)
                w′[1:N′] .= 0
                R1 = align(A, B, w′)
                R2 = align(A, B′, w′)
                @test distance(R, conj(R1)) < 25eps(T)
                @test distance(R, conj(R2)) < 25eps(T)
                @test maximum(abs, (A - R1 .* B)[N′+1:end]) < 40eps(T)
                @test maximum(abs, (A - R2 .* B)[N′+1:end]) < 40eps(T)
            end
        end
    end
end
