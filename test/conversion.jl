@testset verbose=true "Conversion" begin
    @testset "$T" for T in FloatTypes
        using Random
        Random.seed!(4321)
        # as_quat_array
        # as_float_array
        @testset "Euler angles" begin
            N = 5000

            random_angles = [2π .* rand(T, 3) .- π for _ in 1:5000];
            for i in 1:N
                α, β, γ = random_angles[i]
                q1 = from_euler_angles(α, β, γ)
                q2 = exp(imz*α/2) * exp(imy*β/2) * exp(imz*γ/2)
                @test q1 ≈ q2 atol=10eps(T)
            end
            q1 = from_euler_angles.(random_angles)
            q2 = broadcast((αβγ)->(exp(αβγ[1]*imz/2)*exp(αβγ[2]*imy/2)*exp(αβγ[3]*imz/2)), random_angles)
            @test maximum(abs, q1 .- q2) < 10eps(T)

            if isbitstype(T)
                random_rotors = randn_rotor(T, N)
            else
                q = as_quat_array(2*rand(T, 4, N) - 1)
                random_rotors = @. q / abs(q)
            end
            for i in 1:N
                q1 = random_rotors[i]
                q2 = from_euler_angles(to_euler_angles(q1))
                @test q1 ≈ q2 atol=10eps(T)
            end
            q1 = random_rotors
            q2 = from_euler_angles.(to_euler_angles.(random_rotors))
            @test q1 .≈ q2 atol=10eps(T)
        end

        # to_euler_phases!(z::Array{Complex{T}}, R::Quaternion{T})
        # to_euler_phases(R::Quaternion{T})
        # from_euler_phases(z)
    end
end
