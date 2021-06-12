@testset verbose=true "Conversion" begin
    @testset verbose=true "$T" for T in FloatTypes
        using EllipsisNotation
        using Random
        Random.seed!(4321)

        @testset "Array types" begin
            if isbitstype(T)
                dumb_as_quat_array(f) = mapslices(x->Quaternion{T}(x...), f, dims=(1,))[1, ..]
                function dumb_as_float_array(q)
                    f = Array{Float64}(undef, (4, size(q)...))
                    f[1, ..] = getproperty.(q, :w)
                    f[2, ..] = getproperty.(q, :x)
                    f[3, ..] = getproperty.(q, :y)
                    f[4, ..] = getproperty.(q, :z)
                    f
                end
                for dims in [(), (5,), (4, 3), (2, 3, 4)]
                    q = randn(Quaternion{T}, dims)
                    f = as_float_array(q)
                    @test size(f) == (4, size(q)...)
                    @test all(f .== dumb_as_float_array(q))
                    @test all(q .== dumb_as_quat_array(f))
                    @test all(q .== as_quat_array(f))
                end
                q = randn(Quaternion{T})
                f = as_float_array(q)
                @test f isa Vector{T}
                @test eltype(f) === T
                @test size(f) == (4,)
                @test f == q.components
            end
        end

        @testset "Euler angles" begin
            N = 5000

            random_angles = [2π .* rand(T, 3) .- π for _ in 1:5000]
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
                q = [Quaternion((2 .* rand(T, 4) .- 1)...) for _ in 1:N]
                random_rotors = @. q / abs(q)
            end
            for i in 1:N
                q1 = random_rotors[i]
                q2 = from_euler_angles(to_euler_angles(q1))
                @test q1 ≈ q2 atol=80eps(T)
            end
            q1 = random_rotors
            q2 = from_euler_angles.(to_euler_angles.(random_rotors))
            @test maximum(abs, q1 .- q2) < 80eps(T)
        end

        @testset "Euler phases" begin
            N = 5000
            ϵ = (T === Float16 ? 20eps(T) : 10eps(T))

            dumb_to_euler_phases(α, β, γ) = [exp(α*im), exp(β*im), exp(γ*im)]

            random_angles = [[3(π*rand(T))-π, π*rand(T), 3(π*rand(T))-π] for _ in 1:5000]
            for i in 1:N
                α, β, γ = random_angles[i]
                z1 = dumb_to_euler_phases(α, β, γ)

                q1 = from_euler_angles(α, β, γ)
                q2 = from_euler_phases(z1...)
                @test min(abs(q1-q2), abs(q1+q2)) < ϵ
                q3 = from_euler_phases(z1)
                @test q2 == q3

                z2 = to_euler_phases(q1)
                @test z1[2] ≈ z2[2] atol=ϵ
                if abs(β) < ϵ  # Typically only happens with Float16
                    @test z1[2] ≈ one(T) atol=ϵ
                    @test z2[2] ≈ one(T) atol=ϵ
                    @test z1[1]*z1[3] ≈ z2[1]*z2[3] atol=ϵ
                elseif abs(β-π) < ϵ  # Typically only happens with Float16
                    @test z1[2] ≈ -one(T) atol=ϵ
                    @test z2[2] ≈ -one(T) atol=ϵ
                    @test z1[1]*z1[3] ≈ z2[1]*z2[3] atol=ϵ
                else
                    @test z1[1] ≈ z2[1] atol=ϵ
                    @test z1[3] ≈ z2[3] atol=ϵ
                end

                q4 = from_euler_angles(α, zero(T), γ)
                @test iszero(q4.x)
                @test iszero(q4.y)
                z3 = dumb_to_euler_phases(α, zero(T), γ)
                z4 = to_euler_phases(q4)
                @test z3[2] ≈ one(T) atol=ϵ
                @test z4[2] ≈ one(T) atol=ϵ
                @test z3[1]*z3[3] ≈ z4[1]*z4[3] atol=ϵ

                q6 = from_euler_angles(α, T(π), γ)
                @test q6.w ≈ zero(T) atol=ϵ
                @test q6.z ≈ zero(T) atol=ϵ
                z5 = dumb_to_euler_phases(α, T(π), γ)
                z6 = to_euler_phases(q6)
                @test z5[2] ≈ -one(T) atol=ϵ
                @test z6[2] ≈ -one(T) atol=ϵ
                @test z5[1]*conj(z5[3]) ≈ z6[1]*conj(z6[3]) atol=(T === Float16 ? 20ϵ : ϵ)

                q7 = Quaternion(0, q6.x, q6.y, 0)
                z7 = to_euler_phases(q7)
                @test z7[2] == -one(T)
                @test z7[1]*conj(z7[3]) ≈ exp((α-γ)*im) atol=10eps(T)
            end
        end

        @testset "Spherical coordinates" begin
            N = 5000

            random_angles = [2π .* rand(T, 3) .- π for _ in 1:5000]
            for i in 1:N
                α, β, γ = random_angles[i]
                q1 = from_spherical_coordinates(β, α)
                q2 = exp(imz*α/2) * exp(imy*β/2) * exp(imz*γ/2)
                @test q1*imz*inv(q1) ≈ q2*imz*inv(q2) atol=10eps(T)
            end

            if isbitstype(T)
                random_rotors = randn_rotor(T, N)
            else
                q = [Quaternion((2 .* rand(T, 4) .- 1)...) for _ in 1:N]
                random_rotors = @. q / abs(q)
            end
            for i in 1:N
                q1 = random_rotors[i]
                q2 = from_spherical_coordinates(to_spherical_coordinates(q1))
                @test q1*imz*inv(q1) ≈ q2*imz*inv(q2) atol=150eps(T)
            end
        end

        @testset "Rotation matrices" begin
            if isbitstype(T)
                for _ in 1:5000
                    q1 = randn_rotor(T)
                    q2 = from_rotation_matrix(to_rotation_matrix(q1))
                    @test min(abs(q1-q2), abs(q1+q2)) < 50eps(T)
                end
            end
        end
    end
end
