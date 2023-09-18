@testset verbose=true "Gradients" begin
    ϵ = 200eps()

    @testset verbose=true "log" begin
        ∂log∂q(q) = [
            ForwardDiff.derivative(ϵ->log(q+ϵ), 0),
            ForwardDiff.derivative(ϵ->log(q+ϵ*imx), 0),
            ForwardDiff.derivative(ϵ->log(q+ϵ*imy), 0),
            ForwardDiff.derivative(ϵ->log(q+ϵ*imz), 0)
        ]
        for q in randn(RotorF64, 10_000)
            ∂1 = ∂log(q)
            ∂2 = ∂log∂q(q)
            for i in 1:4
                @test ∂1[i] ≈ ∂2[i] rtol=ϵ
            end
            l, ∂3 = log∂log(q)
            @test l ≈ log(q) rtol=4eps() atol=10eps()
            for i in 1:4
                @test ∂3[i] ≈ ∂2[i] rtol=ϵ
            end
        end
        ∂1 = ∂log(rotor(1))
        ∂2 = QuaternionF64[1, 0, 0, imz]
        for i in 1:4
            @test ∂1[i] ≈ ∂2[i] rtol=ϵ
        end
        l, ∂1 = log∂log(rotor(1))
        ∂2 = QuaternionF64[1, 0, 0, imz]
        @test l ≈ QuaternionF64(0) rtol=4eps() atol=10eps()
        for i in 1:4
            @test ∂1[i] ≈ ∂2[i] rtol=ϵ
        end
    end

    @testset verbose=true "exp" begin
        ∂exp∂q(q) = [
            ForwardDiff.derivative(ϵ->exp(q+(ϵ+0imx)), 0),
            ForwardDiff.derivative(ϵ->exp(q+ϵ*imx), 0),
            ForwardDiff.derivative(ϵ->exp(q+ϵ*imy), 0),
            ForwardDiff.derivative(ϵ->exp(q+ϵ*imz), 0)
        ]
        for q in randn(QuatVecF64, 10_000)
            ∂1 = ∂exp(q)
            ∂2 = ∂exp∂q(q)
            for i in 1:4
                @test ∂1[i] ≈ ∂2[i] rtol=ϵ
            end
            e, ∂3 = exp∂exp(q)
            @test e ≈ exp(q) rtol=4eps()
            for i in 1:4
                @test ∂3[i] ≈ ∂2[i] rtol=ϵ
            end
        end
        ∂1 = ∂exp(QuatVecF64(0))
        ∂2 = QuaternionF64[1, imx, imy, imz]
        for i in 1:4
            @test ∂1[i] ≈ ∂2[i] rtol=ϵ
        end
        e, ∂1 = exp∂exp(QuatVecF64(0))
        @test e ≈ QuaternionF64(1) rtol=4eps() atol=10eps()
        for i in 1:4
            @test ∂1[i] ≈ ∂2[i] rtol=ϵ
        end
    end

    @testset verbose=true "slerp" begin
        ∂slerp(q₁, q₂, τ) = (
            [
                ForwardDiff.derivative(ϵ->slerp(q₁+(ϵ+0imx), q₂, τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁+ϵ*imx, q₂, τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁+ϵ*imy, q₂, τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁+ϵ*imz, q₂, τ), 0)
            ],
            [
                ForwardDiff.derivative(ϵ->slerp(q₁, q₂+(ϵ+0imx), τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁, q₂+ϵ*imx, τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁, q₂+ϵ*imy, τ), 0),
                ForwardDiff.derivative(ϵ->slerp(q₁, q₂+ϵ*imz, τ), 0)
            ]
        )
        for (q₁, q₂) in zip(randn(RotorF64, 1_000), randn(RotorF64, 1_000))
            for τ ∈ [0.0; 1.0; rand(5)...]
                ∂1, ∂2 = ∂slerp(q₁, q₂, τ)
                sa, ∂a, ∂b, ∂c = slerp∂slerp(q₁, q₂, τ)
                @test sa ≈ slerp(q₁, q₂, τ) rtol=2ϵ
                for i in 1:4
                    @test ∂1[i] ≈ ∂a[i] rtol=5ϵ atol=50eps()
                    @test ∂2[i] ≈ ∂b[i] rtol=5ϵ atol=50eps()
                end
                @test ∂c ≈ log(q₂/q₁) * slerp(q₁, q₂, τ) rtol=2ϵ
                sb, ∂d = slerp∂slerp∂τ(q₁, q₂, τ)
                @test sa ≈ sb rtol=2ϵ
                @test ∂c ≈ ∂d rtol=2ϵ
            end
        end
    end

    @testset verbose=true "squad" begin
        qs = Rotor{Float64}[1, imx, imy, imz, -imy, -imz, -imx, -rotor(1)]
        ts = Float64.(1:length(qs))
        ∂squad(t) = ForwardDiff.derivative(τ->squad(qs, ts, τ), t)
        for i ∈ 1:length(ts)-1
            qᵢ, qᵢ₊₁ = qs[i], qs[i+1]
            ta, tb = ts[i], ts[i+1]
            for τ ∈ [0.001, 0.1, 0.49, 0.5, 0.51, 0.9, 0.999]
                t = ta + τ*(tb-ta)
                A, B = Quaternionic.squad_control_points(qs, ts, i)
                ∂1 = ∂squad(t)
                s, ∂2 = squad∂squad∂t(qᵢ, A, B, qᵢ₊₁, ta, tb, t)
                @test ∂1 ≈ ∂2 rtol=100ϵ atol=ϵ
                @test s ≈ squad(qs, ts, t) rtol=ϵ atol=ϵ
            end
        end
    end
end
