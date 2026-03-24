@testset verbose=true "Math" begin
    components = [:x, :y, :z]

    basic_functions = [conj, abs2, abs, inv, log, exp, sqrt]

    ### Test equivalence with complex math
    ℍ_to_ℂ(q::Real, component) = q + 0im
    function ℍ_to_ℂ(q::AbstractQuaternion, component)
        c = q[1] + im * getproperty(q, component)
    end
    function ℂ_to_ℍ(c::Complex, component)
        quaternion(real(c), [comp_i == component ? imag(c) : zero(real(c)) for comp_i in components]...)
    end
    @testset "Complex equivalence $T" for T in FloatTypes
        ϵ = (T === Float16 ? 20eps(T) : 10eps(T))
        scalars = [zero(T), one(T), -one(T), one(T)*7//3, -one(T)*5//2]
        for c in [a+b*im for a in scalars for b in scalars]
            # The following constructs a quaternion from the complex number `c` by equating `im`
            # with one of the three quaternionic imaginaries, looping over each option.  We then
            # apply the quaternionic math function to that quaternion, and map back to ℂ; the result
            # should be identical to just applying the complex math function directly.  The only
            # exceptions are for `log` and `sqrt` when the complex number is a negative real number,
            # for which we have to arbitrarily choose the "vector" part of the quaternion, which is
            # chosen to be `imz` in the code.
            for component in components
                q = ℂ_to_ℍ(c, component)
                for unary_function in basic_functions
                    if unary_function ∈ [log, sqrt] && component!=:z && real(c) < zero(T) && imag(c) == zero(T)
                        continue  # We arbitrarily chose the z component for these return values
                    end
                    cval = unary_function(c)
                    qval = ℍ_to_ℂ(unary_function(q), component)
                    if !≈(cval, qval, rtol=ϵ, nans=true)
                        @info "A" unary_function c q cval qval unary_function(q)
                    end
                    @test cval ≈ qval rtol=ϵ nans=true
                    # Repeat the test on QuatVec for `exp` for pure-imaginary input
                    if unary_function ∈ [exp] && iszero(real(c))
                        qval = ℍ_to_ℂ(unary_function(quatvec(q)), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                    end
                    # Repeat the test on Rotor for `log` and `sqrt` for nonzero input
                    if unary_function ∈ [log, sqrt] && !iszero(c)
                        cval = unary_function(c/abs(c))
                        qval = ℍ_to_ℂ(unary_function(rotor(q)), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                    end
                end
                for unary_function in [angle]  # Not identical, but related
                    cval = unary_function(c)
                    qval = ℍ_to_ℂ(unary_function(q), component)
                    @test 2*abs(cval) ≈ qval rtol=ϵ nans=true
                    if !iszero(c)
                        cval = unary_function(c/abs(c))
                        qval = ℍ_to_ℂ(unary_function(rotor(q)), component)
                        @test 2*abs(cval) ≈ qval rtol=ϵ nans=true
                    end
                end
                for binary_function in [(^)]
                    if binary_function ∈ [(^)] && ((component!=:z && real(c) < zero(T) && imag(c) == zero(T)) || iszero(c))
                        continue  # We arbitrarily chose the z component for these return values
                    end
                    for s ∈ (T === Float16 ? 1 : 10) * randn(T, 50)
                        cval = binary_function(c, float(s))
                        qval = ℍ_to_ℂ(binary_function(q, s), component)
                        let ϵ = (T === BigFloat ? 15ϵ : 10ϵ)
                            @test cval ≈ qval rtol=ϵ nans=true
                        end
                        if !iszero(c)
                            cval = binary_function(c/abs(c), s)
                            qval = ℍ_to_ℂ(binary_function(rotor(q), s), component)
                            @test cval ≈ qval rtol=6ϵ nans=true
                        end
                    end
                    for s ∈ (-3, -2, -1, 0, 1, 2, 3)
                        cval = binary_function(c, s)
                        qval = ℍ_to_ℂ(binary_function(q, s), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                        if iszero(real(c))
                            qval = ℍ_to_ℂ(binary_function(quatvec(q), s), component)
                            @test cval ≈ qval rtol=ϵ nans=true
                        end
                        if !iszero(c)
                            cval = binary_function(c/abs(c), s)
                            qval = ℍ_to_ℂ(binary_function(rotor(q), s), component)
                            @test cval ≈ qval rtol=ϵ nans=true
                        end
                    end
                end

                @test exp(log(q)) ≈ q rtol=ϵ nans=true
                @test exp(zero(q)) ≈ one(q) rtol=ϵ nans=true
                @test exp(log(rotor(q))) ≈ rotor(q) rtol=ϵ nans=true
                @test exp(zero(quatvec(q))) ≈ one(q) rtol=ϵ nans=true

                @test exp(2q) ≈ exp(q)^2 rtol=ϵ nans=true
            end
        end
    end

    @testset "Special values for abs $T" for T in FloatTypes
        # abs2 and abs2vec now return the spinor norm (Σzᵢ²) for Complex{T} components,
        # not the Euclidean norm (Σ|zᵢ|²).  Rotor always returns one(real(T)) regardless.
        @test abs2(Quaternion{Complex{T}}(1+2im, 3+4im, false, false)) == Complex{T}(-10, 28)
        @test abs2(QuatVec{Complex{T}}(1+2im, 3+4im, false, false)) == Complex{T}(-7, 24)
        @test abs2(Rotor{Complex{T}}(1+2im, 3+4im, false, false)) == one(T)
        @test abs2vec(Quaternion{Complex{T}}(1+2im, 3+4im, false, false)) == Complex{T}(-7, 24)
        @test abs2vec(QuatVec{Complex{T}}(1+2im, 3+4im, false, false)) == Complex{T}(-7, 24)
        @test abs2vec(Rotor{Complex{T}}(1+2im, 3+4im, false, false)) == Complex{T}(-7, 24)
    end

    @testset "Special values for sqrt $T" for T in FloatTypes
        ϵ = 4eps(T)

        # sqrt(0) = 0
        q = quaternion(zero(T), zero(T), zero(T), zero(T))
        @test q ≈ sqrt(quaternion(zero(T), zero(T), zero(T), zero(T))) rtol=ϵ nans=true
        @test sqrt(q) ≈ quaternion(zero(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true

        # sqrt(1) = 1
        q = quaternion(one(T), zero(T), zero(T), zero(T))
        @test q ≈ sqrt(quaternion(one(T), zero(T), zero(T), zero(T))) rtol=ϵ nans=true
        @test sqrt(q) ≈ quaternion(one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true
        q = rotor(one(T), zero(T), zero(T), zero(T))
        @test q ≈ sqrt(rotor(one(T), zero(T), zero(T), zero(T))) rtol=ϵ nans=true
        @test sqrt(q) ≈ rotor(one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true

        # sqrt(-1) = 𝐤
        q = quaternion(-one(T), zero(T), zero(T), zero(T))
        @test sqrt(q) ≈ quaternion(zero(T), zero(T), zero(T), one(T)) rtol=ϵ nans=true
        @test sqrt(q)^2 ≈ quaternion(-one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true
        q = rotor(-one(T), zero(T), zero(T), zero(T))
        @test sqrt(q) ≈ rotor(zero(T), zero(T), zero(T), one(T)) rtol=ϵ nans=true
        @test sqrt(q)^2 ≈ rotor(-one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true

        # sqrt(4) = 2
        q = quaternion(4one(T), zero(T), zero(T), zero(T))
        @test sqrt(q) ≈ quaternion(2one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true

        # sqrt(-4) = 2𝐤
        q = quaternion(-4one(T), zero(T), zero(T), zero(T))
        @test sqrt(q) ≈ quaternion(zero(T), zero(T), zero(T), 2one(T)) rtol=ϵ nans=true
        @test sqrt(q)^2 ≈ quaternion(-4one(T), zero(T), zero(T), zero(T)) rtol=ϵ nans=true

        # sqrt(-1 ± ε𝐯) = ε/2 ± 𝐯
        for (s,ε,𝐯) ∈ ((s,ε,𝐯)
            for s ∈ (-1,1)
            for ε ∈ (4√eps(T), √eps(T)/4, 4eps(T), eps(T), eps(T)^T(7//3))
            for 𝐯 ∈ (𝐢, 𝐣, 𝐤)
        )
            if ≉(sqrt(-1 + s*ε*𝐯), ε/2 + s*𝐯, rtol=ϵ, nans=true)
                @info "sqrt test failure A" s ε 𝐯 sqrt(-1 + s*ε*𝐯) ε/2 + s*𝐯
            end
            @test sqrt(-1 + s*ε*𝐯) ≈ ε/2 + s*𝐯 rtol=ϵ nans=true
            if ≉(sqrt(-1 + s*ε*𝐯)^2, -1 + s*ε*𝐯, rtol=ϵ, nans=true)
                @info "sqrt test failure B" s ε 𝐯 sqrt(-1 + s*ε*𝐯)^2 -1 + s*ε*𝐯
            end
            @test sqrt(-1 + s*ε*𝐯)^2 ≈ -1 + s*ε*𝐯 rtol=ϵ nans=true
        end

    end

    @testset "Special values for log $T" for T in FloatTypes
        ϵ = (T === Float16 ? 20eps(T) : 10eps(T))

        # log(0) = -Inf
        q = quaternion(zero(T), zero(T), zero(T), zero(T))
        @test q ≈ exp(quaternion(-T(Inf), 0, 0, 0)) rtol=ϵ nans=true
        @test log(q) ≈ quaternion(-T(Inf), 0, 0, 0) rtol=ϵ nans=true

        # log(1) = 0
        q = quaternion(one(T), zero(T), zero(T), zero(T))
        @test q ≈ exp(quaternion(zero(T), 0, 0, 0)) rtol=ϵ nans=true
        @test log(q) ≈ quaternion(zero(T), 0, 0, 0) rtol=ϵ nans=true
        q = rotor(one(T), zero(T), zero(T), zero(T))
        @test q ≈ exp(quatvec(0, zero(T), 0, 0)) rtol=ϵ nans=true
        @test log(q) ≈ quatvec(0, zero(T), 0, 0) rtol=ϵ nans=true

        # log(-1) = π𝐤
        q = quaternion(-one(T), zero(T), zero(T), zero(T))
        @test q ≈ exp(quaternion(0, 0, 0, T(π))) rtol=ϵ nans=true
        @test log(q) ≈ quaternion(0, 0, 0, T(π)) rtol=ϵ nans=true
        q = rotor(-one(T), zero(T), zero(T), zero(T))
        @test q ≈ exp(quatvec(0, 0, 0, T(π))) rtol=ϵ nans=true
        @test log(q) ≈ quatvec(0, 0, 0, T(π)) rtol=ϵ nans=true

        for (s,f,v) ∈ ((s,f,v)
                for s ∈ (-1,1)
                for f ∈ (4√eps(T), √eps(T)/4, 4eps(T), eps(T))
                for v ∈ (𝐢, 𝐣, 𝐤)
            )
            Δ = s*f*v

            # log(1) + Δ = Δ
            q = quaternion(one(T)) + Δ
            @test q ≈ exp(quaternion(log(abs(q))) + Δ) rtol=ϵ nans=true
            @test log(q) ≈ quaternion(log(abs(q))) + Δ rtol=ϵ nans=true
            r = rotor(q)
            @test r ≈ exp(Δ) rtol=ϵ nans=true
            @test log(r) ≈ Δ rtol=ϵ nans=true

            # log(-1) + Δ = (π - |Δ|) Δ/|Δ|
            q = quaternion(-one(T)) + Δ
            @test q ≈ exp(quaternion(log(abs(q))) + (T(π)-f)*s*v) rtol=ϵ nans=true
            @test log(q) ≈ quaternion(log(abs(q))) + (T(π)-f)*s*v rtol=ϵ nans=true
            r = rotor(q)
            @test r ≈ exp((T(π)-f)*s*v) rtol=ϵ nans=true
            @test log(r) ≈ (T(π)-f)*s*v rtol=ϵ nans=true

            q = quaternion(one(T) * 3//2) + Δ
            @test q ≈ exp(quaternion(log(abs(q))) + Δ * 2//3) rtol=ϵ nans=true
            @test log(q) ≈ quaternion(log(abs(q))) + Δ * 2//3 rtol=ϵ nans=true

        end
    end

end
