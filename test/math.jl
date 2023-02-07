@testset verbose=true "Math" begin
    components = [:x, :y, :z]

    basic_functions = [conj, abs2, abs, inv, log, exp, sqrt]

    ### Test equivalence with complex math
    ℍ_to_ℂ(q::Real, component) = q + 0im
    function ℍ_to_ℂ(q::AbstractQuaternion, component)
        c = q[1] + im * getproperty(q, component)
    end
    function ℂ_to_ℍ(c::Complex, component)
        Quaternion(real(c), [comp_i == component ? imag(c) : zero(real(c)) for comp_i in components]...)
    end
    @testset "Complex equivalence $T" for T in FloatTypes
        ϵ = (T === Float16 ? 20eps(T) : 10eps(T))
        scalars = [zero(T), one(T), -one(T)]
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
                    @test cval ≈ qval rtol=ϵ nans=true
                    # Repeat the test on QuatVec for `exp` for pure-imaginary input
                    if unary_function ∈ [exp] && iszero(real(c))
                        qval = ℍ_to_ℂ(unary_function(QuatVec(q)), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                    end
                    # Repeat the test on Rotor for `log` and `sqrt` for nonzero input
                    if unary_function ∈ [log, sqrt] && !iszero(c)
                        cval = unary_function(c/abs(c))
                        qval = ℍ_to_ℂ(unary_function(Rotor(q)), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                    end
                end
                for unary_function in [angle]  # Not identical, but related
                    cval = unary_function(c)
                    qval = ℍ_to_ℂ(unary_function(q), component)
                    @test 2*abs(cval) ≈ qval rtol=ϵ nans=true
                    if !iszero(c)
                        cval = unary_function(c/abs(c))
                        qval = ℍ_to_ℂ(unary_function(Rotor(q)), component)
                        @test 2*abs(cval) ≈ qval rtol=ϵ nans=true
                    end
                end
                for binary_function in [(^)]
                    if binary_function ∈ [(^)] && ((component!=:z && real(c) < zero(T) && imag(c) == zero(T)) || iszero(c))
                        continue  # We arbitrarily chose the z component for these return values
                    end
                    if T !== BigFloat
                        S = 10randn(T, 100)
                        for s ∈ S
                            cval = binary_function(c, float(s))
                            qval = ℍ_to_ℂ(binary_function(q, s), component)
                            let ϵ = (T === Float16 ? 10ϵ/6 : ϵ)
                                @test cval ≈ qval rtol=6ϵ nans=true
                            end
                            if !iszero(c)
                                cval = binary_function(c/abs(c), s)
                                qval = ℍ_to_ℂ(binary_function(Rotor(q), s), component)
                                @test cval ≈ qval rtol=6ϵ nans=true
                            end
                        end
                    end
                    S = [-3, -2, -1, 0, 1, 2, 3]
                    for s ∈ S
                        cval = binary_function(c, s)
                        qval = ℍ_to_ℂ(binary_function(q, s), component)
                        @test cval ≈ qval rtol=ϵ nans=true
                        if iszero(real(c))
                            qval = ℍ_to_ℂ(binary_function(QuatVec(q), s), component)
                            @test cval ≈ qval rtol=ϵ nans=true
                        end
                        if !iszero(c)
                            cval = binary_function(c/abs(c), s)
                            qval = ℍ_to_ℂ(binary_function(Rotor(q), s), component)
                            @test cval ≈ qval rtol=ϵ nans=true
                        end
                    end
                end

                @test exp(log(q)) ≈ q rtol=ϵ nans=true
                @test exp(zero(q)) ≈ one(q) rtol=ϵ nans=true
                @test exp(log(Rotor(q))) ≈ Rotor(q) rtol=ϵ nans=true
                @test exp(zero(QuatVec(q))) ≈ one(q) rtol=ϵ nans=true

            end
        end
    end

end
