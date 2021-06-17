@testset verbose=true "Math" begin
    components = [:x, :y, :z]

    q_to_c(q::Real, component) = q + 0im
    function q_to_c(q::Quaternion, component)
        c = q.w + im * getproperty(q, component)
    end

    function c_to_q(c::Complex, component)
        Quaternion(c.re, [comp_i == component ? c.im : zero(real(c)) for comp_i in components]...)
    end

    basic_functions = [conj, abs2, abs, inv, log, exp, sqrt]

    @testset "$T" for T in FloatTypes
        scalars = [zero(T), one(T), -one(T)]
        for c in [a+b*im for a in scalars for b in scalars]
            for component in components
                q = c_to_q(c, component)
                for unary_function in basic_functions
                    if unary_function ∈ [log, sqrt] && component!=:z && c.re < zero(T) && c.im == zero(T)
                        continue  # We arbitrarily chose the z component for these return values
                    end
                    # println((unary_function, T, c, q))
                    cval = unary_function(c)
                    qval = q_to_c(unary_function(q), component)
                    @test cval ≈ qval rtol=eps(T) nans=true
                end
                for unary_function in [angle]  # Not identical, but related
                    cval = unary_function(c)
                    qval = q_to_c(unary_function(q), component)
                    @test 2*abs(cval) ≈ qval rtol=eps(T) nans=true
                end
                @test exp(log(q)) ≈ q rtol=eps(T)
                @test exp(zero(q)) ≈ one(q) rtol=eps(T)
            end
        end
    end

end
