# The quaternions form an associative normed division algebra

module FundamentalTests
    using Test: @test
    using Quaternionic

    # Algebra
    ## Vector space
    test_vector_associativity(u::Quaternion, v::Quaternion, w::Quaternion) = @test u + (v + w) ≈ (u + v) + w rtol=eps(v)
    test_vector_commutativity(u::Quaternion, v::Quaternion) = @test u + v ≈ v + u rtol=eps(v)
    function test_vector_identity(v::Quaternion)
        @test v + zero(v) == v
        @test zero(v) + v == v
    end
    test_vector_inverse(v::Quaternion) = @test v + (-v) == zero(v)
    test_vector_scalar_identity(v::Quaternion) = @test one(eltype(v)) * v == v
    test_vector_scalar_associativity(a, b, v::Quaternion) = @test a*(b*v) ≈ (a*b)*v rtol=eps(v)
    test_vector_scalar_distributivity1(a, u::Quaternion, v::Quaternion) = @test a*(u+v) ≈ (a*u) + (a*v) rtol=eps(v)
    test_vector_scalar_distributivity2(a, b, v::Quaternion) = @test (a+b)*v ≈ (a*v) + (b*v) rtol=eps(v)
    ## Product
    test_right_distributivity(x::Quaternion, y::Quaternion, z::Quaternion) = @test (x+y)*z ≈ (x*z) + (y*z) rtol=eps(x)
    test_left_distributivity(x::Quaternion, y::Quaternion, z::Quaternion) = @test z*(x+y) ≈ (z*x) + (z*y) rtol=eps(x)
    test_scalar_compatibility(a, b, x::Quaternion, y::Quaternion) = @test (a*x) * (b*y) ≈ (a*b) * (x*y) rtol=eps(x)

    # Division
    function test_identity(v::Quaternion)
        @test one(v) * v == v
        @test v * one(v) == v
    end
    function test_inverse(v::Quaternion)
        if !iszero(v)
            @test v * inv(v) == one(v)
            @test inv(v) * v == one(v)
        end
    end

    # Normed
    test_norm_maps_to_field(q::Quaternion) = @test typeof(abs2(q)) == eltype(q)
    test_norm_nondegenerate(v::Quaternion) = @test iszero(v) ⊻ !iszero(abs2(v))
    test_norm_quadratic(v::Quaternion) = @test abs2(2*v) ≈ 4*abs2(v) rtol=eps(v)
    test_norm_quadratic(a, v::Quaternion) = @test abs2(a*v) ≈ a^2*abs2(v) rtol=eps(v)

    # Associative
    test_associativity(u::Quaternion, v::Quaternion, w::Quaternion) = @test (u * v) * w ≈ u * (v * w) rtol=eps(v)
end


@testset verbose=true "Fundamentals" begin
    @testset "$T" for T in [FloatTypes...; IntTypes...]
        println("Testing fundamentals with type $T")

        # Construct a variety of arguments
        scalars = [zero(T), one(T), -one(T)]#, 2*one(T), -2*one(T)]#, eps(T), -eps(T)]
        quaternions = [Quaternion(a, b, c, d) for a in scalars for b in scalars for c in scalars for d in scalars]

        # Iterate over all tests above
        for n in names(FundamentalTests, all=true)
            if T <: Integer && (n == :test_inverse)
                continue  # Don't try to invert integer quaternions
            end
            f = getproperty(FundamentalTests, n)
            if !isempty(methods(f)) && startswith(string(f), "test_")
                types = methods(f).ms[1].sig.parameters[2:end]
                args = Base.Iterators.product([type===Any ? scalars : quaternions for type in types]...)
                for arg in args
                    f(arg...)
                end
            end
        end
    end
    @testset "$T" for T in SymbolicTypes
        println("Testing fundamentals with type $T")
        chars = Iterators.Stateful(Iterators.cycle("abcdefghijkl"))
        function next_scalar!(chars)
            x = Symbol(popfirst!(chars))
            xvar = @variables $x
            xvar[1]
        end
        function next_quaternion!(chars)
            x = popfirst!(chars)
            xvar = @variables $x[1:4]
            Quaternion(xvar...)
        end

        # Iterate over all tests above
        for n in names(FundamentalTests, all=true)
            f = getproperty(FundamentalTests, n)
            if !isempty(methods(f)) && startswith(string(f), "test_")
                types = methods(f).ms[1].sig.parameters[2:end]
                args = [type===Any ? next_scalar!(chars) : next_quaternion!(chars) for type in types]
                f(args...)
            end
        end
    end
end
