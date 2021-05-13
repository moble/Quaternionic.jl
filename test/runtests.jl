using Quaternionic
using Symbolics
using Test

@variables w x y z a b c d

Types = [
    Float64, Float32, Float16,
    Int128, Int64, Int32, Int16, Int8,
    BigFloat, BigInt,
    Num
]

@testset verbose=true "Basis" begin
    @testset "$T" for T in Types
        # Typed constants
        ϵ = T<:AbstractFloat ? eps(T) : zero(T)
        unit = T(1)
        nil = T(0)

        # Define basis elements
        u = Quaternion(unit, nil, nil, nil)
        i = Quaternion(nil, unit, nil, nil)
        j = Quaternion(nil, nil, unit, nil)
        k = Quaternion(nil, nil, nil, unit)
        basis = [u, i, j, k]

        # Check basis elements
        for (index, element) in enumerate(basis)
            for component in 1:4
                if component == index
                    @test element[component] == unit
                else
                    @test element[component] == nil
                end
            end
        end

        # Because of signed zeros, we can't check for exact equality, so we check for approximate
        # equality with atol=0 and rtol either zero (for exact types) or eps (for floats).

        # Standard expressions
        @test u * u ≈ u rtol=ϵ
        @test i * i ≈ -u rtol=ϵ
        @test j * j ≈ -u rtol=ϵ
        @test k * k ≈ -u rtol=ϵ
        @test i * j * k ≈ -u rtol=ϵ

        # Full multiplication table
        @test u * u ≈ u rtol=ϵ
        @test u * i ≈ i rtol=ϵ
        @test u * j ≈ j rtol=ϵ
        @test u * k ≈ k rtol=ϵ
        @test i * u ≈ i rtol=ϵ
        @test i * i ≈ -u rtol=ϵ
        @test i * j ≈ k rtol=ϵ
        @test i * k ≈ -j rtol=ϵ
        @test j * u ≈ j rtol=ϵ
        @test j * i ≈ -k rtol=ϵ
        @test j * j ≈ -u rtol=ϵ
        @test j * k ≈ i rtol=ϵ
        @test k * u ≈ k rtol=ϵ
        @test k * i ≈ j rtol=ϵ
        @test k * j ≈ -i rtol=ϵ
        @test k * k ≈ -u rtol=ϵ
    end
end
