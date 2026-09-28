# Arithmetic under `@fastmath` must give the same types and values as ordinary arithmetic.
# See issue #46: the generic fallbacks in `Base.FastMath` promote mixed arguments before
# applying the operator, so that `@fastmath 1.0 * imz` returned a `Quaternion` instead of a
# `QuatVec`, and `@fastmath 2.0 * r` returned the `Rotor` `r` instead of `2r`.

@testmodule FastMathCases begin
    using Quaternionic

    """
        matches(fast, slow)

    Return `true` if the result `fast` of an expression under `@fastmath` has the same type
    as the result `slow` of the ordinary expression, and is identical to it: `isequal`
    component by component, so that even the signs of zeros must agree.  (Comparing with
    `==` alone would not suffice, because a `Quaternion` can equal a `QuatVec`.)  If the
    result is a `Rotor`, it must also have unit norm.  That norm is computed from the
    components, because `abs` and `abs2` simply return 1 for any `Rotor`.
    """
    matches(fast, slow) = typeof(fast) === typeof(slow) && isequal(fast, slow)
    function matches(fast::AbstractQuaternion, slow::AbstractQuaternion)
        typeof(fast) === typeof(slow) && isequal(components(fast), components(slow)) &&
            (!(fast isa Rotor) || unit_norm(fast))
    end

    """
        unit_norm(R)

    Return `true` if the components of `R` have a sum of squares equal to 1, up to roundoff.

    The tolerance is set by `Float32` even for a `RotorF64`, because a `RotorF32` raised to
    a `Float64` power is a `RotorF64` that is only as accurate as its `Float32` input.
    """
    unit_norm(R) = isapprox(sum(abs2, components(R)), 1; rtol=10eps(Float32))

    const scalars = (2.0, 3, 1.5f0)
    const quaternions = (
        quaternion(1.0, 2.0, 3.0, 4.0),
        rotor(1.0, 2.0, 3.0, 4.0),
        quatvec(1.0, 2.0, 3.0),
        quaternion(1.0f0, -2.0f0, 3.0f0, 0.5f0),
        rotor(-1.0f0, 2.0f0, 0.5f0, 4.0f0),
        quatvec(-1.0f0, 2.0f0, 0.5f0),
    )
    const values = (scalars..., quaternions...)
end


@testitem "fastmath: binary arithmetic" tags=[:unit, :fast] setup=[FastMathCases] begin
    using .FastMathCases: matches, values
    for a ∈ values, b ∈ values
        (a isa AbstractQuaternion || b isa AbstractQuaternion) || continue
        @test matches(@fastmath(a + b), a + b)
        @test matches(@fastmath(a - b), a - b)
        @test matches(@fastmath(a * b), a * b)
        @test matches(@fastmath(a / b), a / b)
    end
end


@testitem "fastmath: issue #46" tags=[:unit, :fast] setup=[FastMathCases] begin
    using .FastMathCases: unit_norm
    @test @fastmath(1.0 * imz) isa QuatVecF64
    @test @fastmath(imz * 1.0) isa QuatVecF64
    @test @fastmath(imz / 2.0) isa QuatVecF64

    # Scaling a rotor gives a general quaternion with the scaled norm...
    r₁ = rotor(1.0, 2.0, 3.0, 4.0)
    r₂ = rotor(-1.0, 0.5, 2.0, 1.0)
    @test @fastmath(2.0 * r₁) isa QuaternionF64
    @test components(@fastmath(2.0 * r₁)) == 2 * components(r₁)
    @test components(@fastmath(r₁ / 2.0)) == components(r₁) / 2

    # ...but products, quotients, and powers of rotors are rotors with unit norm.
    for R ∈ @fastmath((r₁ * r₂, r₁ / r₂, r₁ * r₂ * r₁, r₁ ^ 2.5, r₁ ^ 3, r₁ ^ -1, -r₁))
        @test R isa RotorF64
        @test unit_norm(R)
    end
end


@testitem "fastmath: chained arithmetic" tags=[:unit, :fast] setup=[FastMathCases] begin
    using .FastMathCases: matches, scalars, quaternions
    # `@fastmath a * b * c` calls `mul_fast(a, b, c)`, which is handled separately from the
    # binary case.  Put a quaternion in each of the first three positions.
    for q ∈ quaternions, s ∈ scalars, t ∈ scalars
        for (a, b, c) ∈ ((q, s, t), (s, q, t), (s, t, q), (q, q, s), (q, s, q), (s, q, q), (q, q, q))
            @test matches(@fastmath(a * b * c), a * b * c)
            @test matches(@fastmath(a + b + c), a + b + c)
            @test matches(@fastmath(a * b * c * s), a * b * c * s)
            @test matches(@fastmath(a + b + c + s), a + b + c + s)
        end
    end
end


@testitem "fastmath: powers" tags=[:unit, :fast] setup=[FastMathCases] begin
    using .FastMathCases: matches, quaternions
    for q ∈ quaternions
        for s ∈ (2.5, 0.5f0, 3, -2)
            @test matches(@fastmath(q ^ s), q ^ s)
        end
        @test matches(@fastmath(q ^ 2), q ^ 2)
        @test matches(@fastmath(q ^ -1), q ^ -1)
    end
end


@testitem "fastmath: unary functions" tags=[:unit, :fast] setup=[FastMathCases] begin
    using .FastMathCases: matches, quaternions
    for q ∈ quaternions
        @test matches(@fastmath(-q), -q)
        @test matches(@fastmath(conj(q)), conj(q))
        @test matches(@fastmath(inv(q)), inv(q))
        @test matches(@fastmath(abs(q)), abs(q))
        @test matches(@fastmath(abs2(q)), abs2(q))
        @test matches(@fastmath(sqrt(q)), sqrt(q))
        # `exp` is not defined for a `Rotor`, nor `log` for a `QuatVec`.
        if !(q isa Rotor)
            @test matches(@fastmath(exp(q)), exp(q))
        end
        if !(q isa QuatVec)
            @test matches(@fastmath(log(q)), log(q))
        end
    end
end
