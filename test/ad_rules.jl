# Tests of the `frule`s, `rrule`s, and projections of `QuaternionicChainRulesCoreExt`.
#
# Each rule is checked in two ways.  `checkrule` (in the `RuleChecks` snippet below) compares
# the pullback and the pushforward with central finite differences evaluated in 256-bit
# `BigFloat` arithmetic (`bigvjp` and `bigjvp` from `ADTestUtils`), at every point of the
# point set where the source is differentiable; this includes the points near the identity
# and near -1, where the Float64 finite differences of ChainRulesTestUtils are too coarse.
# `crtu` runs ChainRulesTestUtils' `test_rrule` and `test_frule` at generic points, which
# also check inference, thunks, zero cotangents, and the types of the tangents.  Points
# where the source is not differentiable in every direction (the negative real axis, and
# the cone points of norms) are checked explicitly.
#
# ChainRules itself is loaded in every item, because the purpose of many of these rules is
# to displace the generic `Number` rules of ChainRules, which would otherwise capture
# quaternion arguments.

@testsnippet RuleChecks begin
    using Quaternionic
    using Test
    using ChainRulesCore
    using ChainRulesCore: rrule, frule, unthunk
    using ChainRulesTestUtils
    import ChainRules
    using Random: Random, MersenneTwister
    using LinearAlgebra: LinearAlgebra, dot, norm, normalize
    using StaticArrays: SVector, @SVector
    using .ADTestUtils

    rng = MersenneTwister(20261002)

    # Reseed both sources of random tangents: `rng`, from which this module draws the
    # tangents it chooses, and the global generator, from which ChainRulesTestUtils draws
    # the rest.  `checkrule` and `crtu` call this first, so that each check uses the same
    # tangents however the test items are filtered or ordered, and a failure can be
    # reproduced by running its item alone.
    function reseed!()
        Random.seed!(rng, 20261002)
        Random.seed!(20261002)
        nothing
    end

    # A random tangent or cotangent of `x`.  Arguments that are not numbers or arrays of
    # numbers (functions, `Val`s, and characters) get `NoTangent()`, as do integers and
    # `Bool`s.
    randtangent(x) = NoTangent()
    randtangent(x::Union{Number,AbstractArray{<:Number}}) = rand_tangent(rng, x)
    randtangent(x::Rational) = randn(rng, float(typeof(x)))

    # Whether `∂` has the type of a tangent of `x` in the conventions of the extension: a
    # `Quaternion` for a `Quaternion` or `Rotor`, a `QuatVec` for a `QuatVec`, each with
    # floating-point components of the precision of `x`; a number of the same kind for a
    # number; and an array of the same kind for an array.
    tangentok(∂, x) = true
    tangentok(∂, x::AbstractQuaternion) = ∂ isa Quaternion{float(typeof(x[1]))}
    tangentok(∂, x::QuatVec) = ∂ isa QuatVec{float(typeof(x[1]))}
    tangentok(∂, x::Real) = ∂ isa Real
    tangentok(∂, x::Complex) = ∂ isa Number && !(∂ isa AbstractQuaternion)
    tangentok(∂, x::AbstractArray) = ∂ isa AbstractArray && size(∂) == size(x)
    tangentok(∂, x::Vector) = ∂ isa Vector && size(∂) == size(x)

    # The discrepancy between a tangent `∂` returned by a rule and the reference `ref` for
    # the primal `x`: the relative error, or `Inf` if the type is wrong.  A reference of
    # `NoTangent()` requires a zero tangent, and a zero tangent counts as the zero of the
    # reference.
    function discrepancy(∂, ref, x)
        ∂ = unthunk(∂)
        ref isa AbstractZero && return ∂ isa AbstractZero ? 0.0 : Inf
        ∂ isa AbstractZero && return Float64(relerr(zero(ref), ref))
        tangentok(∂, x) || return Inf
        return Float64(relerr(∂, ref))
    end

    describe(f) = string(f)
    describe(f::Rotor) = "R(v)"
    describe(f::Type) = string(f)

    """
        checkrule(f, args...; ȳ, ẋs, ref, rtol, pullback, pushforward)

    Check `rrule(f, args...)` and `frule((NoTangent(), ẋs...), f, args...)` against the
    256-bit finite-difference references `bigvjp(ref, args, ȳ)` and `bigjvp(ref, args, ẋs)`,
    to the relative error `rtol`.  The primal must equal `f(args...)` in value and type, and
    each tangent must have the type of the conventions.  By default, `ȳ` and `ẋs` are random,
    and `ref` is `f` itself.  When `f` is a `Rotor`, it is the rotor of `R(v)`, and it is
    differentiated as an argument.
    """
    function checkrule(f, args...; ȳ=nothing, ẋs=nothing, ref=f, rtol=1e-10,
                       pullback=true, pushforward=true)
        reseed!()
        selfdiff = f isa Rotor
        g = selfdiff ? ((h, a...) -> h(a...)) : ref
        gargs = selfdiff ? (f, args...) : args
        Ω₀ = f(args...)
        @testset "$(describe(f))($(join(map(a -> string(typeof(a)), args), ", ")))" begin
            if pullback
                res = rrule(f, args...)
                @test res !== nothing
                if res !== nothing
                    Ω, pb = res
                    @test typeof(Ω) == typeof(Ω₀)
                    @test Ω ≈ Ω₀
                    Δ = ȳ === nothing ? randtangent(Ω₀) : ȳ
                    ∂ = pb(Δ)
                    @test length(∂) == length(gargs) + !selfdiff
                    selfdiff || @test ∂[1] isa NoTangent
                    refs = bigvjp(g, gargs, Δ)
                    ∂args = selfdiff ? ∂ : Base.tail(∂)
                    for i ∈ eachindex(gargs)
                        a = gargs[i]
                        if a isa Integer && !(a isa Bool) && unthunk(∂args[i]) isa Real
                            # The reference cannot perturb an integer, so a real cotangent of
                            # an integer (which `ProjectTo` of an integer gives) is compared
                            # with the reference for the equivalent float.
                            floatrefs = bigvjp(g, Base.setindex(gargs, float(a), i), Δ)
                            @test discrepancy(∂args[i], floatrefs[i], float(a)) ≤ rtol
                        else
                            @test discrepancy(∂args[i], refs[i], a) ≤ rtol
                        end
                    end
                end
            end
            if pushforward
                ẋ = ẋs === nothing ? map(randtangent, gargs) : ẋs
                res = selfdiff ? frule(ẋ, f, args...) : frule((NoTangent(), ẋ...), f, args...)
                @test res !== nothing
                if res !== nothing
                    Ω, Ω̇ = res
                    @test typeof(Ω) == typeof(Ω₀)
                    @test Ω ≈ Ω₀
                    @test discrepancy(Ω̇, bigjvp(g, gargs, ẋ), Ω₀) ≤ rtol
                end
            end
        end
    end

    # Tolerances of ChainRulesTestUtils' finite differences, which are computed in the
    # precision of the primal
    crtutol(::Type{Float32}) = (rtol=5e-3, atol=5e-3)
    crtutol(::Type{Float64}) = (rtol=1e-8, atol=1e-8)
    crtutol(::Type{BigFloat}) = (rtol=1e-12, atol=1e-12)

    # The type of the real numbers in the arguments
    realtypeof(a::Number) = float(real(typeof(a)))
    realtypeof(a::Integer) = Union{}
    realtypeof(a::AbstractQuaternion) = realtypeof(a[1])
    realtypeof(a::AbstractArray{<:Number}) = realtypeof(first(a))
    realtypeof(a) = Union{}
    realtypeof(args::Tuple) = mapreduce(realtypeof, promote_type, args; init=Union{})

    # Attach explicit random tangents to the arguments with complex components
    withtangent(a) = a
    withtangent(a::AbstractQuaternion{<:Complex}) = a ⊢ rand_tangent(rng, a)
    withtangent(a::Complex) = a ⊢ rand_tangent(rng, a)

    """
        crtu(f, args...; check_inferred=true, frule=true, kwargs...)

    Run ChainRulesTestUtils' `test_rrule` (and `test_frule`, unless `frule` is `false`) on
    `f(args...)`, with tolerances for the precision of the arguments and explicit tangents
    for arguments with complex components.
    """
    function crtu(f, args...; check_inferred=true, frule=true, kwargs...)
        reseed!()
        tol = crtutol(realtypeof(args))
        targs = map(withtangent, args)
        test_rrule(f, targs...; check_inferred=check_inferred, tol..., kwargs...)
        frule && test_frule(f, targs...; check_inferred=check_inferred, tol..., kwargs...)
    end

    # A quaternion with complex components built from a four-entry point, with fixed
    # imaginary parts
    QCx(x) = QC([x; [0.2, -0.1, 0.3, 0.05]])

    # The tolerance of the comparison with the BigFloat references at a point.  Near -1 the
    # derivatives of `log`, `sqrt`, and the powers are of order 1e7, and roundoff in the
    # Float64 rules grows with them.
    pointtol(name) = name === :nearminus1 ? 1e-7 : 1e-10
end


@testitem "ChainRules: structural rules" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    # Constructors with explicit types store their arguments without normalizing them, so
    # their references are the same constructors with the element type of the (possibly
    # BigFloat) perturbed arguments.
    numtype(x::Number) = typeof(x)
    numtype(x::AbstractArray) = eltype(x)
    numtype(x::AbstractQuaternion) = typeof(x[1])
    storing(QT) = (a...) -> QT{promote_type(map(numtype, a)...)}(a...)

    for T ∈ (Float64, Float32, BigFloat)
        rtol = T === Float32 ? 1e-5 : T === BigFloat ? 1e-40 : 1e-12
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        @testset "Constructors $T" begin
            for QT ∈ (Quaternion, Rotor, QuatVec)
                for args ∈ ((w, x, y, z), (x, y, z), (w,), (SVector(w, x, y, z),), ([w, x, y, z],),
                            ([x, y, z],), ([w],), (Quaternion{T}(w, x, y, z),), (Rotor{T}(w, x, y, z),))
                    # `Rotor{T}`, `Quaternion{T}`, and `QuatVec{T}` have no methods for vectors
                    # of lengths 3 and 1.
                    if !(args[1] isa AbstractVector && length(args[1]) < 4)
                        checkrule(QT{T}, args...; ref=storing(QT), rtol=rtol)
                    end
                    checkrule(QT, args...; rtol=rtol)
                end
            end
            for f ∈ (quaternion, rotor, quatvec)
                for args ∈ ((w, x, y, z), (x, y, z), (w,), (SVector(w, x, y, z),), ([w, x, y, z],),
                            ([x, y, z],), ([w],), (Quaternion{T}(w, x, y, z),), (Rotor{T}(w, x, y, z),))
                    checkrule(f, args...; rtol=rtol)
                end
                @test_throws DimensionMismatch rrule(f, T[1, 2])
            end
            # Integer arguments of a constructor get the floating-point cotangents of the
            # equivalent floats.
            Δ = quaternion(T(3//10), T(-11//10), T(7//10), T(1//2))
            for (f, ints) ∈ ((quaternion, (1, 2, 3, 4)), (rotor, (1, 2, 3, 4)), (quatvec, (2, 3, 4)), (rotor, (2, 3, 4)))
                ∂ints = rrule(f, ints...)[2](Δ)
                ∂floats = rrule(f, map(float, ints)...)[2](Δ)
                @test all(map((a, b) -> unthunk(a) isa AbstractFloat && unthunk(a) ≈ unthunk(b), ∂ints[2:end], ∂floats[2:end]))
            end
        end
    end

    @testset "Constructors with complex arguments" begin
        z1, z2, z3, z4 = 0.3 + 0.1im, -0.5 + 0.2im, 0.7 - 0.4im, 0.2 + 0.6im
        for f ∈ (quaternion, Quaternion, Quaternion{ComplexF64}, quatvec, QuatVec, Rotor{ComplexF64})
            ref = f isa Type && !(f isa UnionAll) ? storing(f.name.wrapper) : f
            checkrule(f, z1, z2, z3, z4; ref=ref)
            checkrule(f, [z1, z2, z3, z4]; ref=ref)
        end
        checkrule(rotor, z1, z2, z3, z4)
    end

    @testset "Structural access at $name" for (name, xp) ∈ points()
        tol = pointtol(name)
        for q ∈ (Q(xp), R(xp), Ru(xp), V(xp), QCx(xp))
            checkrule(components, q; rtol=tol)
            for i ∈ 1:4
                checkrule(getindex, q, i; rtol=tol)
            end
            checkrule(real, q; rtol=tol)
            checkrule(imag, q; rtol=tol)
            checkrule(vec, q; rtol=tol)
        end
    end

    @testset "rotor and normalize at $name" for (name, xp) ∈ points()
        tol = pointtol(name)
        checkrule(normalize, Q(xp); rtol=tol)
        checkrule(normalize, Ru(xp); rtol=tol)
        checkrule(normalize, QCx(xp); rtol=tol)
        checkrule(rotor, Q(xp); rtol=tol)
        checkrule(rotor, Ru(xp); rtol=tol)
        checkrule(Rotor, Q(xp); rtol=tol)
        checkrule(rotor, xp...; rtol=tol)
        checkrule(rotor, xp; rtol=tol)
        if !iszero(xp[2:4])
            checkrule(rotor, xp[2:4]...; rtol=tol)
            checkrule(normalize, V(xp); rtol=tol)
        end
        iszero(xp[1]) || checkrule(rotor, xp[1]; rtol=tol)
    end

    @testset "ChainRulesTestUtils, $T" for T ∈ (Float64, Float32, BigFloat)
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        sv, v4, v3, v1 = @SVector[w, x, y, z], [w, x, y, z], [x, y, z], [w]
        for QT ∈ (Quaternion, Rotor, QuatVec)
            # This list is ported from the retired `auto_differentiation.jl`.
            for args ∈ ((sv,), (v4,), (w, x, y, z), (x, y, z), (w,))
                crtu(QT{T}, args...)
                crtu(QT, args...)
            end
            for args ∈ ((v3,), (v1,))
                crtu(QT, args...)
            end
        end
        for f ∈ (quaternion, rotor, quatvec), args ∈ ((sv,), (v4,), (v3,), (v1,), (w, x, y, z), (x, y, z), (w,))
            crtu(f, args...)
        end
        for q ∈ (Quaternion{T}(w, x, y, z), Rotor{T}(w, x, y, z), QuatVec{T}(x, y, z))
            crtu(components, q)
            crtu(real, q)
            crtu(imag, q)
            crtu(vec, q)
            for i ∈ 1:4
                crtu(getindex, q, i)
            end
            crtu(normalize, q)
            q isa QuatVec || crtu(rotor, q)
        end
    end
end


@testitem "ChainRules: arithmetic" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    import Base.FastMath: add_fast, sub_fast, mul_fast, div_fast, inv_fast
    p₀ = quaternion(0.3, 0.8, -0.4, 1.1)
    pc = QC([-0.2, 0.1, 0.6, -0.3, 0.4, 0.3, -0.1, 0.2])
    zc = 0.7 + 0.2im

    @testset "Binary operations at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        q, ru, v, qc = Q(x), Ru(x), V(x), QCx(x)
        nonzerovec = !iszero(x[2:4])
        pairs = Any[
            (q, p₀), (p₀, q), (q, 0.7), (0.7, q), (q, 2), (v, p₀), (p₀, v), (v, V0), (V0, v),
            (v, 0.7), (0.7, v), (ru, p₀), (p₀, ru), (ru, R0), (R0, ru), (ru, Ru(x[[2, 1, 4, 3]])),
            (0.7, ru), (ru, 0.7), (v, ru), (ru, v), (q, ru), (qc, pc), (pc, qc), (qc, zc),
            (zc, qc), (qc, 0.7), (q, pc), (pc, q),
        ]
        for f ∈ (+, -, *, /, \, add_fast, sub_fast, mul_fast, div_fast), (a, b) ∈ pairs
            # Division by a `QuatVec` with a zero vector part is division by zero.
            divisor = f ∈ (/, div_fast) ? b : f === (\) ? a : nothing
            divisor isa QuatVec && !nonzerovec && continue
            # The fast-math operations are checked with quaternion pairs and scalars only.
            if f ∈ (add_fast, sub_fast, mul_fast, div_fast)
                (a isa Rotor || b isa Rotor || a isa QuatVec || b isa QuatVec) && continue
            end
            checkrule(f, a, b; rtol=tol)
        end
    end

    @testset "Unary operations at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for q ∈ (Q(x), R(x), Ru(x), V(x), QCx(x))
            iszero(q) && continue
            for f ∈ (-, conj, inv, sub_fast, inv_fast)
                checkrule(f, q; rtol=tol)
            end
        end
    end

    @testset "muladd at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        q, ru = Q(x), Ru(x)
        for (a, b, c) ∈ ((q, p₀, P0), (q, p₀, 0.4), (q, 0.4, p₀), (0.4, q, p₀), (q, 0.4, 0.3),
                         (0.4, q, 0.3), (0.4, 0.3, q), (ru, p₀, R0), (QCx(x), pc, zc))
            checkrule(muladd, a, b, c; rtol=tol)
        end
    end

    @testset "n-ary products and sums at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        q, ru = Q(x), Ru(x)
        # The frules of the n-ary operations are opted out (the 3-argument frule of
        # ChainRules multiplies in the wrong order), so only the pullbacks are checked.
        for f ∈ (*, +, mul_fast, add_fast)
            for args ∈ ((q, p₀, P0), (2.0, q, p₀), (q, 2.0, p₀), (q, p₀, 2.0), (2.0, 3.0, q),
                        (2.0, 3.0, q, p₀), (q, p₀, q, P0), (0.5, q, 2.0, p₀, 3.0), (ru, R0, ru),
                        (QCx(x), pc, zc, 0.3))
                checkrule(f, args...; rtol=tol, pushforward=false)
            end
        end
        # ChainRules' own rule for four or more factors reaches the fold rules through its
        # recursion, also when the first quaternion is the fifth factor or later.
        checkrule(*, 2.0, 3.0, 4.0, 5.0, q; rtol=tol, pushforward=false)
        checkrule(*, 2.0, 3.0, 4.0, 5.0, 6.0, q, p₀; rtol=tol, pushforward=false)
        checkrule(*, 0.3 + 0.1im, 2.0, 0.5im, QCx(x); rtol=tol, pushforward=false)
        checkrule(*, 0.3 + 0.1im, 2.0, 0.5im, q, ru; rtol=tol, pushforward=false)
        # The known gap of section 2 of the specification: when the first quaternion of a
        # sum is its fourth argument or later, ChainRules' unprojected n-ary rule is reached,
        # and the real arguments receive quaternion cotangents.  Their scalar parts are
        # correct.
        Δ = quaternion(0.3, -1.1, 0.7, 0.5)
        ∂ = map(unthunk, rrule(+, 2.0, 3.0, 4.0, q)[2](Δ))
        scalarpart(x) = x isa AbstractQuaternion ? x[1] : x
        @test all(i -> scalarpart(∂[i]) ≈ Δ[1], 2:4)
        @test ∂[5] ≈ Δ
    end

    @testset "ChainRulesTestUtils, $T" for T ∈ (Float64, Float32, BigFloat)
        xs = (T[1.2, -0.7, 0.5, 0.3], T[0.5, 0.3, -0.2, 0.4], T[-0.8, 0.6, 0.5, -0.7])
        q, p, r1, r2, v, t = Q(xs[1]), Q(xs[2]), Ru(xs[2]), Ru(xs[3]), V(xs[3]), T(7//10)
        for f ∈ (+, -, *, /, \)
            for (a, b) ∈ ((q, p), (q, t), (t, q), (v, q), (q, v), (r1, r2), (r1, q), (q, r1),
                          (t, r1), (v, r1), (r1, v))
                crtu(f, a, b)
            end
        end
        for f ∈ (-, conj, inv), a ∈ (q, r1, v)
            crtu(f, a)
        end
        crtu(muladd, q, p, r1)
        crtu(muladd, t, q, p)
        crtu(*, T(2), T(3), q, p; frule=false)
        crtu(+, T(2), q, T(3), p; frule=false)
        crtu(*, q, p, r1, v; frule=false)
        if T === Float64
            qc = QCx(xs[1])
            pc′ = QCx(xs[2])
            for f ∈ (+, -, *, /), (a, b) ∈ ((qc, pc′), (qc, zc), (zc, qc), (q, pc′), (qc, 0.7))
                crtu(f, a, b)
            end
            for f ∈ (-, conj, inv)
                crtu(f, qc)
            end
        end
    end
end


@testitem "ChainRules: exp, log, and sqrt" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    import Base.FastMath: exp_fast, log_fast, sqrt_fast

    @testset "exp at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for q ∈ (Q(x), V(x), Ru(x), QCx(x), L([x[1] / 2; x[2:4]]))
            checkrule(exp, q; rtol=tol)
        end
        checkrule(exp_fast, Q(x); rtol=tol)
    end
    @testset "exp near zero vector parts" begin
        for ϵ ∈ (1e-3, 1e-6, 1e-12, 0.0)
            checkrule(exp, quaternion(0.4, ϵ, -2ϵ, 3ϵ))
            checkrule(exp, quatvec(ϵ, -2ϵ, 3ϵ))
            checkrule(exp, QC([0.4, ϵ, 0, -ϵ, 0.1, 0, 2ϵ, 0]))
        end
        # On both sides of the series threshold (|v⃗|² = 1)
        for a ∈ (0.999, 1.001)
            checkrule(exp, quaternion(0.4, a, 0, 0))
            checkrule(exp, quatvec(a / √2, -a / √2, 0))
        end
    end

    # `log` and `sqrt` of a `Rotor` evaluate formulas that agree only on the unit sphere on
    # the two sides of w = 0, so a finite difference across w = 0 is meaningless unless the
    # stored norm is exactly 1.  The pure-vector point of the point set is normalized in
    # floating point, so it is excluded for rotors, and the rotors at w = 0 are instead
    # checked at points whose stored norm is exactly 1.
    @testset "log and sqrt of rotors at w = 0" begin
        for q ∈ (Rotor{Float64}(0, 1, 0, 0), Rotor{Float64}(0, 0, -1, 0), Rotor{Float64}(0, 0, 0, 1),
                 Rotor{Float64}(0, 0, 0, -1), Rotor{BigFloat}(0, 0, 1, 0))
            checkrule(log, q)
            checkrule(sqrt, q)
        end
    end
    @testset "log at $name" for (name, x) ∈ points(:nocut)
        tol = pointtol(name)
        checkrule(log, Q(x); rtol=tol)
        checkrule(log_fast, Q(x); rtol=tol)
        name === :purevector || checkrule(log, Ru(x); rtol=tol)
        name === :purevector || checkrule(log, R(x); rtol=tol)
    end
    @testset "log near the positive real axis" begin
        # On both sides of the series threshold a² = w²/16
        for a ∈ (0.2499, 0.2501, 1e-4, 1e-9)
            checkrule(log, quaternion(1.0, a, 0, 0))
            checkrule(log, quaternion(2.0, 0, 2a, 0))
            checkrule(log, Rotor{Float64}(1.0, a, 0, 0))
        end
    end
    @testset "sqrt at $name" for (name, x) ∈ points(:nocut)
        tol = pointtol(name)
        checkrule(sqrt, Q(x); rtol=tol)
        checkrule(sqrt_fast, Q(x); rtol=tol)
        iszero(x[2:4]) || checkrule(sqrt, V(x); rtol=tol)
        name === :purevector || checkrule(sqrt, Ru(x); rtol=tol)
        name === :purevector || checkrule(sqrt, R(x); rtol=tol)
    end
    @testset "sqrt at extreme magnitudes" begin
        # The source rescales by powers of two where the squared norm would overflow or
        # underflow.  The finite-difference step of the references is absolute, so instead
        # the rules are checked against the homogeneity √(s q) = √s √q, with s a power of 2,
        # which makes the Jacobian at s q equal to the Jacobian at q divided by √s.
        Δ = quaternion(0.3, -1.1, 0.7, 0.5)
        for q₀ ∈ (quaternion(1.2, -0.7, 0.5, 0.3), quaternion(-1.2, -0.7, 0.5, 0.3)), k ∈ (-1000, -560, -100, 100, 560, 1000)
            s = 2.0^k
            q = s * q₀
            ∂ = unthunk(rrule(sqrt, q)[2](Δ)[2])
            ∂₀ = unthunk(rrule(sqrt, q₀)[2](Δ)[2])
            @test ∂ ≈ ∂₀ / sqrt(s) rtol = 1e-14
            @test frule((NoTangent(), Δ), sqrt, q)[2] ≈ frule((NoTangent(), Δ), sqrt, q₀)[2] / sqrt(s) rtol = 1e-14
        end
    end

    # On the negative real axis, the source returns log|w| + π𝐤 and √(-w) 𝐤, which depend
    # only on w.  The rules differentiate exactly that: the pushforward along w matches the
    # finite difference along w, and the pullback has a zero vector part.
    @testset "On the negative real axis" begin
        e₁ = quaternion(1.0, 0, 0, 0)
        Δ = quaternion(0.3, -1.1, 0.7, 0.5)
        for f ∈ (log, sqrt), q ∈ (quaternion(-1.0, 0, 0, 0), quaternion(-4.0, 0, 0, 0), Rotor{Float64}(-1.0, 0, 0, 0))
            @testset "$f($q)" begin
                Ω, Ω̇ = frule((NoTangent(), e₁), f, q)
                @test Ω ≈ f(q)
                Ω̇ref = bigjvp(f, q, e₁)
                @test relerr(Ω̇, Ω̇ref) < 1e-12
                ∂ = unthunk(rrule(f, q)[2](Δ)[2])
                if ∂ isa AbstractZero
                    @test relerr(zero(Ω̇ref), Ω̇ref) < 1e-12
                else
                    @test ∂ isa Quaternion
                    @test ∂[1] ≈ Ω̇ref ⋅ Δ atol=1e-14
                    @test iszero(∂[2]) && iszero(∂[3]) && iszero(∂[4])
                end
            end
        end
        # log of a `Rotor` on the negative real axis is the constant π𝐤, whatever its norm.
        for q ∈ (Rotor{Float64}(-1.0, 0, 0, 0), Rotor{Float64}(-2.0, 0, 0, 0))
            ∂ = unthunk(rrule(log, q)[2](quatvec(0.3, -1.1, 0.7))[2])
            @test ∂ isa AbstractZero || iszero(∂)
        end
    end

    @testset "ChainRulesTestUtils, $T" for T ∈ (Float64, Float32, BigFloat)
        for xp ∈ (T[1.2, -0.7, 0.5, 0.3], T[0.5, 0.3, -0.2, 0.4], T[-0.8, 0.6, 0.5, -0.7])
            crtu(exp, Q(xp))
            crtu(exp, V(xp))
            crtu(exp, Ru(xp))
            crtu(log, Q(xp))
            crtu(log, Ru(xp))
            crtu(sqrt, Q(xp))
            crtu(sqrt, V(xp))
            crtu(sqrt, Ru(xp))
        end
        T === Float64 && crtu(exp, QCx([1.2, -0.7, 0.5, 0.3]))
        T === Float64 && crtu(exp, L([0.3, -0.7, 0.5, 0.3]))
    end
    # Ported from the retired `auto_differentiation.jl`
    @testset "exp with given output tangents" begin
        w, x, y, z = 1.2, 3.4, 5.6, 7.8
        crtu(exp, Quaternion(w, x, y, z); output_tangent=Quaternion(x, y, z, w))
        crtu(exp, QuatVec(x, y, z); output_tangent=Quaternion(x, y, z, w))
    end
end


@testitem "ChainRules: powers" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    import Base.FastMath: pow_fast

    @testset "Integer powers at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for n ∈ -3:3
            checkrule(^, Q(x), n; rtol=tol)
            checkrule(^, Ru(x), n; rtol=tol)
            checkrule(^, R(x), n; rtol=tol)
            iszero(x[2:4]) || checkrule(^, V(x), n; rtol=tol)
            checkrule(^, QCx(x), n; rtol=tol)
            checkrule(pow_fast, Q(x), n; rtol=tol)
            checkrule(Base.literal_pow, ^, Q(x), Val(n); rtol=tol)
            checkrule(Base.literal_pow, ^, Ru(x), Val(n); rtol=tol)
        end
    end

    @testset "Real powers at $name" for (name, x) ∈ points(:nocut)
        tol = pointtol(name)
        for s ∈ (0.3, -1.7, 2//3)
            checkrule(^, Q(x), s; rtol=tol)
            checkrule(^, Ru(x), s; rtol=tol)
            checkrule(^, R(x), s; rtol=tol)
            iszero(x[2:4]) || checkrule(^, V(x), s; rtol=tol)
        end
        checkrule(pow_fast, Q(x), 0.3; rtol=tol)
        checkrule(pow_fast, Ru(x), 0.3; rtol=tol)
        # Large exponents of rotors near the identity, beyond the series in `y`
        checkrule(^, Ru(x), 1e5; rtol=1e-8)
    end

    @testset "Quaternion powers at $name" for (name, x) ∈ points(:nocut)
        tol = pointtol(name)
        checkrule(^, Q(x), P0; rtol=tol)
        checkrule(^, Q(x), R0; rtol=tol)
        checkrule(^, Q(x), V0; rtol=tol)
        checkrule(^, Ru(x), P0; rtol=tol)
        checkrule(^, 2.0, Q(x); rtol=tol)
        checkrule(^, 0.7, Ru(x); rtol=tol)
        checkrule(^, P0, Q(x); rtol=tol)
    end

    @testset "Rotor powers at -1" begin
        # At exactly -1, `R^s` is cos(πs) + sin(πs) 𝐤 for any norm, so the cotangent of the
        # rotor is zero, and that of `s` is ⟨π𝐤 Ω, Δ⟩.
        for q ∈ (Rotor{Float64}(-1.0, 0, 0, 0), Rotor{Float64}(-2.0, 0, 0, 0)), s ∈ (0.3, -1.7)
            Ω, pb = rrule(^, q, s)
            @test Ω ≈ q^s
            Δ = quaternion(0.3, -1.1, 0.7, 0.5)
            _, ∂q, ∂s = pb(Δ)
            ∂q = unthunk(∂q)
            @test ∂q isa AbstractZero || iszero(∂q)
            @test unthunk(∂s) ≈ (π * 𝐤 * quaternion(Ω)) ⋅ Δ
            @test unthunk(∂s) ≈ bigvjp(t -> q^t, s, Δ)
            # Integer powers at -1 are exact products, differentiable in every direction.
            checkrule(^, q, 3)
        end
    end

    @testset "ChainRulesTestUtils, $T" for T ∈ (Float64, Float32, BigFloat)
        for xp ∈ (T[1.2, -0.7, 0.5, 0.3], T[0.5, 0.3, -0.2, 0.4], T[-0.8, 0.6, 0.5, -0.7])
            for n ∈ (-2, 0, 1, 3)
                crtu(^, Q(xp), n)
                crtu(^, Ru(xp), n)
                crtu(^, V(xp), n)
            end
            crtu(^, Q(xp), T(3//10))
            crtu(^, Ru(xp), T(3//10))
            crtu(^, V(xp), T(3//10))
            crtu(^, Q(xp), Q(T[0.2, 0.3, -0.1, 0.4]))
            crtu(^, T(2), Q(xp))
            crtu(Base.literal_pow, ^, Q(xp), Val(2))
            crtu(Base.literal_pow, ^, Ru(xp), Val(-1))
        end
    end
end


@testitem "ChainRules: norms and geometry" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    import Base.FastMath: angle_fast

    @testset "Norms at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for q ∈ (Q(x), R(x), Ru(x), V(x), QCx(x), L([x[1] / 2; x[2:4]]))
            checkrule(abs2, q; rtol=tol)
            checkrule(abs2vec, q; rtol=tol)
            checkrule(norm, q; rtol=tol)
            iszero(q) || checkrule(abs, q; rtol=tol)
            # `absvec` has a cone-shaped kink where the vector part vanishes.
            iszero(x[2:4]) || checkrule(absvec, q; rtol=tol)
        end
    end
    @testset "Norms of rotors are constant" begin
        # `abs`, `abs2`, and (for real components) `norm` of a `Rotor` return 1 without
        # reading the stored components, so their derivatives are zero, even off the sphere.
        # The old Zygote tests asserted 2R for `abs2`, which was wrong.
        for q ∈ (Ru([1.2, -0.7, 0.5, 0.3]), R([1.2, -0.7, 0.5, 0.3])), f ∈ (abs, abs2, norm)
            ∂ = unthunk(rrule(f, q)[2](1.0)[2])
            @test ∂ isa AbstractZero || iszero(∂)
            Ω̇ = unthunk(frule((NoTangent(), quaternion(0.3, 0.2, -0.1, 0.4)), f, q)[2])
            @test Ω̇ isa AbstractZero || iszero(Ω̇)
        end
    end
    @testset "Subgradients at cone points" begin
        # Where a norm vanishes, the cotangent is zero, as for `hypot`.
        for q ∈ (quaternion(0.7, 0, 0, 0), quatvec(0.0, 0.0, 0.0), Rotor{Float64}(1, 0, 0, 0))
            ∂ = unthunk(rrule(absvec, q)[2](1.0)[2])
            @test ∂ isa AbstractZero || iszero(∂)
        end
        for q ∈ (quaternion(0.0, 0, 0, 0), quatvec(0.0, 0.0, 0.0))
            ∂ = unthunk(rrule(abs, q)[2](1.0)[2])
            @test ∂ isa AbstractZero || iszero(∂)
        end
        for q ∈ (quaternion(0.7, 0, 0, 0), Rotor{Float64}(1, 0, 0, 0), quaternion(-0.7, 0, 0, 0))
            ∂ = unthunk(rrule(angle, q)[2](1.0)[2])
            @test ∂ isa AbstractZero || iszero(∂)
        end
    end

    @testset "R(v) at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for r ∈ (R(x), Ru(x), L([x[1] / 2; x[2:4]]))
            checkrule(r, V0; rtol=tol)
            iszero(x[2:4]) || checkrule(r, V(x[[4, 3, 2, 1]]); rtol=tol)
        end
        checkrule(R0, V(x); rtol=tol)
    end

    @testset "angle at $name" for (name, x) ∈ points(:nonreal)
        tol = pointtol(name)
        for q ∈ (Q(x), R(x), Ru(x))
            checkrule(angle, q; rtol=tol)
            checkrule(angle_fast, q; rtol=tol)
        end
    end

    @testset "distance2 and distance at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for (a, b) ∈ ((R(x), R0), (R0, R(x)), (Ru(x), Ru([0.3, 0.1, -0.9, 0.2])), (Ru(x), Ru(-x)))
            checkrule(distance2, a, b; rtol=tol)
            distance(a, b) > 1e-6 && checkrule(distance, a, b; rtol=tol)
        end
        # At R₁ == R₂, distance2 is smooth with a zero gradient, and distance has a cone.
        checkrule(distance2, R(x), R(x); rtol=tol)
        checkrule(distance2, Ru(x), Ru(2x); rtol=tol)
        ∂ = rrule(distance, R(x), R(x))[2](1.0)
        @test all(d -> unthunk(d) isa AbstractZero || iszero(unthunk(d)), ∂[2:3])
    end
    @testset "distance2 near the series threshold" begin
        # The source switches to a series in x = a²/w² when x ≤ ∜eps/2.
        x₀ = ∜eps(Float64) / 2
        for x ∈ (0.99x₀, 1.01x₀, 1e-12, 1e-30)
            a = sqrt(x)
            checkrule(distance2, Rotor{Float64}(1.0, a, 0, 0), Rotor{Float64}(1.0, 0, 0, 0))
            checkrule(distance2, rotor(0.5, 0.2, -0.3, 0.1), rotor(0.5, 0.2, -0.3, 0.1) * rotor(1.0, 0, a, 0))
        end
    end

    @testset "dot at $name" for (name, x) ∈ points()
        tol = pointtol(name)
        for (a, b) ∈ ((Q(x), P0), (P0, Q(x)), (Ru(x), R0), (V(x), V0), (Q(x), V0), (QCx(x), QCx(x[[2, 3, 4, 1]])))
            checkrule(dot, a, b; rtol=tol)
        end
    end

    @testset "ChainRulesTestUtils, $T" for T ∈ (Float64, Float32, BigFloat)
        for xp ∈ (T[1.2, -0.7, 0.5, 0.3], T[0.5, 0.3, -0.2, 0.4], T[-0.8, 0.6, 0.5, -0.7])
            for q ∈ (Q(xp), Ru(xp), V(xp))
                for f ∈ (abs2, abs2vec, abs, absvec, norm)
                    crtu(f, q)
                end
            end
            crtu(angle, Q(xp))
            crtu(angle, Ru(xp))
            crtu(Ru(xp), V(T[0, 0.3, -0.6, 0.2]))
            crtu(distance2, Ru(xp), Ru(T[0.3, 0.1, -0.9, 0.2]))
            crtu(distance, Ru(xp), Ru(T[0.3, 0.1, -0.9, 0.2]))
            crtu(dot, Q(xp), Ru(T[0.3, 0.1, -0.9, 0.2]))
        end
        if T === Float64
            for f ∈ (abs2, abs, norm, abs2vec, absvec)
                crtu(f, QCx([1.2, -0.7, 0.5, 0.3]))
            end
            crtu(dot, QCx([1.2, -0.7, 0.5, 0.3]), QCx([0.5, 0.3, -0.2, 0.4]))
            crtu(L([0.3, -0.7, 0.5, 0.3]), V0)
        end
    end
end


@testitem "ChainRules: primitives that backends cannot trace" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    # `csqrt` is the complex square root, isolated so that Mooncake can import its rule.
    for z ∈ (0.3 + 0.4im, -0.7 + 0.1im, -0.7 - 0.1im, 2.0 + 0.0im, 0.3f0 + 0.4f0im, complex(big"0.3", big"0.4"))
        checkrule(Quaternionic.csqrt, z; rtol=max(eps(real(typeof(z)))^(3 / 4), 1e-40))
        crtu(Quaternionic.csqrt, z)
    end

    # `dominant_eigenvector_lapack(A, uplo)` reads only the triangle of `A` named by `uplo`.
    # The rules are checked with an asymmetric matrix, so that reading the wrong triangle is
    # detected, at the Bar-Itzhack matrices of a generic rotation, of the identity, and of a
    # rotation by π, whose dominant eigenvalues are simple.
    using LinearAlgebra: Symmetric, I
    using GenericLinearAlgebra
    function baritzhack(R::Rotor)
        M = to_rotation_matrix(R)
        K = zeros(4, 4)
        K[1, 1] = (M[1, 1] - M[2, 2] - M[3, 3]) / 3
        K[2, 2] = (M[2, 2] - M[1, 1] - M[3, 3]) / 3
        K[3, 3] = (M[3, 3] - M[1, 1] - M[2, 2]) / 3
        K[4, 4] = (M[1, 1] + M[2, 2] + M[3, 3]) / 3
        K[1, 2] = K[2, 1] = (M[2, 1] + M[1, 2]) / 3
        K[1, 3] = K[3, 1] = (M[3, 1] + M[1, 3]) / 3
        K[2, 3] = K[3, 2] = (M[3, 2] + M[2, 3]) / 3
        K[1, 4] = K[4, 1] = (M[2, 3] - M[3, 2]) / 3
        K[2, 4] = K[4, 2] = (M[3, 1] - M[1, 3]) / 3
        K[3, 4] = K[4, 3] = (M[1, 2] - M[2, 1]) / 3
        K
    end
    rng′ = MersenneTwister(7)
    for R ∈ (R0, rotor(1.0, 0, 0, 0), rotor(0, 1.0, 0, 0)), uplo ∈ ('U', 'L'), T ∈ (Float64, Float32)
        A = T.(baritzhack(R) + 0.01 * randn(rng′, 4, 4))
        f(A) = Quaternionic.dominant_eigenvector_lapack(A, uplo)
        # The finite differences of ChainRulesTestUtils are too coarse in Float32 for LAPACK.
        # Julia 1.10 cannot infer the type of the cotangent of `A` (the values are correct),
        # so the inference check runs only on later versions.
        T === Float64 && crtu(
            Quaternionic.dominant_eigenvector_lapack, A, uplo ⊢ NoTangent();
            check_inferred=VERSION ≥ v"1.11"
        )
        # The eigenvector is defined up to its sign, so the reference uses the eigenvector
        # with the sign of the primal.
        v = f(A)
        ȳ = randn(rng′, T, 4)
        _, pb = rrule(Quaternionic.dominant_eigenvector_lapack, A, uplo)
        _, Ā, ū = pb(ȳ)
        @test ū isa NoTangent
        Ā = unthunk(Ā)
        @test Ā isa Matrix{T}
        # The entries outside the triangle named by `uplo` have zero cotangents.
        outside = uplo == 'U' ? [i > j for i ∈ 1:4, j ∈ 1:4] : [i < j for i ∈ 1:4, j ∈ 1:4]
        @test all(iszero, Ā[outside])
        # The references take the matrix as a vector of its entries, and they compute the
        # eigenvector of the BigFloat matrix with GenericLinearAlgebra.
        signed(a) = (v′ = Quaternionic.dominant_eigenvector(Symmetric(reshape(a, 4, 4), Symbol(uplo)));
                     sign(v′ ⋅ v) * v′)
        @test relerr(vec(Ā), bigvjp(signed, vec(A), ȳ)) < (T === Float32 ? 1e-4 : 1e-10)
        Ȧ = randn(rng′, T, 4, 4)
        _, v̇ = frule((NoTangent(), Ȧ, NoTangent()), Quaternionic.dominant_eigenvector_lapack, A, uplo)
        @test relerr(v̇, bigjvp(signed, vec(A), vec(Ȧ))) < (T === Float32 ? 1e-4 : 1e-10)
    end
end


@testitem "ChainRules: cotangent forms" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    # Backends pass cotangents of quaternion outputs in many forms: Zygote's `real` rule
    # passes a number, `q.components` gives a `Tangent` or `NamedTuple` of the `SVector`'s
    # data, and arrays, tuples, thunks, `nothing`, and stray `Rotor`s or `QuatVec`s occur.
    # Every pullback must accept each form and give the same result as for the equivalent
    # `Quaternion`, without renormalizing a `Rotor`-typed cotangent.
    a, b, c, d = 0.3, -1.1, 0.7, 0.5
    Δ = quaternion(a, b, c, d)
    forms(T) = [
        Δ => Δ,
        @thunk(Δ) => Δ,
        Rotor{Float64}(a, b, c, d) => Δ,
        quatvec(b, c, d) => quaternion(0.0, b, c, d),
        a => quaternion(a, 0.0, 0.0, 0.0),
        Tangent{T}(components=Tangent{SVector{4,Float64}}(data=(a, b, c, d))) => Δ,
        (components=(data=(a, b, c, d),),) => Δ,
        (a, b, c, d) => Δ,
        (a, nothing, c, ZeroTangent()) => quaternion(a, 0.0, c, 0.0),
        [a, b, c, d] => Δ,
        SVector(a, b, c, d) => Δ,
    ]
    zeroforms = (nothing, ZeroTangent(), NoTangent(), @thunk(ZeroTangent()), (nothing, nothing, nothing, nothing))
    same(x, y) = (x = unthunk(x); y = unthunk(y);
                  (x isa AbstractZero && y isa AbstractZero) || (!(x isa AbstractZero) && !(y isa AbstractZero) && typeof(x) == typeof(y) && x ≈ y))
    iszerocot(x) = (x = unthunk(x); x isa AbstractZero)

    q, p, v = quaternion(1.2, -0.7, 0.5, 0.3), quaternion(0.3, 0.8, -0.4, 1.1), quatvec(0.3, -0.6, 0.2)
    r1, r2 = Rotor{Float64}(0.5, 0.3, -0.2, 0.4), Rotor{Float64}(-0.8, 0.6, 0.5, -0.7)
    cases = Any[
        (exp, (q,)), (exp, (v,)), (log, (q,)), (sqrt, (q,)), (inv, (q,)), (-, (q,)), (conj, (r1,)),
        (*, (q, p)), (/, (q, p)), (*, (2.0, q)), (+, (q, 2.0)), (*, (r1, r2)), (/, (q, r1)),
        (^, (q, 3)), (^, (q, 0.3)), (^, (r1, 0.3)), (^, (r1, 2)), (normalize, (q,)), (rotor, (q,)),
        (quaternion, (1.2, -0.7, 0.5, 0.3)), (Quaternion{Float64}, (1.2, -0.7, 0.5, 0.3)),
        (Rotor{Float64}, (1.2, -0.7, 0.5, 0.3)), (rotor, (1.2, -0.7, 0.5, 0.3)),
        (quaternion, ([1.2, -0.7, 0.5, 0.3],)), (r1, (v,)), (muladd, (q, p, r1)), (*, (q, p, r1)),
    ]
    for (f, args) ∈ cases
        @testset "$f$(map(typeof, args))" begin
            Ω, pb = rrule(f, args...)
            T = typeof(Ω)
            for (form, equivalent) ∈ forms(T)
                # A `QuatVec` output accepts the vector part of any form.
                if Ω isa QuatVec
                    equivalent = quaternion(0.0, equivalent[2], equivalent[3], equivalent[4])
                end
                expected = pb(equivalent)
                @test all(map(same, pb(form), expected))
            end
            for form ∈ zeroforms
                @test all(iszerocot, pb(form))
            end
        end
    end

    # Scalar outputs accept thunks, `nothing`, and zeros.
    for (f, args) ∈ ((abs2, (q,)), (abs, (q,)), (absvec, (q,)), (angle, (q,)), (dot, (q, p)),
                     (distance2, (r1, r2)), (real, (q,)), (getindex, (q, 2)))
        pbs = rrule(f, args...)[2]
        @test all(map(same, pbs(@thunk(0.7)), pbs(0.7)))
        for form ∈ (nothing, ZeroTangent(), NoTangent())
            @test all(iszerocot, pbs(form))
        end
    end
    # A complex scalar output accepts the `Tangent` of a complex number.
    pbc = rrule(abs2, QCx([1.2, -0.7, 0.5, 0.3]))[2]
    @test all(map(same, pbc(Tangent{ComplexF64}(re=0.7, im=-0.2)), pbc(0.7 - 0.2im)))
end


@testitem "ChainRules: projection" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    q, r, v = quaternion(1.2, -0.7, 0.5, 0.3), Rotor{Float64}(0.5, 0.3, -0.2, 0.4), quatvec(0.3, -0.6, 0.2)
    Δ = quaternion(0.3, -1.1, 0.7, 0.5)

    # The tangent of a `Quaternion` or `Rotor` is a `Quaternion`, never a renormalized `Rotor`.
    @test ProjectTo(q)(Δ) === Δ
    @test ProjectTo(r)(Rotor{Float64}(0.3, -1.1, 0.7, 0.5)) === Δ
    @test ProjectTo(r)(Δ) === Δ
    # The tangent of a `QuatVec` is a `QuatVec`, with any scalar part discarded.
    @test ProjectTo(v)(Δ) === quatvec(-1.1, 0.7, 0.5)
    @test ProjectTo(v)(Δ) isa QuatVec{Float64}
    # The element type follows the precision of the primal, and integers get floats.
    @test ProjectTo(quaternion(1.2f0, 0, 0, 0))(Δ) isa Quaternion{Float32}
    @test ProjectTo(quaternion(big"1.2", 0, 0, 0))(big"0.5" * Δ) isa Quaternion{BigFloat}
    @test ProjectTo(quaternion(1, 2, 3, 4))(Δ) === Δ
    @test ProjectTo(quatvec(1, 2, 3))(Δ) === quatvec(-1.1, 0.7, 0.5)
    # Real components project away the imaginary parts of complex cotangents; complex
    # components keep them.
    Δc = quaternion(0.3 + 0.1im, -1.1, 0.7im, 0.5)
    @test ProjectTo(q)(Δc) === quaternion(0.3, -1.1, 0.0, 0.5)
    @test ProjectTo(QCx([1.2, -0.7, 0.5, 0.3]))(Δc) === quaternion(0.3 + 0.1im, -1.1 + 0im, 0.0 + 0.7im, 0.5 + 0im)
    # Quaternions with `Bool` components, such as 𝐢, are not differentiable.
    @test ProjectTo(𝐢) isa ProjectTo{NoTangent}
    @test ProjectTo(imz)(Δ) isa NoTangent
    @test ProjectTo(quaternion(true, false, false, false))(Δ) isa NoTangent
    # A real or complex primal that receives a quaternion cotangent keeps its scalar part.
    @test ProjectTo(0.7)(Δ) === 0.3
    @test ProjectTo(0.7f0)(Δ) === 0.3f0
    @test ProjectTo(0.7 + 0.1im)(Δc) === 0.3 + 0.1im
    @test ProjectTo(0.7)(Tangent{typeof(q)}(components=Tangent{SVector{4,Float64}}(data=(0.3, 1.0, 2.0, 3.0)))) === 0.3
    @test ProjectTo(true)(Δ) isa NoTangent
    # Other forms of quaternion cotangents
    @test ProjectTo(q)(0.7) === quaternion(0.7, 0.0, 0.0, 0.0)
    @test ProjectTo(q)(Tangent{ComplexF64}(re=0.7, im=0.1)) === quaternion(0.7, 0.0, 0.0, 0.0)
    @test ProjectTo(q)((0.3, -1.1, 0.7, 0.5)) === Δ
    @test ProjectTo(q)((components=(data=(0.3, -1.1, 0.7, 0.5),),)) === Δ
    @test ProjectTo(q)(Tangent{typeof(q)}(components=Tangent{SVector{4,Float64}}(data=(0.3, -1.1, 0.7, 0.5)))) === Δ
    @test ProjectTo(q)([0.3, -1.1, 0.7, 0.5]) === Δ
    @test ProjectTo(v)([0.3, -1.1, 0.7, 0.5]) === quatvec(-1.1, 0.7, 0.5)
    @test ProjectTo(q)(ZeroTangent()) isa AbstractZero
    @test ProjectTo(q)(NoTangent()) isa AbstractZero
    # An array of quaternions projects elementwise.
    @test ProjectTo([q, r])([Δ, Rotor{Float64}(0.3, -1.1, 0.7, 0.5)]) == [Δ, Δ]
end


@testitem "ChainRules: integer and Bool primals" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    # ChainRulesTestUtils cannot perturb quaternions with integer components (its finite
    # differences rebuild them with integer components and throw an `InexactError`), so
    # these rules are checked against the BigFloat references only.  Integer quaternions
    # have `Float64` tangents, while quaternions with `Bool` components (𝐢, 𝐣, and 𝐤) and
    # integer or `Bool` scalars have `NoTangent()`.
    qi, pi_, vi = quaternion(1, -2, 3, 4), quaternion(2, 0, 1, -1), quatvec(1, -2, 2)
    q = quaternion(1.2, -0.7, 0.5, 0.3)
    for f ∈ (+, -, *, /), (a, b) ∈ ((qi, pi_), (qi, q), (q, qi), (qi, 2), (2, qi), (qi, 0.7), (vi, qi),
                                    (𝐢, q), (q, 𝐣), (𝐤, qi), (true, q), (q, true), (vi, 𝐢))
        checkrule(f, a, b)
    end
    for f ∈ (-, conj, inv, exp, abs2, abs, abs2vec, absvec, normalize, rotor, components, real, imag, vec)
        checkrule(f, qi)
        checkrule(f, vi)
    end
    for f ∈ (log, sqrt, angle)
        checkrule(f, qi)
    end
    for n ∈ (-2, 2, 3)
        checkrule(^, qi, n)
    end
    checkrule(^, qi, 0.3)
    checkrule(getindex, qi, 3)
    checkrule(dot, qi, pi_)
    checkrule(muladd, qi, 2, q)
    checkrule(*, qi, 2, q; pushforward=false)
    checkrule(rotor(qi), vi)
    # A real parameter combined with the constants 𝐢, 𝐣, and 𝐤 gets its own cotangent, and
    # the constants get `NoTangent()`.
    Δ = quaternion(0.3, -1.1, 0.7, 0.5)
    for f ∈ (*, +, -, /)
        _, ∂t, ∂j = rrule(f, 0.7, 𝐣)[2](Δ)
        @test ∂j isa NoTangent
        @test unthunk(∂t) isa Float64
        _, ∂j, ∂t = rrule(f, 𝐣, 0.7)[2](Δ)
        @test ∂j isa NoTangent
        @test unthunk(∂t) isa Float64
    end
end


@testitem "ChainRules: method ambiguities" tags=[:ad, :chainrules] begin
    import ChainRulesCore, ChainRules
    using Test
    import Zygote
    ext = Base.get_extension(Quaternionic, :QuaternionicChainRulesCoreExt)
    @test ext !== nothing
    ext === nothing || @test isempty(Test.detect_ambiguities(ext))
    zext = Base.get_extension(Quaternionic, :QuaternionicZygoteExt)
    @test zext !== nothing
    zext === nothing || @test isempty(Test.detect_ambiguities(zext))
end


@testitem "ChainRules: method tables" tags=[:ad, :chainrules] begin
    # Ported from the retired `auto_differentiation.jl`: every rule must make sense to
    # ChainRulesCore.
    using ChainRulesTestUtils
    import ChainRules
    test_method_tables()
end


@testitem "ChainRules: gradients through constructors under Zygote" tags=[:ad, :chainrules, :zygote] begin
    # Ported from the retired `auto_differentiation.jl`, with every constructor variant.
    # `abs2` of a `Rotor` is the constant 1, whatever its stored components, so its gradient
    # is zero; the old tests asserted Zygote's former answer of 2R.  A `QuatVec` discards its
    # scalar part, so the gradient with respect to that argument is zero or `nothing`.  As
    # in the old tests, the `Quaternion` and `QuatVec` families are also differentiated with
    # symbolic arguments, which are here variables rather than constants.
    using Zygote
    using StaticArrays: @SVector
    using LinearAlgebra: I
    import Symbolics
    zeroish(g) = g === nothing || iszero(g)
    small(T) = g -> g === nothing || abs(g) < 10eps(T)
    same(g, e) = g ≈ e
    same(g::Symbolics.Num, e) = iszero(Symbolics.simplify(g - e; expand=true))
    same(g::Nothing, e) = false
    for T ∈ (BigFloat, Float64, Float32, Float16, Symbolics.Num)
        @testset "$T" begin
            if T === Symbolics.Num
                w, x, y, z = Symbolics.@variables w x y z
            else
                w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
            end
            for f ∈ (
                (a, b, c, d) -> abs2(Quaternion{T}(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(Quaternion{T}([a, b, c, d])),
                (a, b, c, d) -> abs2(Quaternion{T}(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(Quaternion{T}(a, b, c, d)),
                (a, b, c, d) -> abs2(quaternion(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(quaternion([a, b, c, d])),
                (a, b, c, d) -> abs2(quaternion(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(quaternion(a, b, c, d)),
                (a, b, c, d) -> abs2(Quaternion(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(Quaternion([a, b, c, d])),
                (a, b, c, d) -> abs2(Quaternion(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(Quaternion(a, b, c, d)),
            )
                @test all(map(same, Zygote.gradient(f, w, x, y, z), (2w, 2x, 2y, 2z)))
            end
            for f ∈ (
                (a, b, c, d) -> abs2(Quaternion{T}(b, c, d)),
                (a, b, c, d) -> abs2(quaternion(b, c, d)),
                (a, b, c, d) -> abs2(quaternion([b, c, d])),
                (a, b, c, d) -> abs2(Quaternion(b, c, d)),
            )
                ∇ = Zygote.gradient(f, w, x, y, z)
                @test isnothing(∇[1]) && all(map(same, ∇[2:4], (2x, 2y, 2z)))
            end
            for f ∈ (
                (a, b, c, d) -> abs2(Quaternion{T}(a)),
                (a, b, c, d) -> abs2(quaternion(a)),
                (a, b, c, d) -> abs2(quaternion([a])),
                (a, b, c, d) -> abs2(Quaternion(a)),
            )
                ∇ = Zygote.gradient(f, w, x, y, z)
                @test all(isnothing, ∇[2:4]) && same(∇[1], 2w)
            end

            for f ∈ (
                (a, b, c, d) -> abs2(QuatVec{T}(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(QuatVec{T}([a, b, c, d])),
                (a, b, c, d) -> abs2(QuatVec{T}(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(QuatVec{T}(a, b, c, d)),
                (a, b, c, d) -> abs2(quatvec(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(quatvec([a, b, c, d])),
                (a, b, c, d) -> abs2(quatvec(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(QuatVec(@SVector[a, b, c, d])),
                (a, b, c, d) -> abs2(QuatVec([a, b, c, d])),
                (a, b, c, d) -> abs2(QuatVec(Quaternion{T}(@SVector[a, b, c, d]))),
                (a, b, c, d) -> abs2(quatvec(a, b, c, d)),
                (a, b, c, d) -> abs2(QuatVec(a, b, c, d)),
                (a, b, c, d) -> abs2(QuatVec{T}(b, c, d)),
                (a, b, c, d) -> abs2(quatvec(b, c, d)),
                (a, b, c, d) -> abs2(quatvec([b, c, d])),
                (a, b, c, d) -> abs2(QuatVec(b, c, d)),
            )
                ∇ = Zygote.gradient(f, w, x, y, z)
                @test zeroish(∇[1]) && all(map(same, ∇[2:4], (2x, 2y, 2z)))
            end
            for f ∈ (
                (a, b, c, d) -> abs2(QuatVec{T}(a)),
                (a, b, c, d) -> abs2(quatvec([a])),
                (a, b, c, d) -> abs2(quatvec(a)),
                (a, b, c, d) -> abs2(QuatVec(a)),
            )
                @test all(zeroish, Zygote.gradient(f, w, x, y, z))
            end

            # `abs2` of any `Rotor` is constant.  (The rotors are not checked symbolically,
            # because normalization takes square roots.)
            if T !== Symbolics.Num
                for f ∈ (
                    (a, b, c, d) -> abs2(Rotor{T}(@SVector[a, b, c, d])),
                    (a, b, c, d) -> abs2(Rotor{T}([a, b, c, d])),
                    (a, b, c, d) -> abs2(Rotor{T}(Rotor{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(Rotor{T}(Quaternion{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(Rotor{T}(a, b, c, d)),
                    (a, b, c, d) -> abs2(rotor(@SVector[a, b, c, d])),
                    (a, b, c, d) -> abs2(rotor([a, b, c, d])),
                    (a, b, c, d) -> abs2(rotor(Rotor{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(rotor(Quaternion{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(rotor(a, b, c, d)),
                    (a, b, c, d) -> abs2(Rotor(@SVector[a, b, c, d])),
                    (a, b, c, d) -> abs2(Rotor([a, b, c, d])),
                    (a, b, c, d) -> abs2(Rotor(Rotor{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(Rotor(Quaternion{T}(@SVector[a, b, c, d]))),
                    (a, b, c, d) -> abs2(Rotor(a, b, c, d)),
                    (a, b, c, d) -> abs2(Rotor{T}(b, c, d)),
                    (a, b, c, d) -> abs2(rotor(b, c, d)),
                    (a, b, c, d) -> abs2(rotor([b, c, d])),
                    (a, b, c, d) -> abs2(Rotor(b, c, d)),
                    (a, b, c, d) -> abs2(Rotor([b, c, d])),
                    (a, b, c, d) -> abs2(Rotor{T}(a)),
                    (a, b, c, d) -> abs2(rotor(a)),
                    (a, b, c, d) -> abs2(rotor([a])),
                    (a, b, c, d) -> abs2(Rotor(a)),
                    (a, b, c, d) -> abs2(Rotor([a])),
                )
                    @test all(zeroish, Zygote.gradient(f, w, x, y, z))
                end
                # The sum of the squares of the normalized components is 1 to within roundoff,
                # so its gradient is zero to within roundoff.
                for f ∈ (
                    (a, b, c, d) -> sum(abs2, components(rotor(a, b, c, d))),
                    (a, b, c, d) -> sum(abs2, components(rotor(b, c, d))),
                    (a, b, c, d) -> sum(abs2, components(rotor([a]))),
                )
                    @test all(small(T), Zygote.gradient(f, w, x, y, z))
                end
            end

            # Jacobians of the components
            if T ∈ (BigFloat, Float64)
                J = Zygote.jacobian((w, x, y, z) -> components(Quaternion(w, x, y, z)), w, x, y, z)
                @test all(J .== eachrow(Matrix{T}(I, 4, 4)))
                J = Zygote.jacobian((w, x, y, z) -> components(QuatVec(w, x, y, z)), w, x, y, z)
                @test all(J .== (T[0, 0, 0, 0], T[0, 1, 0, 0], T[0, 0, 1, 0], T[0, 0, 0, 1]))
                n = √(w^2 + x^2 + y^2 + z^2)
                J = Zygote.jacobian((w, x, y, z) -> components(Rotor(w, x, y, z)), w, x, y, z)
                c = [w, x, y, z]
                @test all(J[i] ≈ ([j == i for j ∈ 1:4] .- c[i] .* c ./ n^2) ./ n for i ∈ 1:4)
            end
        end
    end
end


@testitem "ChainRules: series thresholds, folds of complex rotors, and irrational bases" tags=[:ad, :chainrules] setup=[ADTestUtils, RuleChecks] begin
    # Near the identity, `R^s` evaluates the cosine and sinc of the angle sθ as series in
    # y = (sθ)², below `trigseriestolerance`, and in closed form above it.  For θ = 0.07, the
    # exponents 1.49 and 1.51 lie on either side of the threshold.  The rules are checked on
    # and off the unit sphere.
    θ = 0.07
    for Rθ ∈ (rotor(cos(θ), (sin(θ) .* NHAT)...), Rotor{Float64}(1.1cos(θ), (1.1sin(θ) .* NHAT)...))
        for s ∈ (1.49, 1.51)
            checkrule(^, Rθ, s)
        end
    end

    # Products of complex rotors are renormalized only when the sum of the squared
    # magnitudes of their components is below 16, which the fold rules for three or more
    # factors and the rules for `muladd` mirror.  `Rbig` is just below the threshold, and its
    # products are above it.  The frules of the folds are opted out.
    Rbig = Rotor{ComplexF64}(2.9 + 0.4im, 1.3 - 1.2im, -1.5 + 0.05im, 0.2 + 1.3im)
    Rsmall = Rotor{ComplexF64}(0.9 + 0.1im, 0.2 - 0.1im, -0.3 + 0.05im, 0.1 + 0.2im)
    checkrule(*, Rbig, Rbig, Rbig; pushforward=false)
    checkrule(*, Rsmall, Rsmall, Rsmall; pushforward=false)
    checkrule(*, Rsmall, Rbig, Rsmall, Rbig; pushforward=false)
    checkrule(muladd, Rbig, Rbig, Rsmall)
    checkrule(muladd, Rsmall, Rsmall, Rbig)
    @test frule((NoTangent(), Rbig, Rbig, Rbig), *, Rbig, Rbig, Rbig) === nothing

    # An irrational base, as in `π^q` and `ℯ^q`, once made the rules throw a `MethodError`.
    q = quaternion(0.3, -0.2, 0.5, 0.1)
    Δ = quaternion(0.3, -1.1, 0.7, 0.5)
    @testset "$b^q" for b ∈ (π, ℯ)
        Ω, pb = rrule(^, b, q)
        @test Ω ≈ b^q
        ∂ = pb(Δ)
        @test unthunk(∂[2]) ≈ bigvjp(t -> t^q, Float64(b), Δ)
        @test relerr(∂[3], bigvjp(p -> b^p, q, Δ), 0) < 1e-12
        Ω₂, Ω̇ = frule((NoTangent(), ZeroTangent(), Δ), ^, b, q)
        @test Ω₂ ≈ b^q
        @test relerr(Ω̇, bigjvp(p -> b^p, q, Δ), 0) < 1e-12
    end

    # `Quaternionic.value` is non-differentiable, and its rule must not be opted out.
    @test rrule(Quaternionic.value, 1.0) !== nothing
end


@testitem "ChainRules: Zygote through irrational bases, complex exp, and literal powers" tags=[:ad, :chainrules, :zygote] setup=[ADTestUtils] begin
    using Quaternionic
    using Test
    using .ADTestUtils
    import Zygote, ForwardDiff

    # Gradients and Hessians of powers with an irrational base, against ForwardDiff on the
    # source and against the BigFloat references
    x4 = [0.3, -0.2, 0.5, 0.1]
    @testset "$b^q" for b ∈ (π, ℯ)
        f = x -> sc(b^quaternion(x...))
        @test relerr(Zygote.gradient(f, x4)[1], ForwardDiff.gradient(f, x4), 0) < 1e-13
        @test relerr(Zygote.gradient(f, x4)[1], bigfd(f, x4), 0) < 1e-13
        @test relerr(Zygote.hessian(f, x4), ForwardDiff.hessian(f, x4), 0) < 1e-12
        @test relerr(Zygote.hessian(f, x4), bigfdH(f, x4), 0) < 1e-12
    end

    # Hessians of `exp` of quaternions with complex components at a pure boost and at a
    # quaternion whose complex vector part has a negative real square, through the rules
    @testset "exp(QC(x)) at $x" for x ∈ ([0.1, 0, 0, 0, 0, 1.3, -0.8, 0.4], [1.0, 0, 0, 0, 0.2, 0, -0.3, 0])
        f = x -> sc(exp(QC(x)))
        @test relerr(Zygote.hessian(f, x), bigfdH(f, x), 0) < 1e-12
    end

    # A literal power of a complex `Rotor`, whose rule returns `nothing` because the power
    # is opted out, so that the Zygote extension differentiates `q^p` instead
    x8 = [0.7, 0.3, -0.2, 0.4, 0.13, 0.1, 0.05, -0.2]
    fpow = x -> sc(rotor(QC(x))^2)
    @test relerr(Zygote.gradient(fpow, x8)[1], bigfd(fpow, x8), 0) < 1e-12
end
