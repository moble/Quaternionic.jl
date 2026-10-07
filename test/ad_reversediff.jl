# Tests of ReverseDiff on quaternions and of `QuaternionicReverseDiffExt`.  ReverseDiff
# differentiates quaternions of `TrackedReal`s natively, by recording the operations on
# their components.  The extension makes the source take the same value-dependent branches
# as it does for plain numbers, broadcasts a `TrackedReal` with an array of quaternions
# elementwise, and records the starting vector of the dominant eigenvector in
# `from_rotation_matrix` and `align` so that replayed tapes recompute it.  Every reference
# is a central finite difference evaluated with 256-bit `BigFloat`s (from `ADTestUtils`).

@testitem "ReverseDiff: value and iszerovalue see through tracking" tags=[:ad, :reversediff] begin
    using ReverseDiff, ForwardDiff
    using Test: detect_ambiguities
    import Quaternionic: value, iszerovalue

    t = ReverseDiff.track(0.3)
    @test value(t) === 0.3
    @test iszerovalue(ReverseDiff.track(0.0))
    @test !iszerovalue(t)
    # The value is also found through nested tracking, as in `ReverseDiff.hessian`, and
    # through a dual number inside the tracking.
    @test value(ReverseDiff.track(ReverseDiff.track(0.3))) === 0.3
    @test value(ReverseDiff.track(ForwardDiff.Dual(0.3, 1.0))) === 0.3
    @test iszerovalue(ReverseDiff.track(ForwardDiff.Dual(0.0, 1.0)))
    @test iszerovalue(quaternion(ReverseDiff.track(0.0), ReverseDiff.track(0.0), 0.0, 0.0))

    # The extension's methods are not ambiguous with ReverseDiff's own.
    ext = Base.get_extension(Quaternionic, :QuaternionicReverseDiffExt)
    @test ext isa Module
    ambiguities = detect_ambiguities(ext, ReverseDiff; recursive=true)
    @test isempty(filter(a -> a[1].module === ext || a[2].module === ext, ambiguities))
end

@testitem "ReverseDiff: gradients at the point set" tags=[:ad, :reversediff] setup=[ADTestUtils] begin
    using ReverseDiff
    using LinearAlgebra: normalize
    using .ADTestUtils

    # Each case maps a four-entry point `x` to a value, which `sc` turns into a real loss.
    # The kinds of points are those of `ADTestUtils.points`.
    cases = [
        ("Q * P0", x -> Q(x) * P0, :all),
        ("P0 * Q", x -> P0 * Q(x), :all),
        ("Q / P0", x -> Q(x) / P0, :all),
        ("P0 / Q", x -> P0 / Q(x), :all),
        ("Q + 2.5", x -> Q(x) + 2.5, :all),
        ("x[1] * P0", x -> x[1] * P0, :all),
        ("inv(Q)", x -> inv(Q(x)), :all),
        ("conj(Q)", x -> conj(Q(x)), :all),
        ("abs(Q)", x -> abs(Q(x)), :all),
        ("abs2(Q)", x -> abs2(Q(x)), :all),
        ("absvec(Q)", x -> absvec(Q(x)), :nonreal),
        ("normalize(Q)", x -> normalize(Q(x)), :all),
        ("exp(Q)", x -> exp(Q(x)), :all),
        ("exp(V)", x -> exp(V(x)), :all),
        ("log(Q)", x -> log(Q(x)), :nocut),
        ("sqrt(Q)", x -> sqrt(Q(x)), :nocut),
        ("Q^2", x -> Q(x)^2, :all),
        ("Q^-3", x -> Q(x)^-3, :all),
        ("Q^0.3", x -> Q(x)^0.3, :nocut),
        ("V^2", x -> V(x)^2, :all),
        ("R", x -> R(x), :all),
        ("R * R0", x -> R(x) * R0, :all),
        ("R / R0", x -> R(x) / R0, :all),
        ("Ru * R0", x -> Ru(x) * R0, :offsphere),
        ("log(R)", x -> log(R(x)), :nocut),
        ("sqrt(R)", x -> sqrt(R(x)), :nocut),
        ("R^0.3", x -> R(x)^0.3, :nocut),
        ("R(V0)", x -> R(x)(V0), :all),
        ("Ru(V0)", x -> Ru(x)(V0), :offsphere),
        ("slerp(R, R0, 0.3)", x -> slerp(R(x), R0, 0.3), :all),
        ("slerp(R0, R, 0.3)", x -> slerp(R0, R(x), 0.3), :all),
        ("distance2(R, R0)", x -> distance2(R(x), R0), :all),
        ("distance(R, R0)", x -> distance(R(x), R0), :all),
        ("distance2(R, R)", x -> distance2(R(x), R(x)), :all),
        ("Q ⋅ P0", x -> Q(x) ⋅ P0, :all),
        ("components(Q)", x -> components(Q(x)), :all),
        ("to_rotation_matrix(R)", x -> to_rotation_matrix(R(x)), :all),
        ("L", x -> L(x), :all),
        ("L * L", x -> L(x) * L(x), :all),
    ]
    # The Lorentz cases at `:identity` and `:minus1`, where the rotor's vector part is zero,
    # are regression tests: ReverseDiff records the derivative of `abs(::Complex)`, which is
    # `hypot`, as NaN at zero, so `_hypot` in src/quaternion.jl takes its scale from the
    # values of the components only.
    @testset "$name" for (name, f, kind) ∈ cases
        loss = x -> sc(f(x))
        for (pname, x) ∈ points(kind)
            g = ReverseDiff.gradient(loss, x)
            ref = bigfd(loss, x)
            @test all(isfinite, g)
            @test relerr(g, ref) < 1e-12
        end
    end
    # Base's `exp`, `cos`, and `sin` of a `Complex` take shortcuts when the real or imaginary
    # part is exactly zero, which ReverseDiff differentiates incorrectly.  The package
    # computes them from real functions instead, and these cases are regression tests for
    # that.  `x8` gives a real scalar part, `x8v` a vector part with a negative real square,
    # and `x8g` the generator of a boost; `Bp(x)` is a pure boost of rapidity `x[1]`.
    x8 = [0.7, 0.3, -0.2, 0.4, 0.0, 0.1, 0.05, -0.2]
    x8v = [0.7, 0.0, 0.0, 0.0, 0.2, 0.3, -0.2, 0.4]
    x8g = vcat(zeros(5), 0.6 .* NHAT)
    Bp(x) = Boost(x[1], NHAT)
    zerobranches = [
        ("exp(QC), real scalar part", x -> sc(exp(QC(x))), x8),
        ("exp(QC), imaginary absvec", x -> sc(exp(QC(x))), x8v),
        ("exp of a boost generator", x -> sc(exp(QC(x))), x8g),
        ("power of a boost", x -> sc(Bp(x)^x[2]), [1.0, 0.3]),
        ("large power of a small boost", x -> sc(Bp(x)^x[2]), [0.1, 30.0]),
    ]
    @testset "$name" for (name, loss, x) ∈ zerobranches
        @test relerr(ReverseDiff.gradient(loss, x), bigfd(loss, x)) < 1e-12
    end
end

@testitem "ReverseDiff: Hessians at special points" tags=[:ad, :reversediff] setup=[ADTestUtils] begin
    using ReverseDiff
    using .ADTestUtils

    # `ReverseDiff.hessian` nests `TrackedReal`s, so the value-based thresholds of the
    # source must see through both levels.  The slerp Hessian at q₁ == q₂ was wrong by 1.06
    # before `Quaternionic.value(::TrackedReal)` existed.
    x₁ = [0.4, -0.2, 0.7, 0.5]
    cases = [
        ("slerp at q₁ == q₂", x -> sc(slerp(R(x, 1), R(x, 5), x[9])), [x₁; x₁; 0.3]),
        ("slerp near q₁ == q₂", x -> sc(slerp(R(x, 1), R(x, 5), x[9])), [x₁; x₁ .+ 1e-9; 0.3]),
        ("R^0.3 at the identity", x -> sc(R(x)^0.3), point(:identity)),
        ("R^0.3 near the identity", x -> sc(R(x)^0.3), point(:nearidentity)),
        ("R^s at the identity", x -> sc(R(x)^x[5]), [point(:identity); 0.3]),
        ("distance2(R, R)", x -> distance2(R(x, 1), R(x, 5)), [x₁; x₁]),
        ("log(R) at the identity", x -> sc(log(R(x))), point(:identity)),
        ("exp(V) at zero", x -> sc(exp(V(x))), zeros(4)),
        ("absvec(Q) near the real axis", x -> absvec(Q(x)), point(:nearidentity)),
        ("sqrt(Q) at the identity", x -> sc(sqrt(Q(x))), point(:identity)),
    ]
    @testset "$name" for (name, f, x) ∈ cases
        H = ReverseDiff.hessian(f, x)
        @test all(isfinite, H)
        @test relerr(H, bigfdH(f, x)) < 1e-10
    end
end

@testitem "ReverseDiff: a TrackedReal with arrays of quaternions" tags=[:ad, :reversediff] setup=[ADTestUtils] begin
    using ReverseDiff
    using StaticArrays: SVector
    using .ADTestUtils

    # ReverseDiff's own `materialize` methods for a `TrackedReal` and an array of numbers
    # call `broadcast` again, which overflowed the stack for arrays of quaternions.
    qs = [quaternion(0.3, -0.5, 0.2, 0.9), quaternion(-1.1, 0.4, 0.7, -0.2), P0]
    Rs = [R0, rotor(0.1, 0.9, -0.3, 0.2)]
    M = [qs[1] qs[2]; qs[3] P0]
    sqs = SVector(qs[1], qs[2])
    cases = [
        ("t * qs", t -> t * qs),
        ("qs * t", t -> qs * t),
        ("qs / t", t -> qs / t),
        ("t .* qs", t -> t .* qs),
        ("qs .* t", t -> qs .* t),
        ("qs ./ t", t -> qs ./ t),
        ("t ./ qs", t -> t ./ qs),
        ("t .+ qs", t -> t .+ qs),
        ("qs .- t", t -> qs .- t),
        ("t .- qs", t -> t .- qs),
        ("qs .^ t", t -> qs .^ t),
        ("t .^ qs", t -> t .^ qs),
        ("Rs .* t", t -> Rs .* t),
        ("t .* M", t -> t .* M),
        ("view(qs) .* t", t -> view(qs, 1:2) .* t),
        ("t .* sqs", t -> t .* sqs),
        ("sqs ./ t", t -> sqs ./ t),
    ]
    @testset "$name" for (name, f) ∈ cases
        loss = x -> sc(f(x[1]))
        x = [0.7]
        g = ReverseDiff.gradient(loss, x)
        @test relerr(g, bigfd(loss, x)) < 1e-13
    end
    # The elementwise results have the same values as plain broadcasting, and the same
    # container types, so that static arrays give static results.
    t = ReverseDiff.track(0.7)
    @test ReverseDiff.value.(components(first(t .* qs))) == components(0.7 * qs[1])
    @test size(t .* M) == size(M)
    @test t .* M isa Matrix
    @test qs ./ t isa Vector
    @test t .* sqs isa SVector{2}
    @test sqs ./ t isa SVector{2}
end

@testitem "ReverseDiff: from_rotation_matrix and align" tags=[:ad, :reversediff] setup=[ADTestUtils] begin
    using ReverseDiff
    using LinearAlgebra: normalize
    import Quaternionic: value
    using .ADTestUtils: sc, bigfd, bigfdH, relerr

    # At an exact rotation, the matrix whose dominant eigenvector gives the rotor has a
    # triply degenerate eigenvalue.  The derivatives were NaN at the identity and at
    # pure-vector rotors, and wrong at second order everywhere, before version 4.4.5.  At
    # w = 0, the sign of the result is chosen by `positive_hemisphere`, so the loss is
    # multiplied by the sign that makes the function smooth for the finite differences.
    rotors = [
        rotor(1.0, 0.0, 0.0, 0.0),
        normalize(rotor(0.0, 0.48, -0.6, 0.64)),
        rotor(1.2, -0.7, 0.5, 0.3),
        rotor(1.0, 1e-9, -2e-9, 3e-9),
        rotor(0.0, 1.0, 0.0, 0.0),
        rotor(1e-8, 0.6, -0.48, 0.64),
    ]
    @testset "$R₀" for R₀ ∈ rotors
        f = x -> (S = from_rotation_matrix(reshape(x, 3, 3));
                  sign(value(real(S ⋅ R₀))) * sc(S))
        x₀ = vec(Matrix(to_rotation_matrix(R₀)))
        # The matrix `x₁` is near the rotation, but not orthogonal.
        x₁ = x₀ .+ 1e-3 .* sin.(1:9)
        for x ∈ (x₀, x₁)
            g = ReverseDiff.gradient(f, x)
            @test all(isfinite, g)
            @test relerr(g, bigfd(f, x)) < 1e-13
            H = ReverseDiff.hessian(f, x)
            @test all(isfinite, H)
            @test relerr(H, bigfdH(f, x)) < 1e-12
        end
    end

    # The function `align`, given vectors, uses the same dominant eigenvector.
    a⃗ = [quatvec(sin(k), cos(2k), sin(3k + 1)) for k ∈ 1:5]
    @testset "align at $R₀" for R₀ ∈ (rotor(0.4, -0.2, 0.7, 0.5), rotor(1.0, 0.0, 0.0, 0.0))
        f = x -> sc(align(a⃗, [quatvec(x[3k-2], x[3k-1], x[3k]) for k ∈ 1:5]))
        x₀ = reduce(vcat, [collect(vec(R₀(a))) for a ∈ a⃗])
        for x ∈ (x₀, x₀ .+ 1e-2 .* cos.(1:15))
            g = ReverseDiff.gradient(f, x)
            @test all(isfinite, g)
            @test relerr(g, bigfd(f, x)) < 1e-13
        end
    end
end

@testitem "ReverseDiff: recorded and compiled tapes" tags=[:ad, :reversediff] setup=[ADTestUtils] begin
    using ReverseDiff
    using .ADTestUtils

    # A tape records the branches taken at the recording point, so it is replayed here at
    # generic points, on the same side of every threshold as the recording point.
    x₀ = [1.2, -0.7, 0.5, 0.3]
    replays = ([1.1, -0.6, 0.55, 0.25], [0.9, 0.1, -0.3, 0.6], x₀)
    cases = [
        ("exp(Q)", x -> sc(exp(Q(x)))),
        ("log(Q)", x -> sc(log(Q(x)))),
        ("sqrt(Q)", x -> sc(sqrt(Q(x)))),
        ("Q^0.3", x -> sc(Q(x)^0.3)),
        ("R^0.3", x -> sc(R(x)^0.3)),
        ("slerp(R, R0, 0.3)", x -> sc(slerp(R(x), R0, 0.3))),
        ("distance2(R, R0)", x -> distance2(R(x), R0)),
        ("R(V0)", x -> sc(R(x)(V0))),
        ("from_rotation_matrix", x -> sc(from_rotation_matrix(to_rotation_matrix(Q(x))))),
        ("align", x -> sc(align([V0, quatvec(0.5, 0.1, -0.4), quatvec(-0.2, 0.3, 0.9)],
                                [R(x)(V0), R(x)(quatvec(0.5, 0.1, -0.4)), R(x)(quatvec(-0.2, 0.3, 0.9))]))),
    ]
    @testset "$name" for (name, f) ∈ cases
        tape = ReverseDiff.GradientTape(f, x₀)
        compiled = ReverseDiff.compile(tape)
        for x ∈ replays
            ref = bigfd(f, x)
            @test relerr(ReverseDiff.gradient!(similar(x), tape, x), ref) < 1e-13
            @test relerr(ReverseDiff.gradient!(similar(x), compiled, x), ref) < 1e-13
        end
    end

    # The sign chosen by `positive_hemisphere` is frozen on the tape.  A replay at a rotation
    # in the other hemisphere of the recording (here, `real(R(x₀) ⋅ R(x)) < 0`) therefore
    # returns the derivative of the negated rotor, not an unrelated value.
    f = x -> sc(from_rotation_matrix(to_rotation_matrix(Q(x))))
    compiled = ReverseDiff.compile(ReverseDiff.GradientTape(f, x₀))
    x = [0.3, 0.8, -0.2, 0.5]
    @test real(R(x₀) ⋅ R(x)) < 0
    @test relerr(ReverseDiff.gradient!(similar(x), compiled, x), -bigfd(f, x)) < 1e-13

    # A `HessianTape` records the inner gradient tape onto an outer tape, so the starting
    # vector of the dominant eigenvector must be recomputed on the outer tape as well.
    # Before it was, replays away from the recording point were wrong by up to 37.  The
    # input matrices are slightly non-orthogonal, and every rotation is in the hemisphere
    # of the recording point.
    a⃗ = [quatvec(sin(k), cos(2k), sin(3k + 1)) for k ∈ 1:5]
    hessiancases = [
        (
            "from_rotation_matrix",
            x -> sc(from_rotation_matrix(reshape(x, 3, 3))),
            R -> vec(Matrix(to_rotation_matrix(R))) .+ 1e-3 .* sin.(1:9),
        ),
        (
            "align",
            x -> sc(align(a⃗, [quatvec(x[3k-2], x[3k-1], x[3k]) for k ∈ 1:5])),
            R -> reduce(vcat, [collect(vec(R(a))) for a ∈ a⃗]) .+ 1e-2 .* cos.(1:15),
        ),
    ]
    @testset "HessianTape of $name" for (name, f, input) ∈ hessiancases
        tape = ReverseDiff.HessianTape(f, input(rotor(1.2, -0.7, 0.5, 0.3)))
        compiled = ReverseDiff.compile(tape)
        for Rₓ ∈ (rotor(1.2, -0.7, 0.5, 0.3), rotor(0.9, 0.4, -0.3, 0.2), rotor(1, 0, 0, 0))
            x = input(Rₓ)
            ref = bigfdH(f, x)
            @test relerr(ReverseDiff.hessian!(tape, x), ref) < 1e-12
            @test relerr(ReverseDiff.hessian!(compiled, x), ref) < 1e-12
        end
    end
end
