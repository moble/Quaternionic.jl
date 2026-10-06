# Tests of Enzyme, which differentiates the package's source natively, and of
# `QuaternionicEnzymeExt`, which adds `Enzyme.onehot` for quaternion arguments, rules for
# the LAPACK call behind `from_rotation_matrix` and `align`, and rules for `abs` and
# `absvec` that work around a bug in Enzyme's reverse rule for `hypot`.  The references are
# central finite differences in 256-bit `BigFloat` arithmetic (`bigfd` from `ADTestUtils`),
# except for the Hessians, which are compared with ForwardDiff's.
#
# Enzyme's reverse mode is tested at batch width 1 (the gradient of a scalar loss) except in
# the items on batched reverse mode and on `dominant_eigenvector_lapack`.  At batch widths
# greater than 1, which DifferentiationInterface uses for Jacobians, Enzyme 0.13.205 crashes
# the compiler on `squad` and `unflip` (an upstream bug), so those functions are tested only
# at batch width 1.  Enzyme's reverse rule for `hypot` with three or more arguments also
# fails at those batch widths (another upstream bug).  The extension's rules for `abs` and
# `absvec` avoid it, and `sqrt` of a quaternion with a negative scalar part computes the
# norm of its vector part with `abs` for that reason.  Forward mode pushes every tangent
# of a vector argument at once (batch width 4 or more), which works.

@testmodule EnzymeHelpers begin
    using Quaternionic
    import Enzyme

    export enzyme_gradients, enzyme_hessian

    """
        enzyme_gradients(g, x)

    Return the gradients of the real function `g` at the vector `x` computed by Enzyme in
    reverse mode (batch width 1) and in forward mode, as two `Vector`s.
    """
    function enzyme_gradients(g, x)
        reverse = Enzyme.gradient(Enzyme.Reverse, g, x)[1]
        forward = Enzyme.gradient(Enzyme.Forward, g, x)[1]
        collect(reverse), collect(forward)
    end

    """
        enzyme_hessian(g, x; runtime_activity=false)

    Return the Hessian of the real function `g` at `x`, computed by Enzyme in forward mode
    over reverse mode.
    """
    function enzyme_hessian(g, x; runtime_activity=false)
        fwd = runtime_activity ? Enzyme.set_runtime_activity(Enzyme.Forward) : Enzyme.Forward
        rev = runtime_activity ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.Reverse
        Enzyme.jacobian(fwd, y -> Enzyme.gradient(rev, g, y)[1], x)[1]
    end
end


@testitem "Enzyme: onehot and quaternion arguments" tags=[:ad, :enzyme] setup=[ADTestUtils] begin
    import Enzyme
    using .ADTestUtils

    # A loss that involves every component, the product with a constant, and `absvec`, which
    # has a rule in the extension.
    f(q) = sc(q * P0) + absvec(q)
    x = point(:generic)
    for (name, make) in (("Quaternion", Q), ("Rotor", Ru), ("QuatVec", V))
        q = make(x)
        ref = bigfd(y -> f(make(y)), x)
        # `V` ignores `x[1]`, and the tangents of a `QuatVec` are its three vector components.
        expected = q isa QuatVec ? ref[2:4] : ref
        n = length(expected)

        # The tangents are unit quaternions of the argument's own type.
        tangents = Enzyme.onehot(q)
        @test length(tangents) == n
        @test all(t -> typeof(t) === typeof(q), tangents)
        @test [collect(q isa QuatVec ? vec(t) : components(t)) for t in tangents] ==
              [Float64.(1:n .== i) for i in 1:n]

        # Forward mode, unchunked and in chunks of 1 and 3.
        for chunk in (nothing, Val(1), Val(3))
            g = Enzyme.gradient(Enzyme.Forward, f, q; chunk=chunk)[1]
            @test g isa Tuple
            @test length(g) == n
            @test relerr(collect(g), expected) < 1e-14
        end

        # Reverse mode returns a shadow of the argument's own type, holding the raw gradient
        # components; for a `Rotor` it is not normalized.
        gr = Enzyme.gradient(Enzyme.Reverse, f, q)[1]
        @test typeof(gr) === typeof(q)
        @test relerr(collect(q isa QuatVec ? vec(gr) : components(gr)), expected) < 1e-14
    end
end


@testitem "Enzyme: quaternion algebra" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    using .ADTestUtils, .EnzymeHelpers
    t = 2.5
    cases = [
        ("Q * P0", x -> Q(x) * P0, :all),
        ("P0 * Q", x -> P0 * Q(x), :all),
        ("Q / P0", x -> Q(x) / P0, :all),
        ("P0 / Q", x -> P0 / Q(x), :all),
        ("Q \\ P0", x -> Q(x) \ P0, :all),
        ("inv(Q)", x -> inv(Q(x)), :all),
        ("conj(Q)", x -> conj(Q(x)), :all),
        ("-Q", x -> -Q(x), :all),
        ("Q + P0", x -> Q(x) + P0, :all),
        ("P0 - Q", x -> P0 - Q(x), :all),
        ("t * Q", x -> t * Q(x), :all),
        ("Q / t", x -> Q(x) / t, :all),
        ("x[1] * Q", x -> x[1] * Q(x), :all),
        ("t / Q", x -> t / Q(x), :all),
        ("abs(Q)", x -> abs(Q(x)), :all),
        ("abs2(Q)", x -> abs2(Q(x)), :all),
        ("absvec(Q)", x -> absvec(Q(x)), :nonreal),
        ("abs2vec(Q)", x -> abs2vec(Q(x)), :all),
        ("normalize(Q)", x -> normalize(Q(x)), :all),
        ("Q ⋅ P0", x -> Q(x) ⋅ P0, :all),
        ("V * V0", x -> V(x) * V0, :nonreal),
        ("V × V0", x -> V(x) × V0, :nonreal),
        ("abs(V)", x -> abs(V(x)), :nonreal),
        ("V^2", x -> V(x)^2, :all),
    ]
    for (name, f, kind) in cases, (pointname, x) in points(kind)
        g = y -> sc(f(y))
        ref = bigfd(g, x)
        reverse, forward = enzyme_gradients(g, x)
        @testset "$name at $pointname" begin
            @test relerr(reverse, ref) < 1e-13
            @test relerr(forward, ref) < 1e-13
        end
    end

    # At a zero result, the rules for `abs` and `absvec` give a zero gradient, as Enzyme's
    # own rules for `hypot` do (the derivative does not exist there).
    for (name, g, x) in (("abs(Q) at zero", x -> abs(Q(x)), zeros(4)),
                         ("absvec(Q) at the identity", x -> absvec(Q(x)), [1.0, 0.0, 0.0, 0.0]),
                         ("abs(V) at zero", x -> abs(V(x)), [0.7, 0.0, 0.0, 0.0]))
        reverse, forward = enzyme_gradients(g, x)
        @testset "$name" begin
            @test reverse == zeros(4)
            @test forward == zeros(4)
        end
    end
end


@testitem "Enzyme: n-ary operations, constants, structure, and broadcasting" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    using .ADTestUtils, .EnzymeHelpers
    cases = [
        ("muladd", x -> muladd(Q(x), P0, x[2]), :all),
        ("t*t*t*t*q", x -> x[1] * x[2] * x[3] * x[4] * P0, :all),
        ("t*t*t*q*p", x -> x[1] * x[2] * x[3] * Q(x) * P0, :all),
        ("q*p*q*p", x -> Q(x) * P0 * Q(x) * P0, :all),
        ("t+2t+3t+q", x -> x[1] + 2x[2] + 3x[3] + Q(x), :all),
        ("t+q+p+t", x -> (x[1] + Q(x) + P0 + x[4]) * P0, :all),
        ("fastmath", x -> @fastmath(x[1] * Q(x) + x[2]), :all),
        ("fastmath division", x -> @fastmath(Q(x) / P0 - x[2] / Q(x)), :all),
        ("exp(θ𝐣)", x -> exp(x[1] * 𝐣 + x[2] * 𝐤), :all),
        ("θ+𝐣", x -> (x[1] + 𝐣) * Q(x), :all),
        ("θ-2𝐤", x -> x[1] - 2𝐤, :all),
        ("q*𝐢", x -> Q(x) * 𝐢 + 𝐤 * Q(x), :all),
        ("imag", x -> imag(Q(x) * P0), :all),
        ("vec", x -> vec(Q(x) * P0), :all),
        ("components", x -> components(R(x) * P0), :all),
        ("getindex", x -> (Q(x) * P0)[2] - 3 * (Q(x) * P0)[4], :all),
        ("to_float_array(q)", x -> to_float_array(Q(x) * P0), :all),
        ("to_float_array(qs)", x -> to_float_array([Q(x), P0 * Q(x)]), :all),
        ("quaternion(vector)", x -> quaternion(x) * P0, :all),
        ("rotor(vector)", x -> rotor(x) * P0, :all),
        ("qs .* ts", x -> [Q(x), P0] .* [x[2], 2.0], :all),
        ("t .* qs", x -> x[1] .* [Q(x), P0], :all),
        ("t * qs", x -> x[1] * [Q(x), P0], :all),
        ("qs / t", x -> [Q(x), P0] / x[1], :nonzeroscalar),
        ("qs ./ q", x -> [P0, Q(x)] ./ Q(x), :all),
        ("qs .^ 2", x -> [Q(x), P0] .^ 2, :all),
        ("qs .^ -1", x -> [Q(x), P0] .^ -1, :all),
        ("abs2.(qs)", x -> sum(abs2.([Q(x), P0])), :all),
        ("conj.(qs)", x -> conj.([Q(x), P0 * Q(x)]), :all),
        ("sum", x -> sum([Q(x), P0 * Q(x)]), :all),
        ("prod", x -> prod([Q(x), P0, Q(x)]), :all),
        ("map", x -> map(q -> q * P0, [Q(x), R(x)]), :all),
    ]
    for (name, f, kind) in cases
        pts = kind === :nonzeroscalar ? points(exclude=(:purevector,)) : points(kind)
        for (pointname, x) in pts
            g = y -> sc(f(y))
            ref = bigfd(g, x)
            reverse, forward = enzyme_gradients(g, x)
            @testset "$name at $pointname" begin
                @test relerr(reverse, ref) < 1e-13
                @test relerr(forward, ref) < 1e-13
            end
        end
    end
end


@testitem "Enzyme: elementary functions and powers" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    using .ADTestUtils, .EnzymeHelpers
    cases = [
        ("exp(Q)", x -> exp(Q(x)), :all),
        ("log(Q)", x -> log(Q(x)), :nocut),
        ("sqrt(Q)", x -> sqrt(Q(x)), :nocut),
        ("Q^0.3", x -> Q(x)^0.3, :nocut),
        ("Q^P0", x -> Q(x)^P0, :nocut),
        ("Q^t", x -> Q(x)^(x[2] + 0.5), :nocut),
        ("2^Q", x -> 2.0^Q(x), :all),
        ("R0^Q", x -> R0^Q(x), :all),
        ("exp(V)", x -> exp(V(x)), :all),
        ("sqrt(V)", x -> sqrt(V(x)), :nonreal),
        [("Q^$n", x -> Q(x)^n, :all) for n in -3:3]...,
    ]
    for (name, f, kind) in cases, (pointname, x) in points(kind)
        g = y -> sc(f(y))
        ref = bigfd(g, x)
        reverse, forward = enzyme_gradients(g, x)
        # Near -1, the derivatives of `log`, `sqrt`, and non-integer powers are large, and
        # the error grows with them.
        tol = pointname === :nearminus1 ? 1e-9 : 1e-13
        @testset "$name at $pointname" begin
            @test relerr(reverse, ref) < tol
            @test relerr(forward, ref) < tol
        end
    end
end


@testitem "Enzyme: rotors" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    using .ADTestUtils, .EnzymeHelpers
    cases = [
        ("R", x -> R(x), :all),
        ("R * R0", x -> R(x) * R0, :all),
        ("R0 / R", x -> R0 / R(x), :all),
        ("R / R0", x -> R(x) / R0, :all),
        ("conj(R)", x -> conj(R(x)), :all),
        ("R(V0)", x -> R(x)(V0), :all),
        ("log(R)", x -> log(R(x)), :nocut),
        ("sqrt(R)", x -> sqrt(R(x)), :nocut),
        ("R^0.3", x -> R(x)^0.3, :nocut),
        ("R^-2", x -> R(x)^-2, :all),
        ("angle(R)", x -> angle(R(x)), :nonreal),
        ("distance2(R, R0)", x -> distance2(R(x), R0), :all),
        ("distance(R, R0)", x -> distance(R(x), R0), :all),
        ("slerp(R0, R, 0.3)", x -> slerp(R0, R(x), 0.3), :nocut),
        ("slerp(R, R0, 0.3)", x -> slerp(R(x), R0, 0.3), :nocut),
        ("slerp(R0, R, τ)", x -> slerp(R0, R(x), x[1]), :nocut),
        ("slerp(R, R, 0.3)", x -> slerp(R(x), R(x), 0.3), :all),
        ("R^3", x -> R(x)^3, :all),
        ("distance2(Q, P0)", x -> distance2(Q(x), P0), :all),
        ("R0(V)", x -> R0(V(x)), :all),
        ("R * V0 * conj(R)", x -> R(x) * V0 * conj(R(x)), :all),
        ("to_rotation_matrix(R)", x -> to_rotation_matrix(R(x)), :all),
        # Rotors off the unit sphere, whose components are stored without normalization
        ("Ru * R0", x -> Ru(x) * R0, :offsphere),
        ("Ru / R0", x -> Ru(x) / R0, :offsphere),
        ("log(Ru)", x -> log(Ru(x)), :offsphere),
        ("Ru^0.3", x -> Ru(x)^0.3, :offsphere),
        ("sqrt(Ru)", x -> sqrt(Ru(x)), :offsphere),
        ("Ru^2", x -> Ru(x)^2, :offsphere),
        ("Ru^-1", x -> Ru(x)^-1, :offsphere),
        ("Ru * Ru", x -> Ru(x) * Ru(x), :offsphere),
        ("V0 / Ru", x -> V0 / Ru(x), :offsphere),
        ("absvec(Ru)", x -> absvec(Ru(x)), :offsphere),
    ]
    for (name, f, kind) in cases, (pointname, x) in points(kind)
        g = y -> sc(f(y))
        ref = bigfd(g, x)
        reverse, forward = enzyme_gradients(g, x)
        tol = pointname === :nearminus1 ? 1e-9 : 1e-13
        @testset "$name at $pointname" begin
            @test relerr(reverse, ref) < tol
            @test relerr(forward, ref) < tol
        end
    end
end


@testitem "Enzyme: Lorentz rotors and complex components" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    using .ADTestUtils, .EnzymeHelpers
    L0 = L([0.4, 0.2, -0.3, 0.5])
    B = Boost(0.2, NHAT)
    # `C(x)` is a quaternion with complex components built from the four entries of a point.
    C(x) = QC(vcat(x, [x[2] / 3, -x[3] / 2, x[4] / 5 + 0.2, x[1] / 4 - 0.1]))
    # At `:minus1`, both `L(x)` and `C(x)` have negative real eigenvalues (as complex 2×2
    # matrices), so they lie on the branch cut of `log` and `sqrt`.
    nocut = (:minus1,)
    cases = [
        ("L * L0", x -> L(x) * L0, ()),
        ("L0 * L", x -> L0 * L(x), ()),
        ("L * R0", x -> L(x) * R0, ()),
        ("L / B", x -> L(x) / B, ()),
        ("L * L", x -> L(x) * L(x), ()),
        ("L^2 + L^-1 + L^3", x -> L(x)^2 + L(x)^-1 + L(x)^3, ()),
        ("conj(L)", x -> conj(L(x)), ()),
        ("inv(L)", x -> inv(L(x)), ()),
        ("L(V0)", x -> L(x)(V0), ()),
        ("L(vector)", x -> sum(L(x)([1.0, 0.2, -0.3, 0.5])), ()),
        ("log(L)", x -> log(L(x)), nocut),
        ("sqrt(L)", x -> sqrt(L(x)), nocut),
        ("abs(quaternion(L))", x -> abs(quaternion(L(x))), ()),
        ("C * P0 + P0 * C", x -> C(x) * P0 + P0 * C(x), ()),
        ("C / P0 + P0 / C", x -> C(x) / P0 + P0 / C(x), ()),
        ("exp(C)", x -> exp(C(x)), ()),
        ("log(C)", x -> log(C(x)), nocut),
        ("sqrt(C)", x -> sqrt(C(x)), nocut),
        ("C^2", x -> C(x)^2, ()),
        ("abs2(C)", x -> abs2(C(x)), ()),
        ("(t + 0.5im) * C", x -> (x[2] + 0.5im) * C(x), ()),
    ]
    for (name, f, exclude) in cases, (pointname, x) in points(:all; exclude=exclude)
        g = y -> sc(f(y))
        ref = bigfd(g, x)
        reverse, forward = enzyme_gradients(g, x)
        tol = pointname === :nearminus1 ? 1e-9 : 1e-13
        @testset "$name at $pointname" begin
            @test relerr(reverse, ref) < tol
            @test relerr(forward, ref) < tol
        end
    end
    # At `:identity` and `:minus1`, the vector part of `C(x)` above is imaginary, so `exp`
    # computes `cos` and `sin` of an imaginary number.  Base's `sin`, `cos`, and `exp` of a
    # `Complex` take shortcuts when the real or imaginary part is exactly zero, and Enzyme
    # differentiates those shortcuts, giving wrong derivatives (scratchpad reproduction
    # `upstream/enzyme_5.jl`).  The package therefore computes them from real functions
    # (`Quaternionic.cexp`, `ccos`, and `csin`), and these cases are regression tests for
    # that.  `x8` gives a real scalar part, `x8v` a vector part with a negative real square,
    # and `x8g` the generator of a boost; `Bp(x)` is a pure boost of rapidity `x[1]`, whose
    # non-integer powers reach the same functions.
    x8 = [0.7, 0.3, -0.2, 0.4, 0.0, 0.1, 0.05, -0.2]
    x8v = [0.7, 0.0, 0.0, 0.0, 0.2, 0.3, -0.2, 0.4]
    x8g = vcat(zeros(5), 0.6 .* NHAT)
    VC(x) = quatvec(complex(x[2], x[6]), complex(x[3], x[7]), complex(x[4], x[8]))
    Bp(x) = Boost(x[1], NHAT)
    zerobranches = [
        ("exp(QC), real scalar part", x -> sc(exp(QC(x))), x8),
        ("exp(QC), imaginary absvec", x -> sc(exp(QC(x))), x8v),
        ("exp of a boost generator", x -> sc(exp(QC(x))), x8g),
        ("exp of a boost generator as a QuatVec", x -> sc(exp(VC(x))), x8g),
        ("power of a boost", x -> sc(Bp(x)^x[2]), [1.0, 0.3]),
        ("large power of a small boost", x -> sc(Bp(x)^x[2]), [0.1, 30.0]),
        ("L^s, pure boost", x -> sc(L(x)^0.3) + sc(L(x)^2.7), point(:identity)),
    ]
    @testset "$name" for (name, g, x) in zerobranches
        ref = bigfd(g, x)
        reverse, forward = enzyme_gradients(g, x)
        @test relerr(reverse, ref) < 1e-13
        @test relerr(forward, ref) < 1e-13
    end
end


@testitem "Enzyme: batched reverse mode" tags=[:ad, :enzyme] setup=[ADTestUtils] begin
    # Jacobians in reverse mode pull back several cotangents at once.  Enzyme's own reverse
    # rule for `hypot` with three or more arguments fails at batch widths greater than 1, so
    # without the extension's rules for `abs` and `absvec`, every function that builds a
    # `Rotor` failed here.
    import Enzyme
    import DifferentiationInterface as DI
    using .ADTestUtils
    backend = DI.AutoEnzyme(mode=Enzyme.Reverse, function_annotation=Enzyme.Const)
    cases = [
        ("R", x -> R(x), :all),
        ("R * R0", x -> R(x) * R0, :all),
        ("normalize(Q)", x -> normalize(Q(x)), :all),
        ("abs(Q)", x -> abs(Q(x)), :all),
        ("abs(Q), absvec(Q)", x -> (abs(Q(x)), absvec(Q(x))), :nonreal),
        ("abs(V)", x -> abs(V(x)), :nonreal),
        ("absvec(Ru)", x -> absvec(Ru(x)), :offsphere),
        ("log(R)", x -> log(R(x)), :nocut),
        ("sqrt(Q)", x -> sqrt(Q(x)), :nocut),
        ("sqrt(R)", x -> sqrt(R(x)), :nocut),
        ("R^0.3", x -> R(x)^0.3, :nocut),
        ("slerp(R0, R, 0.3)", x -> slerp(R0, R(x), 0.3), :nocut),
        ("L * R0", x -> L(x) * R0, :all),
    ]
    for (name, f, kind) in cases, (pointname, x) in points(kind)
        h = y -> realcoords(f(y))
        ref = bigfd(h, x)
        tol = pointname === :nearminus1 ? 1e-9 : 1e-13
        @testset "$name at $pointname" begin
            @test relerr(DI.jacobian(h, backend, x), ref) < tol
        end
    end

    # The same rule, called directly with two cotangents of a scalar result
    x = point(:generic)
    dx = (zero(x), zero(x))
    Enzyme.autodiff(Enzyme.Reverse, y -> abs(Q(y)), Enzyme.Active, Enzyme.BatchDuplicated(x, dx))
    @test relerr(dx[1], x / abs(Q(x))) < 1e-15
    @test relerr(dx[2], x / abs(Q(x))) < 1e-15

    # At a zero result, the gradients are zero, as at batch width 1.
    for (g, x) in ((y -> abs(Q(y)), zeros(4)), (y -> absvec(Q(y)), [1.0, 0.0, 0.0, 0.0]))
        dx = (zeros(4), zeros(4))
        Enzyme.autodiff(Enzyme.Reverse, g, Enzyme.Active, Enzyme.BatchDuplicated(x, dx))
        @test dx[1] == zeros(4)
        @test dx[2] == zeros(4)
        @test DI.jacobian(y -> [g(y)], backend, x) == zeros(1, 4)
    end
end


@testitem "Enzyme: from_rotation_matrix and align" tags=[:ad, :enzyme] setup=[ADTestUtils] begin
    import Enzyme
    import DifferentiationInterface as DI
    using LinearAlgebra: I
    using .ADTestUtils

    # Rotation matrices of rotors, at a generic point, at the identity, where the
    # Bar-Itzhack matrix has a triply degenerate eigenvalue, and at a pure-vector rotor.
    # There, `from_rotation_matrix` chooses the sign of its result by a rule that changes
    # under perturbation, so the loss is multiplied by the sign that makes the result agree
    # with the input rotor; the derivative of that sign is zero.
    function g1(x)
        Rx = R(x)
        Rf = from_rotation_matrix(to_rotation_matrix(Rx))
        sc(Rf) * sign(Rf ⋅ Rx)
    end
    # Non-orthogonal matrices
    g2(x) = sc(from_rotation_matrix(reshape(x[1:9], 3, 3) + I))
    # `align` of rotated vectors, which needs runtime activity in Enzyme, because its closures
    # mix active and constant arrays
    vs = [quatvec(1.0, 0.2, -0.3), quatvec(-0.4, 0.9, 0.1), quatvec(0.3, 0.3, 1.0)]
    g3(x) = sc(align([R(x)(v) for v in vs], vs))
    g4(x) = sc(align([R(x)(v) for v in vs], vs, [1.0, 0.5, 2.0]))
    for (name, g, x, runtime_activity) in (
        ("from_rotation_matrix at generic", g1, point(:generic), false),
        ("from_rotation_matrix at identity", g1, point(:identity), false),
        ("from_rotation_matrix at purevector", g1, point(:purevector), false),
        ("from_rotation_matrix of a non-orthogonal matrix", g2,
         0.1 .* collect(1.0:9.0) .- 0.4, false),
        ("align(::Vector{QuatVec})", g3, point(:generic), true),
        ("align(::Vector{QuatVec}, w)", g4, point(:generic), true),
    )
        ref = bigfd(g, x)
        rev = runtime_activity ? Enzyme.set_runtime_activity(Enzyme.Reverse) : Enzyme.Reverse
        fwd = runtime_activity ? Enzyme.set_runtime_activity(Enzyme.Forward) : Enzyme.Forward
        @testset "$name" begin
            @test relerr(Enzyme.gradient(rev, g, x)[1], ref) < 1e-13
            @test relerr(collect(Enzyme.gradient(fwd, g, x)[1]), ref) < 1e-13
        end
    end

    # Single precision, through the `Float32` methods of the eigenvector rules
    x32 = Float32.(point(:generic))
    ref32 = bigfd(g1, Float64.(x32))
    g32 = Enzyme.gradient(Enzyme.Reverse, g1, x32)[1]
    @test eltype(g32) === Float32
    @test relerr(g32, ref32) < 1e-5
    @test relerr(collect(Enzyme.gradient(Enzyme.Forward, g1, x32)[1]), ref32) < 1e-5

    # Jacobians, which push forward or pull back several tangents at once
    h(x) = realcoords(from_rotation_matrix(to_rotation_matrix(R(x))))
    x = point(:generic)
    ref = bigfd(h, x)
    for mode in (Enzyme.Forward, Enzyme.Reverse)
        backend = DI.AutoEnzyme(mode=mode, function_annotation=Enzyme.Const)
        @test relerr(DI.jacobian(h, backend, x), ref) < 1e-13
    end
end


@testitem "Enzyme: rules for dominant_eigenvector_lapack" tags=[:ad, :enzyme] setup=[ADTestUtils] begin
    import Enzyme
    using LinearAlgebra: Symmetric, eigen, dot, triu, tril
    using .ADTestUtils: relerr

    # The loss of the test with runtime activity below, whose matrix argument is active only
    # for `t < 0`
    function h(t, B, D, w, uplo)
        C = t[1] > 0 ? B : B + t[1] * D
        dot(w, Quaternionic.dominant_eigenvector_lapack(C, uplo)) * t[1]
    end

    # A symmetric matrix with a simple dominant eigenvalue, stored in both triangles, with
    # different (ignored) entries in the other triangle
    S = [4.0 0.3 -0.2 0.1; 0.3 1.0 0.5 -0.4; -0.2 0.5 0.7 0.2; 0.1 -0.4 0.2 -1.0]
    w = [0.3, -1.1, 0.7, 0.45]
    for uplo in ('U', 'L')
        A = uplo == 'U' ? triu(S) + tril(fill(9.0, 4, 4), -1) :
                          tril(S) + triu(fill(9.0, 4, 4), 1)
        v₀ = Quaternionic.dominant_eigenvector_lapack(A, uplo)
        f = B -> dot(w, Quaternionic.dominant_eigenvector_lapack(B, uplo))

        # The reference uses the generic eigensolver in `BigFloat`, with the eigenvector's
        # sign chosen to agree with `vref` (by default `v₀`).
        fbig = (a, vref=v₀) -> begin
            E = eigen(Symmetric(reshape(a, 4, 4), Symbol(uplo)))
            v = E.vectors[:, argmax(E.values)]
            dot(w, dot(v, vref) < 0 ? -v : v)
        end
        ref = reshape(ADTestUtils.bigfd(fbig, vec(A)), 4, 4)
        other = uplo == 'U' ? tril(trues(4, 4), -1) : triu(trues(4, 4), 1)
        @test all(iszero, ref[other])

        Ā = zero(A)
        Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, Enzyme.Duplicated(A, Ā))
        @test relerr(Ā, ref) < 1e-13
        @test all(iszero, Ā[other])

        # Forward mode along one tangent, and along two at once
        Ȧ = reshape(collect(1.0:16.0), 4, 4) ./ 10
        ḟ = Enzyme.autodiff(Enzyme.Forward, f, Enzyme.Duplicated(A, Ȧ))[1]
        @test ḟ ≈ dot(ref, Ȧ) rtol=1e-13
        ḟs = Enzyme.autodiff(Enzyme.Forward, f, Enzyme.BatchDuplicated(A, (Ȧ, 2Ȧ)))[1]
        @test ḟs[1] ≈ dot(ref, Ȧ) rtol=1e-13
        @test ḟs[2] ≈ 2dot(ref, Ȧ) rtol=1e-13

        # Reverse mode with two cotangents at once
        Ās = (zero(A), zero(A))
        Enzyme.autodiff(Enzyme.Reverse, f, Enzyme.Active, Enzyme.BatchDuplicated(A, Ās))
        @test relerr(Ās[1], ref) < 1e-13
        @test relerr(Ās[2], ref) < 1e-13

        # With runtime activity, Enzyme passes a matrix that is inactive at run time as a
        # `Duplicated` whose shadow is the matrix itself.  The rules must treat it as
        # constant, and in particular must not accumulate a cotangent into it.  Here the
        # matrix is inactive for `t > 0` and active (through `B + t * D`) for `t < 0`.
        D = [0.0 0.2 0.1 -0.3; 0.2 0.5 -0.1 0.0; 0.1 -0.1 0.0 0.4; -0.3 0.0 0.4 -0.2]
        A₀ = copy(A)
        for t in (0.5, -0.5)
            # The sign of the eigenvector is that which LAPACK chooses at `t`.
            vt = Quaternionic.dominant_eigenvector_lapack(t > 0 ? A : A + t * D, uplo)
            href = ADTestUtils.bigfd(s -> (s > 0 ? fbig(vec(big.(A)), vt) : fbig(vec(A + s * D), vt)) * s, t)
            dt = zeros(1)
            Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), h, Enzyme.Active,
                            Enzyme.Duplicated([t], dt), Enzyme.Const(A), Enzyme.Const(D),
                            Enzyme.Const(w), Enzyme.Const(uplo))
            @test dt[1] ≈ href rtol=1e-13
            @test A == A₀
            dts = (zeros(1), zeros(1))
            Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Reverse), h, Enzyme.Active,
                            Enzyme.BatchDuplicated([t], dts), Enzyme.Const(A), Enzyme.Const(D),
                            Enzyme.Const(w), Enzyme.Const(uplo))
            @test dts[1][1] ≈ href rtol=1e-13
            @test dts[2][1] ≈ href rtol=1e-13
            @test A == A₀
            ḣ = Enzyme.autodiff(Enzyme.set_runtime_activity(Enzyme.Forward), h,
                                Enzyme.Duplicated([t], [1.0]), Enzyme.Const(A), Enzyme.Const(D),
                                Enzyme.Const(w), Enzyme.Const(uplo))[1]
            @test ḣ ≈ href rtol=1e-13
            @test A == A₀
        end
    end
end


@testitem "Enzyme: forward-over-reverse Hessians" tags=[:ad, :enzyme] setup=[ADTestUtils, EnzymeHelpers] begin
    import ForwardDiff
    using .ADTestUtils, .EnzymeHelpers
    cases = [
        ("slerp(R0, R, 0.3)", x -> sc(slerp(R0, R(x), 0.3)), :nocut),
        ("slerp(R, R, 0.3)", x -> sc(slerp(R(x), R(x), 0.3)), :all),
        ("R^0.3", x -> sc(R(x)^0.3), :nocut),
        ("exp(Q)", x -> sc(exp(Q(x))), :all),
        ("log(Q)", x -> sc(log(Q(x))), :nocut),
        ("distance2(R, R0)", x -> distance2(R(x), R0), :all),
        ("distance2(R, R * R1)", x -> distance2(R(x), R(x) * rotor(1.0, 1e-3, 0.0, 0.0)), :all),
        ("abs(Q)", x -> abs(Q(x)), :all),
        ("absvec(Q)", x -> absvec(Q(x)), :nonreal),
    ]
    for (name, g, kind) in cases, (pointname, x) in points(kind; exclude=(:nearminus1,))
        ref = ForwardDiff.hessian(g, x)
        @testset "$name at $pointname" begin
            @test relerr(enzyme_hessian(g, x), ref) < 1e-12
        end
    end

    # The Hessian of `sqrt` needs runtime activity.
    g(x) = sc(sqrt(Q(x)))
    for (pointname, x) in points(:nocut; exclude=(:nearminus1,))
        @testset "sqrt(Q) at $pointname" begin
            H = enzyme_hessian(g, x; runtime_activity=true)
            @test relerr(H, ForwardDiff.hessian(g, x)) < 1e-12
        end
    end
end


@testitem "Enzyme: squad and unflip at batch width 1" tags=[:ad, :enzyme] setup=[ADTestUtils] begin
    # Enzyme's batched reverse mode crashes the compiler on these functions (an upstream
    # bug), so reverse mode is tested only through the gradient of a scalar loss.
    import Enzyme
    using .ADTestUtils
    tin = [0.0, 0.7, 1.5, 2.6]
    x0 = [1.0, 0.1, -0.2, 0.05, 0.9, 0.3, -0.4, 0.1, 0.6, 0.5, -0.5, 0.2, 0.2, 0.7, -0.6, 0.3]
    cases = [
        ("squad", x -> sc(squad([R(x, 1), R(x, 5), R(x, 9), R(x, 13)], tin, x[17])), [x0; 1.1]),
        ("unflip(::Vector{Quaternion})",
         x -> sc(unflip([Q(x, 1), -Q(x, 5), Q(x, 9), -Q(x, 13)])), x0),
        ("unflip(::Vector{Rotor})",
         x -> sc(unflip([R(x, 1), -R(x, 5), R(x, 9), -R(x, 13)])), x0),
    ]
    for (name, g, x) in cases
        ref = bigfd(g, x)
        @testset "$name" begin
            @test relerr(Enzyme.gradient(Enzyme.Reverse, g, x)[1], ref) < 1e-13
            @test relerr(collect(Enzyme.gradient(Enzyme.Forward, g, x)[1]), ref) < 1e-13
        end
    end
end
