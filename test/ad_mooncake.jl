# Tests of Mooncake, which differentiates the package's source natively, and of
# `QuaternionicMooncakeExt`, which imports the ChainRules rules for the two internal
# functions that Mooncake cannot trace: `Quaternionic.csqrt` (complex square roots, used by
# complex quaternions and `Lorentz` rotors) and `Quaternionic.dominant_eigenvector_lapack`
# (used by `from_rotation_matrix` and `align`).  Every gradient is computed in both reverse
# and forward mode, through DifferentiationInterface, and compared with a central finite
# difference evaluated in 256-bit `BigFloat` arithmetic (`bigfd` from `ADTestUtils`).

@testsnippet MooncakeSweep begin
    using Quaternionic
    using Test
    import Mooncake
    import DifferentiationInterface as DI
    using .ADTestUtils: bigfd, relerr

    const MOONCAKE_BACKENDS = (
        "reverse" => DI.AutoMooncake(config=nothing),
        "forward" => DI.AutoMooncakeForward(config=nothing),
    )

    """
        sweep(name, f, pts; tol=1e-10, broken=false)

    Check the gradients of the real function `f` computed by Mooncake in reverse and in
    forward mode at each point `x` of the `name => x` pairs in `pts` against `bigfd`, to the
    relative error `tol`.  Each point and mode is a separate test set, so that an error at
    one does not hide the results at the others.  With `broken=true`, each check is a
    `@test_broken`, which passes when Mooncake either throws or returns a wrong gradient,
    and reports an unexpected pass once the gradient is correct.
    """
    function sweep(name, f, pts; tol=1e-10, broken=false)
        @testset "$name" begin
            for (pname, x) ∈ pts
                ref = bigfd(f, x)
                for (mode, backend) ∈ MOONCAKE_BACKENDS
                    @testset "$name at $pname ($mode)" begin
                        if broken
                            @test_broken relerr(DI.gradient(f, backend, x), ref) < tol
                        else
                            g = DI.gradient(f, backend, x)
                            @test g isa AbstractVector{<:Real}
                            @test relerr(g, ref) < tol
                        end
                    end
                end
            end
        end
    end

    """
        sweep(cases; kwargs...)

    Run `sweep(name, f, pts; kwargs...)` for each `(name, f, pts)` in `cases`.
    """
    function sweep(cases; kwargs...)
        for (name, f, pts) ∈ cases
            sweep(name, f, pts; kwargs...)
        end
    end
end


@testitem "Mooncake: arithmetic and Rotor algebra" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    pts = points()
    sweep([
        ("q*p", x -> sc(Q(x) * P0), pts),
        ("p/q", x -> sc(P0 / Q(x)), pts),
        ("q\\p", x -> sc(Q(x) \ P0), pts),
        ("inv q", x -> sc(inv(Q(x))), pts),
        ("q-t", x -> sc(Q(x) - x[2]) * sc(x[3] - Q(x)), pts),
        ("t*t*t*q*p", x -> sc(x[1] * x[2] * x[3] * Q(x) * P0), pts),
        ("t+q+p+t", x -> sc((x[1] + Q(x) + P0 + x[4]) * P0), pts),
        ("fastmath", x -> sc(@fastmath x[1] * Q(x) + x[2] / Q(x)), pts),
        ("θ+𝐣", x -> sc((x[1] + 𝐣) * Q(x) + 𝐤 * Q(x)), pts),
        ("R*R0", x -> sc(R(x) * R0), pts),
        ("R0/R", x -> sc(R0 / R(x)), pts),
        ("R*R", x -> sc(R(x) * R(x)), pts),
        ("conj R", x -> sc(conj(R(x)) * P0), pts),
        # Off the unit sphere, a `Rotor` stores its components without normalizing them, and
        # the derivative is that of what the source computes from them.
        ("Ru*Ru", x -> sc(Ru(x) * Ru(x)), pts),
        ("P0/Ru", x -> sc(P0 / Ru(x)), pts),
        ("V0/Ru", x -> sc(V0 / Ru(x)), pts),
        ("inv Ru", x -> sc(inv(Ru(x))), pts),
    ])
end


@testitem "Mooncake: exp, log, sqrt, and powers" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    pts, nocut = points(), points(:nocut)
    # `log` and `sqrt` of a `Rotor` are discontinuous at w = 0 off the unit sphere, where
    # finite differences do not apply, so the pure vector is excluded for `Ru`.
    nocutu = points(:nocut; exclude=(:purevector,))
    sweep([
        ("exp q", x -> sc(exp(Q(x))), pts),
        ("exp v", x -> sc(exp(V(x))), pts),
        ("log q", x -> sc(log(Q(x))), nocut),
        ("log R", x -> sc(log(R(x))), nocut),
        ("log Ru", x -> sc(log(Ru(x))), nocutu),
        ("sqrt q", x -> sc(sqrt(Q(x))), nocut),
        ("sqrt R", x -> sc(sqrt(R(x))), nocut),
        ("sqrt v", x -> sc(sqrt(V(x))), points(:nonreal)),
        ("q^n", x -> sc(Q(x)^(-3) + Q(x)^0 + Q(x)^1 + Q(x)^4), pts),
        ("v^2", x -> sc(V(x)^2), pts),
        ("q^0.3", x -> sc(Q(x)^0.3), nocut),
        ("q^t", x -> sc(Q(x)^(x[2] + 0.5)), nocut),
        ("R^0.3", x -> sc(R(x)^0.3), nocut),
        ("Ru^0.3", x -> sc(Ru(x)^0.3), nocutu),
        ("R^-2", x -> sc(R(x)^-2), pts),
    ])
end


@testitem "Mooncake: norms, slerp, distance2, and rotations" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    using LinearAlgebra: normalize
    pts, nonreal = points(), points(:nonreal)
    sweep([
        ("abs", x -> abs(Q(x)), pts),
        ("abs2", x -> abs2(Q(x)), pts),
        ("absvec", x -> absvec(Q(x)), nonreal),
        ("normalize", x -> sc(normalize(Q(x))), pts),
        ("angle R", x -> angle(R(x)), nonreal),
        ("q⋅p", x -> sc(Q(x) ⋅ P0 + Q(x) ⋅ Q(x)), pts),
        ("abs2 R", x -> abs2(R(x)) + abs(Ru(x)), pts),
        ("slerp(R0,R)", x -> sc(slerp(R0, R(x), 0.3)), pts),
        ("slerp(R,R)", x -> sc(slerp(R(x), R(x), 0.3)), pts),
        ("slerp τ", x -> sc(slerp(R0, R(x), x[1])), pts),
        ("distance2", x -> distance2(R(x), R0), pts),
        ("distance2 same", x -> distance2(R(x), R(x) * rotor(1.0, 1e-3, 0, 0)), pts),
        ("distance", x -> distance(R(x), R0), pts),
        ("R(v)", x -> sc(R(x)(V0)), pts),
        ("R(V(x))", x -> sc(R0(V(x))), pts),
    ])
end


@testitem "Mooncake: structure, conversion, and arrays" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    pts = points()
    W = [3, -11, 7, 5] ./ 10
    sweep([
        ("imag and w", x -> imag(Q(x))[1] * Q(x).w, pts),
        ("components", x -> sum(components(R(x) * P0) .* W), pts),
        ("fields", x -> (Q(x) * P0).x + 2 * (V(x) * P0).w, pts),
        ("to_float_array", x -> sum(to_float_array(Q(x) * P0) .* W), pts),
        ("quaternion(vector)", x -> sc(quaternion(x) * P0), pts),
        ("rotor(vector)", x -> sc(rotor(x) * P0), pts),
        ("to_euler_phases", x -> sc(to_euler_phases(R(x))), points(:nonreal)),
        ("bcast .*", x -> sc([Q(x), P0] .* [P0, Q(x)]), pts),
        ("t .* qs", x -> sc(x[1] .* [Q(x), P0]), pts),
        ("abs2.(qs)", x -> sum(abs2.([Q(x), P0])), pts),
        ("map", x -> sc(map(q -> q * P0, [Q(x), R(x)])), pts),
        ("prod", x -> sc(prod([Q(x), P0, Q(x)])), pts),
    ])
end


@testitem "Mooncake: complex components and Lorentz rotors" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    # Most cases here take a complex square root (products of `Lorentz` rotors renormalize,
    # for example), which goes through `Quaternionic.csqrt`.  Without the imported rule,
    # Mooncake refuses to differentiate Base's `sqrt(::Complex)`, because it reinterprets
    # the bits of floating-point numbers.
    pts = points()
    # `C(x)` is a quaternion with complex components built from the four entries of a point.
    C(x) = QC(vcat(x, [x[2] / 3, -x[3] / 2, x[4] / 5 + 0.2, x[1] / 4 - 0.1]))
    B = Boost(0.2, NHAT)
    # At `:minus1`, both `L(x)` and `C(x)` have negative real eigenvalues (as complex 2×2
    # matrices), so they lie on the branch cut of `log`, `sqrt`, and non-integer powers.
    nocut = points(exclude=(:minus1,))
    # `Bv(x)` is a boost given by a velocity, and `phases(x)` perturbs the Euler phases `z`
    # of a generic rotor.
    Bv(x) = Boost(quatvec(x[2] / 3, x[3] / 4, x[4] / 5))
    z = to_euler_phases(rotor(0.6, -0.3, 0.5, 0.4))
    phases(x) = (z[1] * (1 + x[1] / 10), z[2] * complex(1, x[2] / 10), z[3] * complex(1 + x[3] / 10, x[4] / 10))
    # Base's `exp`, `cos`, and `sin` of a complex number take special branches when the real
    # or the imaginary part of the argument is exactly zero, and Mooncake does not
    # differentiate those branches correctly: some of them throw, and others silently give
    # wrong gradients.  The package therefore computes these functions from real functions
    # (`Quaternionic.cexp`, `ccos`, and `csin`).  `exp` of a quaternion with complex
    # components would reach the first branch through its scalar part, and the other two
    # through `absvec`, which is purely imaginary when the vector part has a negative real
    # square, as for the generator of a boost.  Non-integer powers of a pure boost would
    # reach them too: through `exp`, or, for a boost close to the identity raised to a large
    # power, through the closed form that `^` for a `Rotor` uses there.  The cases below
    # that reach these points are regression tests for that change.
    #
    # `C(x)` has a real scalar part and a purely imaginary `absvec` at `:identity` and
    # `:minus1`, and `L(x)` is a pure boost at `:identity`.
    imaginaryabsvec = (:identity, :minus1)
    nopureboost = points(exclude=(:identity, :minus1))
    # `x8` gives a real scalar part and a generic vector part, `x8v` a scalar part with a
    # nonzero imaginary part and a vector part with a negative real square, and `x8g` the
    # generator of a boost; `x8c` is generic.
    x8 = [0.7, 0.3, -0.2, 0.4, 0.0, 0.1, 0.05, -0.2]
    x8v = [0.7, 0.0, 0.0, 0.0, 0.2, 0.3, -0.2, 0.4]
    x8g = vcat(zeros(5), 0.6 .* NHAT)
    x8c = [0.7, 0.3, -0.2, 0.4, 0.13, 0.1, 0.05, -0.2]
    # `VC(x)` is a `QuatVec` with complex components, whose real parts are `x[2:4]` and
    # whose imaginary parts are `x[6:8]`.
    VC(x) = quatvec(complex(x[2], x[6]), complex(x[3], x[7]), complex(x[4], x[8]))
    # `Bp(x)` is a pure boost of rapidity `x[1]`.  For a small power of it, `exp` evaluates
    # `cos` and `sin(a)/a` as series, which are polynomials; larger values of `|s η|` reach
    # the closed forms.
    Bp(x) = Boost(x[1], NHAT)
    sweep([
        ("Lorentz prod", x -> sc(L(x) * B), pts),
        ("Lorentz L*L", x -> sc(L(x) * L(x)), pts),
        ("velocity boost", x -> sc(Bv(x) * R(x)), pts),
        ("complex log", x -> sc(log(L(x))), nocut),
        ("complex sqrt", x -> sc(sqrt(L(x))), nocut),
        ("Lorentz L^s", x -> sc(L(x)^0.3) + sc(L(x)^2.7), nopureboost),
        ("Lorentz L^t", x -> sc(L(x)^(x[2] + 0.5)), nopureboost),
        ("Lorentz small power of a boost", x -> sc(Bp(x)^0.05), ["η=0.6" => [0.6], "η=-1.0" => [-1.0]]),
        ("complex abs", x -> sc(abs(quaternion(L(x)))), pts),
        ("complex inv", x -> sc(inv(L(x))), pts),
        ("conj and ℂconj", x -> sc(conj(L(x))) + sc(Quaternionic.ℂconj(L(x))), pts),
        ("Lorentz action", x -> sum(L(x)([1.0, 0.2, -0.3, 0.5])), pts),
        ("ga_components", x -> sum(ga_components(L(x)) .* (1:8) ./ 8), pts),
        ("RB", x -> sc(Quaternionic.RB(L(x))), pts),
        ("Rv", x -> sc(Quaternionic.Rv(L(x))), pts),
        ("KAN", x -> sc(Quaternionic.KAN(L(x))), pts),
        ("complex q*p", x -> sc(C(x) * P0 + P0 * C(x)), pts),
        ("complex exp", x -> sc(exp(C(x))), points(exclude=imaginaryabsvec)),
        ("complex exp, 8 coordinates", x -> sc(exp(QC(x))), ["x8c" => x8c]),
        ("complex exp of a QuatVec", x -> sc(exp(VC(x))), ["x8c" => x8c]),
        ("complex log q", x -> sc(log(C(x))), nocut),
        ("complex sqrt q", x -> sc(sqrt(C(x))), nocut),
        ("complex abs2", x -> sc(abs2(C(x))), pts),
        ("from_euler_phases", x -> sc(from_euler_phases(phases(x)...)), pts),
    ])
    sweep([
        ("complex exp", x -> sc(exp(C(x))), [p => point(p) for p ∈ imaginaryabsvec]),
        ("complex exp, real scalar part", x -> sc(exp(QC(x))), ["x8" => x8]),
        ("complex exp, imaginary absvec", x -> sc(exp(QC(x))), ["x8v" => x8v]),
        ("exp of a boost generator", x -> sc(exp(QC(x))), ["x8g" => x8g]),
        ("exp of a boost generator as a QuatVec", x -> sc(exp(VC(x))), ["x8g" => x8g]),
        ("Lorentz L^s, pure boost", x -> sc(L(x)^0.3) + sc(L(x)^2.7), [:identity => point(:identity)]),
        ("power of a boost", x -> sc(Bp(x)^x[2]), ["η=1.0, s=0.3" => [1.0, 0.3], "η=-0.8, s=2.7" => [-0.8, 2.7]]),
        ("large power of a small boost", x -> sc(Bp(x)^x[2]), ["η=0.1, s=30" => [0.1, 30.0]]),
        # Near and at the identity, `log` of a `Lorentz` rotor uses the series `logseries`,
        # whose complex square root must go through `csqrt`.  The derivatives are
        # 1.0000009999946 and 1.
        ("log of a Lorentz rotor near the identity",
            x -> (l = log(Lorentz(1.0, x[1] * im, 2e-3, 0)); imag(l[2]) + real(l[3])), ["x=1e-3" => [1e-3]]),
        ("log of a Lorentz rotor at the identity",
            x -> (l = log(Lorentz(1.0, x[1] * im, 0, 0)); imag(l[2]) + real(l[3])), ["x=0" => [0.0]]),
    ])
end


@testitem "Mooncake: from_rotation_matrix and align" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    # The sign of the result of `from_rotation_matrix` and `align` is arbitrary.  Multiplying
    # by the sign that matches the result to `R(x)` gives a smooth loss; the derivative of
    # `sign` is zero.
    signed(Ra, Rb) = sign(components(Ra) ⋅ components(Rb))
    function frm(x)
        Rf = from_rotation_matrix(to_rotation_matrix(R(x)))
        return signed(Rf, R(x)) * sc(Rf)
    end
    a⃗ = [V0, quatvec(0.2, 0.9, -0.4), quatvec(-1.1, 0.1, 0.5)]
    function aligned(x)
        b⃗ = [R(x)(a) for a ∈ a⃗]
        Ra = align(a⃗, b⃗)
        return signed(Ra, R(x)) * sc(Ra)
    end
    # `frm9(y)` differentiates through a matrix that is not exactly orthogonal, built from
    # all 9 entries of `y`.
    M0 = to_rotation_matrix(R0)
    frm9(y) = (Rf = from_rotation_matrix(reshape(y, 3, 3)); signed(Rf, R0) * sc(Rf))
    sweep([
        ("from_rotation_matrix", frm, points()),
        ("align", aligned, points()),
        ("from_rotation_matrix, non-orthogonal", frm9,
         ["perturbed" => vec(collect(M0)) .+ [0.01, -0.02, 0.005, 0.0, 0.015, -0.01, 0.02, 0.0, -0.005]]),
    ]; tol=1e-9)
end


@testitem "Mooncake: Float32" tags=[:ad, :mooncake] setup=[ADTestUtils, MooncakeSweep] begin
    using .ADTestUtils
    # The imported rules are defined for `Float32` as well as `Float64`.  Mooncake's
    # gradients in `Float32` are compared with the `Float64` reference.
    signed(Ra, Rb) = sign(components(Ra) ⋅ components(Rb))
    frm(x) = (Rf = from_rotation_matrix(to_rotation_matrix(R(x))); signed(Rf, R(x)) * sc(Rf))
    # `L` uses the `Float64` direction `NHAT`, which would promote everything to `Float64`.
    T(x) = eltype(x)
    L32(x) = Boost(x[1], T(x).(NHAT)) * Lorentz(rotor(one(T(x)), x[2], x[3], x[4]))
    cases = [
        ("Lorentz L*L", x -> sc(L32(x) * L32(x))),
        ("complex sqrt", x -> sc(sqrt(L32(x)))),
        ("from_rotation_matrix", frm),
    ]
    for (name, f) ∈ cases, (pname, x) ∈ points(:nocut; exclude=(:nearminus1,))
        ref = bigfd(f, x)
        for (mode, backend) ∈ MOONCAKE_BACKENDS
            @testset "$name at $pname ($mode)" begin
                g = DI.gradient(f, backend, Float32.(x))
                @test g isa Vector{Float32}
                @test relerr(g, ref) < 1e-4
            end
        end
    end
end


@testitem "Mooncake: imported rules" tags=[:ad, :mooncake] begin
    using Quaternionic
    using Random: Xoshiro
    using LinearAlgebra: Symmetric, I
    import Mooncake
    using Mooncake.TestUtils: test_rule
    rng = Xoshiro(20261002)
    # `A` is a generic symmetric matrix with a simple largest eigenvalue, as the
    # Bar-Itzhack matrix of a rotation is (its eigenvalues are 3, -1, -1, and -1).
    B = randn(Xoshiro(1), 4, 4)
    A = B + B' + 6I
    for T ∈ (Float32, Float64), mode ∈ (Mooncake.ReverseMode, Mooncake.ForwardMode)
        @testset "csqrt($z) ($T, $mode)" for z ∈ Complex{T}[0.3 + 0.4im, -2 + 1e-3im, 1e-6 - 3im]
            test_rule(rng, Quaternionic.csqrt, z; is_primitive=true, mode=mode)
        end
        @testset "dominant_eigenvector_lapack ($T, $uplo, $mode)" for uplo ∈ ('U', 'L')
            test_rule(
                rng, Quaternionic.dominant_eigenvector_lapack, Matrix{T}(A), uplo;
                is_primitive=true, mode=mode
            )
        end
    end
end


@testitem "Mooncake: test_rule through the source" tags=[:ad, :mooncake] begin
    using Quaternionic
    using Random: Xoshiro
    using LinearAlgebra: normalize
    import Mooncake
    using Mooncake.TestUtils: test_rule
    rng = Xoshiro(20261002)
    # Functions with `Quaternion` and `Rotor` arguments, which Mooncake differentiates by
    # tracing the source.  `QuatVec` arguments cannot be tested this way, because the
    # finite differences of `test_rule` perturb the scalar slot of a `QuatVec`, which every
    # constructor discards.
    q = quaternion(1.2, -0.7, 0.5, 0.3)
    p = quaternion(0.3, 0.5, -0.4, 0.8)
    Ru = Rotor{Float64}(0.5, 0.3, -0.2, 0.4)  # This `Rotor` is off the unit sphere.
    R2 = Rotor{Float64}(normalize([0.2, -0.6, 0.3, 0.7]))
    qc = Quaternion{ComplexF64}(1 + 2im, 0.3, -0.5im, 0.2 + 0.1im)
    cases = [
        ("exp", (exp, q)), ("log", (log, q)), ("sqrt", (sqrt, q)), ("q*p", (*, q, p)),
        ("q/p", (/, q, p)), ("Ru*R", (*, Ru, R2)), ("R^s", (^, R2, 0.3)),
        ("slerp", (slerp, rotor(q), R2, 0.3)), ("distance2", (distance2, rotor(q), R2)),
        ("abs2", (abs2, q)), ("complex product", (*, qc, qc)), ("complex log", (log, qc)),
        ("Lorentz", (t -> Boost(quatvec(t, 0.2, -0.1)) * Lorentz(R2), 0.3)),
        ("from_rotation_matrix", (from_rotation_matrix, to_rotation_matrix(R2))),
    ]
    for (name, args) ∈ cases, mode ∈ (Mooncake.ReverseMode, Mooncake.ForwardMode)
        @testset "$name ($mode)" begin
            test_rule(rng, args...; is_primitive=false, mode=mode)
        end
    end
end


@testitem "Mooncake: friendly tangents" tags=[:ad, :mooncake] setup=[ADTestUtils] begin
    using .ADTestUtils: bigfd, relerr, Q
    using Quaternionic
    import Mooncake
    # With `friendly_tangents=true`, the gradient with respect to a quaternion argument is a
    # `Quaternion` (also for a `Rotor`), or a `QuatVec` for a `QuatVec`.  The hooks exist
    # only in Mooncake 0.5 and later.
    if isdefined(Mooncake, :AsCustomised)
        q = quaternion(1.2, -0.7, 0.5, 0.3)
        Ru = Rotor{Float64}(0.5, 0.3, -0.2, 0.4)
        v = quatvec(0.3, -0.2, 0.5)
        qc = Quaternion{ComplexF64}(1 + 2im, 0.3, -0.5im, 0.2 + 0.1im)
        g(q, r, v) = sum(components(q * r * v * conj(q))) + abs2(q)
        gc(q) = real(sum(components(q * q)))
        garray(qs) = sum(abs2, qs) + sum(components(qs[1] * qs[2]))
        function gradients(f, args...; friendly)
            config = Mooncake.Config(friendly_tangents=friendly)
            cache = Mooncake.prepare_gradient_cache(f, args...; config=config)
            return Mooncake.value_and_gradient!!(cache, f, args...)[2][2:end]
        end
        structural(t::Mooncake.Tangent) = collect(t.fields.components.fields.data)

        ∂q, ∂R, ∂v = gradients(g, q, Ru, v; friendly=true)
        @test ∂q isa QuaternionF64
        @test ∂R isa QuaternionF64
        @test ∂v isa QuatVecF64
        sq, sR, sv = gradients(g, q, Ru, v; friendly=false)
        @test collect(components(∂q)) == structural(sq)
        @test collect(components(∂R)) == structural(sR)
        # The scalar slot of the structural tangent of a `QuatVec` may be nonzero, but it is
        # not a coordinate of the `QuatVec`, so the friendly tangent drops it.
        @test collect(components(∂v))[2:4] == structural(sv)[2:4]
        @test ∂v[1] == 0
        @test relerr(collect(components(∂q)), bigfd(x -> g(Q(x), Ru, v), collect(components(q)))) < 1e-12

        (∂qc,) = gradients(gc, qc; friendly=true)
        @test ∂qc isa Quaternion{ComplexF64}
        @test collect(components(∂qc)) == structural(gradients(gc, qc; friendly=false)[1])

        qs = [q, quaternion(0.1, 0.2, 0.3, 0.4)]
        (∂qs,) = gradients(garray, qs; friendly=true)
        @test ∂qs isa Vector{QuaternionF64}
        @test [collect(components(t)) for t ∈ ∂qs] ==
              [structural(t) for t ∈ gradients(garray, qs; friendly=false)[1]]

        # Integer components have no tangent, and are left to Mooncake's default.
        (∂qi,) = gradients(q -> 2.0 * abs2(q), quaternion(1, 2, 3, 4); friendly=true)
        @test !(∂qi isa AbstractQuaternion)
    else
        @info "Mooncake $(pkgversion(Mooncake)) has no friendly-tangent hooks; skipping"
    end
end
