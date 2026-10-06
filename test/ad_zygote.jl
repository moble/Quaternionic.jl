# Tests of Zygote through the package's ChainRules rules and `QuaternionicZygoteExt`.  Every
# gradient is compared with a central finite difference evaluated in 256-bit `BigFloat`
# arithmetic (`bigfd` from `ADTestUtils`), and ForwardDiff, which differentiates the source
# natively, is checked against the same reference.  The cases are those of the sweep that
# validated the AD design, at the point set of `ADTestUtils`.

@testsnippet ZygoteSweep begin
    using Quaternionic
    using Zygote, ForwardDiff
    using Test
    using .ADTestUtils: bigfd, bigfdH, relerr

    """
        zygote_gradient(f, x)

    Return Zygote's gradient of `f` at `x`, with `nothing` (Zygote's zero) replaced by a
    vector of zeros.
    """
    function zygote_gradient(f, x)
        g = Zygote.gradient(f, x)[1]
        return g === nothing ? zero(x) : g
    end

    """
        sweep(name, f, pts; tol=1e-10, forwarddiff=true)

    Check the gradients of the real function `f` computed by Zygote (and by ForwardDiff,
    unless `forwarddiff` is `false`) at each point `x` of the `name => x` pairs in `pts`
    against `bigfd`, to the relative error `tol`.  Each point is a separate test set, so
    that an error at one point does not hide the results at the others.
    """
    function sweep(name, f, pts; tol=1e-10, forwarddiff=true)
        @testset "$name" begin
            for (pname, x) ∈ pts
                @testset "$name at $pname" begin
                    ref = bigfd(f, x)
                    g = zygote_gradient(f, x)
                    @test g isa AbstractVector{<:Real}
                    @test relerr(g, ref) < tol
                    if forwarddiff
                        @test relerr(ForwardDiff.gradient(f, x), ref) < tol
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


@testitem "Zygote: arithmetic" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    sweep([
        ("q*p", x -> sc(Q(x) * P0), pts),
        ("p*q", x -> sc(P0 * Q(x)), pts),
        ("q/p", x -> sc(Q(x) / P0), pts),
        ("p/q", x -> sc(P0 / Q(x)), pts),
        ("q\\p", x -> sc(Q(x) \ P0), pts),
        ("p\\q", x -> sc(P0 \ Q(x)), pts),
        ("inv q", x -> sc(inv(Q(x))), pts),
        ("q+p", x -> sc(Q(x) + P0), pts),
        ("q-p", x -> sc(Q(x) - P0), pts),
        ("-q", x -> sc(-Q(x)), pts),
        ("conj q", x -> sc(conj(Q(x))), pts),
        ("q/t", x -> sc(P0 / x[1] + Q(x) / x[1]), points(exclude=(:purevector,))),
        ("t/q", x -> sc(x[1] / P0 * Q(x)), pts),
        ("t/q 2", x -> sc(x[2] / Q(x)), pts),
        ("real(t*q)", x -> real(x[1] * P0) + x[2]^2, pts),
        ("t*q", x -> sc(x[1] * Q(x) * x[2]), pts),
        ("q-t", x -> sc(Q(x) - x[2]) * sc(x[3] - Q(x)), pts),
        ("muladd", x -> sc(muladd(Q(x), P0, x[2])), pts),
        ("t*t*t*t*q", x -> sc(x[1] * x[2] * x[3] * x[4] * P0), pts),
        ("t*t*t*q*p", x -> sc(x[1] * x[2] * x[3] * Q(x) * P0), pts),
        ("2*3*4*q*p", x -> sc(2.0 * 3.0 * 4.0 * Q(x) * P0), pts),
        ("q*p*q*p", x -> sc(Q(x) * P0 * Q(x) * P0), pts),
        ("t+2t+3t+q", x -> sc(x[1] + 2x[2] + 3x[3] + Q(x)), pts),
        ("t+q+p+t", x -> sc((x[1] + Q(x) + P0 + x[4]) * P0), pts),
        ("fastmath", x -> sc(@fastmath x[1] * Q(x) + x[2]), pts),
        ("fastmath division", x -> sc(@fastmath Q(x) / P0 - x[2] / Q(x)), pts),
        ("exp(θ𝐣)", x -> sc(exp(x[1] * 𝐣 + x[2] * 𝐤)), pts),
        ("θ+𝐣", x -> sc((x[1] + 𝐣) * Q(x)), pts),
        ("θ-2𝐤", x -> sc(x[1] - 2𝐤), pts),
        ("q*𝐢", x -> sc(Q(x) * 𝐢 + 𝐤 * Q(x)), pts),
    ])
end


@testitem "Zygote: Rotor algebra" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    sweep([
        ("R*R0", x -> sc(R(x) * R0), pts),
        ("R0*R", x -> sc(R0 * R(x)), pts),
        ("R0/R", x -> sc(R0 / R(x)), pts),
        ("R/R0", x -> sc(R(x) / R0), pts),
        ("inv R", x -> sc(inv(R(x))), pts),
        ("conj R", x -> sc(conj(R(x))), pts),
        ("R*R", x -> sc(R(x) * R(x)), pts),
        ("R*P0", x -> sc(R(x) * P0), pts),
        ("P0/R", x -> sc(P0 / R(x)), pts),
        # Off the unit sphere, a `Rotor` stores its components without normalizing them, and
        # the derivative is that of what the source computes from them.
        ("Ru*R0", x -> sc(Ru(x) * R0), pts),
        ("R0*Ru", x -> sc(R0 * Ru(x)), pts),
        ("Ru*Ru", x -> sc(Ru(x) * Ru(x)), pts),
        ("P0/Ru", x -> sc(P0 / Ru(x)), pts),
        ("Ru/R0", x -> sc(Ru(x) / R0), pts),
        ("R0/Ru", x -> sc(R0 / Ru(x)), pts),
        ("2/Ru", x -> sc(2.0 / Ru(x)), pts),
        ("V0/Ru", x -> sc(V0 / Ru(x)), pts),
        ("inv Ru", x -> sc(inv(Ru(x))), pts),
    ])
end


@testitem "Zygote: exp, log, sqrt, and powers" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts, nocut = points(), points(:nocut)
    # `log` and `sqrt` of a `Rotor` are discontinuous at w = 0 off the unit sphere, where
    # finite differences do not apply, so the pure vector is excluded for `Ru`.
    nocutu = points(:nocut; exclude=(:purevector,))
    sweep([
        ("exp q", x -> sc(exp(Q(x))), pts),
        ("exp v", x -> sc(exp(V(x))), pts),
        ("exp R", x -> sc(exp(R(x))), pts),
        ("log q", x -> sc(log(Q(x))), nocut),
        ("log R", x -> sc(log(R(x))), nocut),
        ("log Ru", x -> sc(log(Ru(x))), nocutu),
        ("sqrt q", x -> sc(sqrt(Q(x))), nocut),
        ("sqrt R", x -> sc(sqrt(R(x))), nocut),
        ("sqrt Ru", x -> sc(sqrt(Ru(x))), nocutu),
        ("sqrt v", x -> sc(sqrt(V(x))), points(:nonreal)),
        ("q^2", x -> sc(Q(x)^2), pts),
        ("q^3", x -> sc(Q(x)^3), pts),
        ("q^-2", x -> sc(Q(x)^-2), pts),
        ("q^n", x -> sc(Q(x)^(-3) + Q(x)^0 + Q(x)^1 + Q(x)^4), pts),
        ("v^2", x -> sc(V(x)^2), pts),
        ("q^0.3", x -> sc(Q(x)^0.3), nocut),
        ("q^t", x -> sc(Q(x)^(x[2] + 0.5)), nocut),
        ("R^0.3", x -> sc(R(x)^0.3), nocut),
        ("Ru^0.3", x -> sc(Ru(x)^0.3), nocutu),
        ("R^3", x -> sc(R(x)^3), pts),
        ("R^-2", x -> sc(R(x)^-2), pts),
        ("Ru^2", x -> sc(Ru(x)^2), pts),
        ("Ru^-1", x -> sc(Ru(x)^-1), pts),
        ("q^P0", x -> sc(Q(x)^P0), nocut),
        ("2^q", x -> sc(2.0^Q(x)), pts),
        ("R0^q", x -> sc(R0^Q(x)), pts),
    ])
end


@testitem "Zygote: norms, normalization, and angles" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    using LinearAlgebra: normalize
    pts, nonreal = points(), points(:nonreal)
    sweep([
        ("abs", x -> abs(Q(x)), pts),
        ("abs2", x -> abs2(Q(x)), pts),
        ("abs v", x -> abs(V(x)), nonreal),
        ("abs2 v", x -> abs2(V(x)), pts),
        ("absvec", x -> absvec(Q(x)), nonreal),
        ("abs2vec", x -> abs2vec(Q(x)), pts),
        ("normalize", x -> sc(normalize(Q(x))), pts),
        ("rotor", x -> sc(rotor(Q(x))), pts),
        ("Rotor", x -> sc(Rotor(Q(x))), pts),
        ("angle R", x -> angle(R(x)), nonreal),
        ("angle q", x -> angle(Q(x)), nonreal),
        ("q⋅p", x -> Q(x) ⋅ P0 + Q(x) ⋅ Q(x), pts),
        ("abs2 R", x -> abs2(R(x)) + abs(Ru(x)), pts),
    ])
end


@testitem "Zygote: slerp, distance2, and rotations" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    sweep([
        ("slerp(R0,R)", x -> sc(slerp(R0, R(x), 0.3)), pts),
        ("slerp(R,R0)", x -> sc(slerp(R(x), R0, 0.3)), pts),
        ("slerp(R,R)", x -> sc(slerp(R(x), R(x), 0.3)), pts),
        ("slerp τ", x -> sc(slerp(R0, R(x), x[1])), pts),
        ("distance2", x -> distance2(R(x), R0), pts),
        ("distance2 same", x -> distance2(R(x), R(x) * rotor(1.0, 1e-3, 0, 0)), pts),
        ("distance2 q", x -> distance2(Q(x), P0), pts),
        ("distance", x -> distance(R(x), R0), pts),
        ("R(v)", x -> sc(R(x)(V0)), pts),
        ("R(V(x))", x -> sc(R0(V(x))), pts),
        ("RvR̄", x -> sc(R(x) * V0 * conj(R(x))), pts),
    ])
end


@testitem "Zygote: structure and conversion" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    W = [3, -11, 7, 5] ./ 10
    sweep([
        ("imag", x -> sum(imag(Q(x) * P0)), pts),
        ("imag and w", x -> imag(Q(x))[1] * Q(x).w, pts),
        ("vec", x -> sum(vec(Q(x) * P0) .* W[2:4]), pts),
        ("components", x -> sum(Q(x).components .* W), pts),
        ("components()", x -> sum(components(R(x) * P0) .* W), pts),
        ("getindex", x -> (Q(x) * P0)[2] - 3 * (Q(x) * P0)[4], pts),
        ("fields", x -> (Q(x) * P0).x + 2 * (V(x) * P0).w, pts),
        ("to_float_array", x -> sum(to_float_array(Q(x) * P0) .* W), pts),
        ("quaternion(vector)", x -> sc(quaternion(x) * P0), pts),
        ("quatvec(vector)", x -> sc(quatvec(x[2:4]) * P0), pts),
        ("rotor(vector)", x -> sc(rotor(x) * P0), pts),
    ])
end


@testitem "Zygote: broadcasting" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    sweep([
        ("bcast .*", x -> sc([Q(x), P0] .* [P0, Q(x)]), pts),
        ("bcast ./", x -> sc([Q(x), P0] ./ [P0, Q(x)]), pts),
        ("q .* qs", x -> sc(Q(x) .* [P0, Q(x)]), pts),
        ("qs ./ q", x -> sc([P0, Q(x)] ./ Q(x)), pts),
        ("qs .* ts", x -> sc([Q(x), P0] .* [x[2], 2.0]), pts),
        ("ts ./ qs", x -> sc([x[2], 2.0] ./ [P0, Q(x)]), pts),
        ("q .* ts", x -> sc(Q(x) .* [x[2], 2.0]), pts),
        ("t .* qs", x -> sc(x[1] .* [Q(x), P0]), pts),
        ("qs ./ t", x -> sc([Q(x), P0] ./ x[1]), points(exclude=(:purevector,))),
        ("t * qs", x -> sc(x[1] * [Q(x), P0]), pts),
        ("qs / t", x -> sc([Q(x), P0] / x[1]), points(exclude=(:purevector,))),
        ("qs .* true", x -> sc([Q(x), P0] .* true) + sc(false .* [Q(x)]), pts),
        ("qs .^ 2", x -> sc([Q(x), P0] .^ 2), pts),
        ("qs .^ -1", x -> sc([Q(x), P0] .^ -1), pts),
        ("abs2.(qs)", x -> sum(abs2.([Q(x), P0])), pts),
        ("imag.(qs)", x -> sc(imag.([Q(x), P0 * Q(x)])), pts),
        ("qs .+ t", x -> sc([Q(x), P0] .+ x[2]), pts),
        ("conj.(qs)", x -> sc(conj.([Q(x), P0 * Q(x)])), pts),
        ("real.(qs)", x -> sum(real.([Q(x), P0 * Q(x)])), pts),
        ("map", x -> sc(map(q -> q * P0, [Q(x), R(x)])), pts),
        ("sum", x -> sc(sum([Q(x), P0 * Q(x)])), pts),
        ("prod", x -> sc(prod([Q(x), P0, Q(x)])), pts),
        # A broadcast whose operands are all scalars computes the scalar call.
        ("p .* q", x -> sc(P0 .* Q(x)), pts),
        ("q .* p .* q", x -> sc(Q(x) .* P0 .* Q(x)), pts),
        ("p ./ q", x -> sc(P0 ./ Q(x)), pts),
        ("t ./ q", x -> sc(x[2] ./ Q(x)), pts),
        ("q ./ t", x -> sc(Q(x) ./ x[1]), points(exclude=(:purevector,))),
        ("q .^ 2", x -> sc(Q(x) .^ 2), pts),
        ("q .^ 3", x -> sc(Q(x) .^ 3), pts),
        ("q .^ -1", x -> sc(Q(x) .^ -1), pts),
        ("abs2.(Ru)", x -> abs2.(Ru(x)) + abs2.(Q(x) * P0), pts),
        ("imag.(q)", x -> sum(imag.(Q(x) * P0)), pts),
        ("broadcast(*, t, p, q)", x -> sc(broadcast(*, x[1], P0, Q(x))), pts),
        ("broadcast(*, R, R0, Ru)", x -> sc(broadcast(*, R(x), R0, Ru(x))), pts),
        ("L .* B", x -> sc(L(x) .* Boost(0.2, NHAT) ./ L(x)), pts),
        ("L .^ 2", x -> sc(L(x) .^ 2 + L(x) .^ -1), pts),
    ])
end


@testitem "Zygote: complex components and Lorentz rotors" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    pts = points()
    # `C(x)` is a quaternion with complex components built from the four entries of a point.
    C(x) = QC(vcat(x, [x[2] / 3, -x[3] / 2, x[4] / 5 + 0.2, x[1] / 4 - 0.1]))
    B = Boost(0.2, NHAT)
    # At `:minus1`, both `L(x)` and `C(x)` have negative real eigenvalues (as complex 2×2
    # matrices), so they lie on the branch cut of `log` and `sqrt`.
    nocut = points(exclude=(:minus1,))
    sweep([
        ("Lorentz prod", x -> sc(L(x) * B), pts),
        ("Lorentz quot", x -> sc(L(x) / B), pts),
        ("Lorentz L*L", x -> sc(L(x) * L(x)), pts),
        ("Lorentz literal powers", x -> sc(L(x)^2 + L(x)^-1 + L(x)^3), pts),
        ("big boost prod", x -> sc(Boost(3.0 + x[1], NHAT) * Lorentz(rotor(1.0, x[2], x[3], x[4]))), pts),
        ("complex log", x -> sc(log(L(x))), nocut),
        ("complex sqrt", x -> sc(sqrt(L(x))), nocut),
        ("complex abs", x -> sc(abs(quaternion(L(x)))), pts),
        ("complex inv", x -> sc(inv(L(x))), pts),
        ("Lorentz action", x -> sum(L(x)([1.0, 0.2, -0.3, 0.5])), pts),
        ("complex q*p", x -> sc(C(x) * P0 + P0 * C(x)), pts),
        ("complex q/p", x -> sc(C(x) / P0 + P0 / C(x)), pts),
        ("complex inv q", x -> sc(inv(C(x))), pts),
        ("complex exp", x -> sc(exp(C(x))), pts),
        ("complex log q", x -> sc(log(C(x))), nocut),
        ("complex sqrt q", x -> sc(sqrt(C(x))), nocut),
        ("complex q^2", x -> sc(C(x)^2), pts),
        ("complex abs2", x -> sc(abs2(C(x))), pts),
        ("complex conj", x -> sc(conj(C(x)) * C(x)), pts),
        ("complex t*q", x -> sc((x[2] + 0.5im) * C(x)), pts),
        ("complex real and imag", x -> sc(real(C(x))) + sc(sum(imag(C(x)))), pts),
    ])
end


@testitem "Zygote: from_rotation_matrix and align" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    using ChainRulesCore: ignore_derivatives
    # The sign of the result of `from_rotation_matrix` and `align` is arbitrary.  Multiplying
    # by the sign that matches the result to `R(x)` gives a smooth loss.  The sign is
    # treated as a constant.
    signed(Ra, Rb) = ignore_derivatives() do
        c, d = components(Ra), components(Rb)
        sign(Quaternionic.value(c[1] * d[1] + c[2] * d[2] + c[3] * d[3] + c[4] * d[4]))
    end
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
    sweep([
        ("from_rotation_matrix", frm, points()),
        ("align", aligned, points()),
    ]; tol=1e-9)
end


@testitem "Zygote: Hessians" tags=[:ad, :zygote] setup=[ADTestUtils, ZygoteSweep] begin
    using .ADTestUtils
    using LinearAlgebra: normalize
    # `Zygote.hessian` is ForwardDiff over Zygote, so these run every pullback on `Dual`s.
    pts = points(:nocut; exclude=(:nearminus1,))
    nonreal = points(:nonreal; exclude=(:nearminus1,))
    cases = [
        ("q*p", x -> sc(Q(x) * P0 * Q(x)), pts),
        ("q/p", x -> sc(Q(x) / P0 + P0 / Q(x)), pts),
        ("inv q", x -> sc(inv(Q(x))), pts),
        ("t*q", x -> sc(x[1] * Q(x) + x[2] / Q(x)), pts),
        ("exp q", x -> sc(exp(Q(x))), pts),
        ("exp v", x -> sc(exp(V(x))), pts),
        ("log q", x -> sc(log(Q(x))), pts),
        ("sqrt q", x -> sc(sqrt(Q(x))), pts),
        ("q^3", x -> sc(Q(x)^3), pts),
        ("q^0.3", x -> sc(Q(x)^0.3), pts),
        ("R^0.3", x -> sc(R(x)^0.3), pts),
        ("R*R0", x -> sc(R(x) * R0 * R(x)), pts),
        ("slerp(R0,R)", x -> sc(slerp(R0, R(x), 0.3)), pts),
        ("slerp(R,R)", x -> sc(slerp(R(x), R(x), 0.3)), pts),
        ("distance2", x -> distance2(R(x), R0), pts),
        ("abs", x -> abs(Q(x)), pts),
        ("angle R", x -> angle(R(x)), nonreal),
        ("normalize", x -> sc(normalize(Q(x))), pts),
        ("R(v)", x -> sc(R(x)(V0)), pts),
        ("Lorentz prod", x -> sc(L(x) * L(x)), pts),
    ]
    for (name, f, ps) ∈ cases
        @testset "$name" begin
            for (pname, x) ∈ ps
                @testset "$name at $pname" begin
                    ref = bigfdH(f, x)
                    @test relerr(Zygote.hessian(f, x), ref) < 1e-8
                    @test relerr(ForwardDiff.hessian(f, x), ref) < 1e-8
                end
            end
        end
    end
end


@testitem "Zygote: gradient types" tags=[:ad, :zygote] setup=[ADTestUtils] begin
    using .ADTestUtils: sc, bigvjp, relerr, P0, R0, V0
    using Zygote
    # A `Rotor` argument gets a `Quaternion` gradient, never a `Rotor` (which would be
    # renormalized), and it is the ambient gradient of what the source computes.
    f(r) = sc(r * P0)
    g = Zygote.gradient(f, R0)[1]
    @test g isa QuaternionF64
    @test relerr(g, bigvjp(f, R0, 1.0)) < 1e-12
    Ru0 = Rotor{Float64}(0.5, 0.3, -0.2, 0.4)
    g = Zygote.gradient(r -> sc(log(r)), Ru0)[1]
    @test g isa QuaternionF64
    @test relerr(g, bigvjp(r -> sc(log(r)), Ru0, 1.0)) < 1e-12
    # A `QuatVec` argument gets a `QuatVec` gradient.
    h(v) = sc(exp(v)) + sc(v * P0)
    g = Zygote.gradient(h, V0)[1]
    @test g isa QuatVecF64
    @test relerr(g, bigvjp(h, V0, 1.0)) < 1e-12
    # A `Quaternion` argument gets a `Quaternion` gradient, and a `Float32` quaternion gets a
    # `Float32` gradient.
    g = Zygote.gradient(q -> sc(q * P0), P0)[1]
    @test g isa QuaternionF64
    g = Zygote.gradient(q -> sc(q * q), quaternion(1.0f0, 2.0f0, -0.5f0, 0.25f0))[1]
    @test g isa Quaternion{Float32}
    # A real argument gets a real gradient, also when it is added to a quaternion.
    for k ∈ (t -> sc(t * P0), t -> sc(t + P0), t -> sc(P0 - t), t -> sc(t / P0),
                t -> sc(1.0 + 2.0 + t + P0), t -> sc(1.0 + 2.0 + 3.0 + t + P0))
        gk = Zygote.gradient(k, 0.3)[1]
        @test gk isa Float64
        @test gk ≈ bigvjp(k, 0.3, 1.0) rtol=1e-12
    end
    # The constants 𝐢, 𝐣, and 𝐤 have `Bool` components and no gradient, but functions that
    # use them can be differentiated with respect to other arguments.
    @test Zygote.gradient(c -> sc(c * P0), 𝐣)[1] === nothing
    g = Zygote.gradient(θ -> sc(exp(θ * 𝐣)), 0.7)[1]
    @test g ≈ (-3sin(0.7) + 7cos(0.7)) / 10 rtol=1e-14
    # `abs2` and `abs` of a `Rotor` are the constant 1, so their gradients are zero, even
    # off the unit sphere.
    for r ∈ (R0, Ru0)
        @test let gr = Zygote.gradient(abs2, r)[1]
            gr === nothing || iszero(gr)
        end
        @test let gr = Zygote.gradient(abs, r)[1]
            gr === nothing || iszero(gr)
        end
    end
    # The `imag` rule of Zygote's own number rules would give a complex vector here.
    q = quaternion(1.2, -0.7, 0.5, 0.3)
    @test Zygote.gradient(q -> imag(q)[1] * q.w, q)[1] ≈ quaternion(-0.7, 1.2, 0, 0)
end


@testitem "Zygote: constructors and components" tags=[:ad, :zygote] begin
    # These checks are ported from the retired `auto_differentiation.jl`.
    using Zygote
    using StaticArrays: @SVector
    using LinearAlgebra: I
    iszeroish(g) = g === nothing || iszero(g)
    for T ∈ (BigFloat, Float64, Float32)
        w, x, y, z = T(12//10), T(34//10), T(56//10), T(78//10)
        J = Zygote.jacobian((w, x, y, z) -> components(Quaternion(w, x, y, z)), w, x, y, z)
        @test all(J .== eachcol(Matrix{T}(I, 4, 4)))
        J = Zygote.jacobian((w, x, y, z) -> components(QuatVec(w, x, y, z)), w, x, y, z)
        @test iszero(J[1]) && all(J[2:4] .== eachcol(Matrix{T}(I, 4, 4))[2:4])
        n = √(w^2 + x^2 + y^2 + z^2)
        J = Zygote.jacobian((w, x, y, z) -> components(Rotor(w, x, y, z)), w, x, y, z)
        q = [w, x, y, z]
        @test all(J[i] ≈ ([j == i for j ∈ 1:4] .- q[i] .* q ./ n^2) ./ n for i ∈ 1:4)

        # The gradient of `abs2` is checked through each constructor of `Quaternion`.
        for f ∈ (
            (a, b, c, d) -> abs2(Quaternion{T}(@SVector[a, b, c, d])),
            (a, b, c, d) -> abs2(Quaternion{T}([a, b, c, d])),
            (a, b, c, d) -> abs2(Quaternion{T}(a, b, c, d)),
            (a, b, c, d) -> abs2(quaternion(@SVector[a, b, c, d])),
            (a, b, c, d) -> abs2(quaternion([a, b, c, d])),
            (a, b, c, d) -> abs2(quaternion(a, b, c, d)),
            (a, b, c, d) -> abs2(Quaternion([a, b, c, d])),
            (a, b, c, d) -> abs2(Quaternion(a, b, c, d)),
        )
            @test all(Zygote.gradient(f, w, x, y, z) .≈ (2w, 2x, 2y, 2z))
        end
        for f ∈ ((a, b, c, d) -> abs2(quaternion(b, c, d)), (a, b, c, d) -> abs2(Quaternion([b, c, d])))
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test ∇[1] === nothing && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        ∇ = Zygote.gradient((a, b, c, d) -> abs2(quaternion(a)), w, x, y, z)
        @test ∇[1] ≈ 2w && all(isnothing, ∇[2:4])

        # `abs2` of a `Rotor` is the constant 1, whatever its stored components.
        for f ∈ (
            (a, b, c, d) -> abs2(Rotor{T}(a, b, c, d)),
            (a, b, c, d) -> abs2(Rotor{T}([a, b, c, d])),
            (a, b, c, d) -> abs2(rotor(a, b, c, d)),
            (a, b, c, d) -> abs2(rotor([a, b, c, d])),
            (a, b, c, d) -> abs2(Rotor(b, c, d)),
            (a, b, c, d) -> abs2(rotor(a)),
        )
            @test all(iszeroish, Zygote.gradient(f, w, x, y, z))
        end
        # The norm of the normalized components is 1 to within roundoff, so its gradient is
        # zero to within roundoff.
        for f ∈ ((a, b, c, d) -> sum(abs2, components(rotor(a, b, c, d))),
                 (a, b, c, d) -> sum(abs2, components(rotor(b, c, d))))
            @test maximum(g -> g === nothing ? zero(T) : abs(g), Zygote.gradient(f, w, x, y, z)) < 10eps(T)
        end

        # The gradient of `abs2` is checked through each constructor of `QuatVec`, which
        # discards a scalar part.
        for f ∈ (
            (a, b, c, d) -> abs2(QuatVec{T}(@SVector[a, b, c, d])),
            (a, b, c, d) -> abs2(QuatVec{T}([a, b, c, d])),
            (a, b, c, d) -> abs2(QuatVec{T}(a, b, c, d)),
            (a, b, c, d) -> abs2(quatvec([a, b, c, d])),
            (a, b, c, d) -> abs2(quatvec(a, b, c, d)),
            (a, b, c, d) -> abs2(QuatVec(a, b, c, d)),
            (a, b, c, d) -> abs2(quatvec(b, c, d)),
            (a, b, c, d) -> abs2(QuatVec([b, c, d])),
        )
            ∇ = Zygote.gradient(f, w, x, y, z)
            @test iszeroish(∇[1]) && all(∇[2:4] .≈ (2x, 2y, 2z))
        end
        @test all(iszeroish, Zygote.gradient((a, b, c, d) -> abs2(quatvec(a)), w, x, y, z))
    end
end


@testitem "Zygote: exp of a rotation generator" tags=[:ad, :zygote] begin
    # Ported from the retired `auto_differentiation.jl`: rotating 𝐢 by θ about 𝐣 and
    # projecting back onto 𝐢 gives cos θ.
    using Zygote
    for T ∈ (BigFloat, Float64, Float32, Float16)
        f(θ) = exp((θ / 2) * 𝐣)(𝐢) ⋅ 𝐢
        for θ ∈ LinRange(0, 2T(π), 25)
            @test f(θ) ≈ cos(θ) atol=10eps(T)
            @test Zygote.gradient(f, θ)[1] ≈ -sin(θ) atol=10eps(T)
        end
    end
end
