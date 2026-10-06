# Tests of the derivatives of `exp` and `log`, and of powers of a `Rotor`, by every AD
# backend that the package supports, in every mode that the backend has.
#
# Each function is differentiated at a set of points chosen to exercise every branch of the
# source: a generic point, the identity, vector parts of 1e-4, 1e-8, 1e-12, and exactly 0,
# norms of 1e-3 and 1e3, a point near the negative real axis (w < 0 with a vector part of
# 1e-8), and pure vectors.  The full Jacobian of the real coordinates of the output, the
# gradient of the scalar loss `sc`, and (where the backend has second derivatives) the
# Hessian of that loss are compared with central finite differences in 256-bit `BigFloat`
# arithmetic (`bigfd` and `bigfdH` from `ADTestUtils`).  The errors are measured relative to
# the norm of the reference, with no floor, so that the small derivatives at the point of
# norm 1e3 are checked to their full relative accuracy.

@testsnippet ExpLogCases begin
    using Quaternionic
    using Test
    using .ADTestUtils: bigfd, bigfdH, sc, realcoords, relerr, Q, R, V

    # Every point has five entries.  The first four are the components of the quaternion,
    # from which `Q`, `R`, and `V` build their values (`V` ignores the first), and the fifth
    # is the exponent of the power of a rotor, which the other functions ignore.
    const n̂ = [0.6, -0.48, 0.64]
    const EXPLOG_POINTS = [
        "generic" => [1.2, -0.7, 0.5, 0.3, 0.37],
        "identity" => [1.0, 0.0, 0.0, 0.0, 0.37],
        "vector part 1e-4" => [0.9; 1e-4 * n̂; 0.37],
        "vector part 1e-8" => [0.9; 1e-8 * n̂; 0.37],
        "vector part 1e-12" => [0.9; 1e-12 * n̂; 0.37],
        "vector part 0" => [0.9, 0.0, 0.0, 0.0, 0.37],
        "norm 1e-3" => [6e-4, 5e-4, -4e-4, 4.8e-4, 0.37],
        "norm 1e3" => [60.0, 600.0, -400.0, 690.0, 0.37],
        "near the negative real axis" => [-0.8; 1e-8 * n̂; 0.37],
        "pure unit vector" => [0.0; n̂; 0.37],
        "pure vector of norm 2" => [0.0; 2n̂; -1.3],
    ]

    # The rotor functions take the normalized `R(x)`, because `log` and powers of a `Rotor`
    # assume unit norm.
    const EXPLOG_CASES = [
        "exp(::Quaternion)" => x -> exp(Q(x)),
        "exp(::QuatVec)" => x -> exp(V(x)),
        "log(::Quaternion)" => x -> log(Q(x)),
        "log(::Rotor)" => x -> log(R(x)),
        "Rotor^s" => x -> R(x)^x[5],
    ]

    """
        checkexplog(modes; T=Float64, rtol, cases=EXPLOG_CASES, points=EXPLOG_POINTS, broken)

    Check each mode of a backend for each case at each point, with components of type `T`.
    Each mode is `label => (kind, d)`, where `kind` is `:jacobian`, `:gradient`, or
    `:hessian`, and `d(F, x)` returns the Jacobian of the vector function `F`, the gradient
    of the real function `F`, or its Hessian, at `x`.  The Jacobian is that of the real
    coordinates of the output, and the gradient and Hessian are those of the loss `sc`.  The
    tolerance `rtol` is a function of `kind` and `T`.  The references are computed once for
    each case and point, and only for the kinds that are checked.  The check of a mode is a
    `@test_broken` where `broken(case, point, mode)` is `true`, for the names of the case,
    the point, and the mode.
    """
    function checkexplog(modes; T=Float64, rtol=defaulttol, cases=EXPLOG_CASES, points=EXPLOG_POINTS,
                         broken=(c, p, m) -> false)
        kinds = unique(first(last(m)) for m in modes)
        for (cname, f) in cases
            F = x -> realcoords(f(x))
            g = x -> sc(f(x))
            @testset "$cname at $pname" for (pname, x64) in points
                x = T.(x64)
                refs = Dict(
                    k => (k === :jacobian ? bigfd(f, x) : k === :gradient ? bigfd(g, x) : bigfdH(g, x))
                    for k in kinds
                )
                for (label, (kind, d)) in modes
                    @testset "$label" begin
                        if broken(cname, pname, label)
                            @test_broken relerr(d(kind === :jacobian ? F : g, x), refs[kind], 0) < rtol(kind, T)
                        else
                            @test relerr(d(kind === :jacobian ? F : g, x), refs[kind], 0) < rtol(kind, T)
                        end
                    end
                end
            end
        end
    end

    # The tolerances relative to the norm of the reference.  Near the negative real axis,
    # the derivatives of `log` and the power are of order 1e8 (and the second derivatives
    # of order 1e16), and the roundoff in the vector part grows with them.
    defaulttol(kind, ::Type{Float64}) = kind === :hessian ? 1e-10 : 1e-12
    defaulttol(kind, ::Type{Float32}) = kind === :hessian ? 1e-3 : 1e-4
end


@testitem "exp and log: ForwardDiff" tags=[:slow, :ad, :forwarddiff] setup=[ADTestUtils, ExpLogCases] begin
    import ForwardDiff
    modes = [
        "jacobian" => (:jacobian, ForwardDiff.jacobian),
        "gradient" => (:gradient, ForwardDiff.gradient),
        "hessian" => (:hessian, ForwardDiff.hessian),
        "jacobian of gradient" => (:hessian, (g, x) -> ForwardDiff.jacobian(y -> ForwardDiff.gradient(g, y), x)),
    ]
    checkexplog(modes)
    checkexplog(modes[1:2]; T=Float32)
end


@testitem "exp and log: ReverseDiff" tags=[:slow, :ad, :reversediff] setup=[ADTestUtils, ExpLogCases] begin
    import ReverseDiff
    # A tape records the branches taken at the point where it is recorded, so each tape is
    # recorded at the point where it is evaluated.
    compiledgradient(g, x) = ReverseDiff.gradient!(similar(x), ReverseDiff.compile(ReverseDiff.GradientTape(g, x)), x)
    compiledjacobian(F, x) = ReverseDiff.jacobian!(ReverseDiff.compile(ReverseDiff.JacobianTape(F, x)), x)
    compiledhessian(g, x) = ReverseDiff.hessian!(ReverseDiff.compile(ReverseDiff.HessianTape(g, x)), x)
    modes = [
        "jacobian" => (:jacobian, ReverseDiff.jacobian),
        "gradient" => (:gradient, ReverseDiff.gradient),
        "hessian" => (:hessian, ReverseDiff.hessian),
        "compiled jacobian tape" => (:jacobian, compiledjacobian),
        "compiled gradient tape" => (:gradient, compiledgradient),
        "compiled hessian tape" => (:hessian, compiledhessian),
    ]
    checkexplog(modes)
    checkexplog(modes[1:2]; T=Float32)
end


@testitem "exp and log: Zygote" tags=[:slow, :ad, :zygote] setup=[ADTestUtils, ExpLogCases] begin
    import Zygote
    modes = [
        "jacobian" => (:jacobian, (F, x) -> Zygote.jacobian(F, x)[1]),
        "gradient" => (:gradient, (g, x) -> Zygote.gradient(g, x)[1]),
        "hessian (forward over reverse)" => (:hessian, Zygote.hessian),
    ]
    checkexplog(modes)
    checkexplog(modes[1:2]; T=Float32)
end


@testitem "exp and log: ChainRules rules" tags=[:slow, :ad, :chainrules] setup=[ADTestUtils, ExpLogCases] begin
    using ChainRulesCore: frule, rrule, NoTangent, unthunk
    using ChainRulesTestUtils: test_frule, test_rrule
    using LinearAlgebra: I
    import Random

    # The functions with rules, applied to their primal arguments, which are built from a
    # point by `args`.  A `Rotor` argument is normalized, and its tangents are `Quaternion`s.
    rules = [
        ("exp(::Quaternion)", exp, x -> (Q(x),)),
        ("exp(::QuatVec)", exp, x -> (V(x),)),
        ("log(::Quaternion)", log, x -> (Q(x),)),
        ("log(::Rotor)", log, x -> (R(x),)),
        ("Rotor^s", ^, x -> (R(x), x[5])),
    ]

    # The coordinates of an argument, from which the references rebuild it, and the basis
    # of its tangents
    coordsof(a::QuatVec) = collect(vec(a))
    coordsof(a::AbstractQuaternion) = collect(components(a))
    coordsof(a::Real) = [a]
    basis(a::QuatVec) = [quatvec((1:3 .== k)...) for k in 1:3]
    basis(a::AbstractQuaternion) = [quaternion((1:4 .== k)...) for k in 1:4]
    basis(a::Real) = [one(a)]
    zerotangent(a::QuatVec) = zero(a)
    zerotangent(a::AbstractQuaternion) = zero(quaternion(a))
    zerotangent(a::Real) = zero(a)
    # Rebuild an argument of the type of `a` from coordinates.  A `Rotor` is rebuilt with
    # `rotor`, which normalizes, so that the reference is the derivative along the unit
    # sphere.  The pushforward of a rotor tangent is compared after projection onto the
    # sphere, as `projector` below gives it.
    rebuild(a::QuatVec, c) = quatvec(c...)
    rebuild(a::Rotor, c) = rotor(c...)
    rebuild(a::Quaternion, c) = quaternion(c...)
    rebuild(a::Real, c) = c[1]
    projector(a::Rotor) = I - collect(components(a)) * collect(components(a))'
    projector(a) = I

    @testset "$name at $pname" for (name, f, args) in rules, (pname, x) in EXPLOG_POINTS
        a = args(x)
        y = f(a...)
        # The reference Jacobian with respect to each argument, and its projection
        Js = map(eachindex(a)) do i
            c = coordsof(a[i])
            J = bigfd(cᵢ -> f(Base.setindex(a, rebuild(a[i], cᵢ), i)...), c)
            J = J isa AbstractVector ? reshape(J, :, 1) : J
            J * projector(a[i])
        end
        # The pushforwards of the basis tangents of each argument
        for i in eachindex(a)
            columns = map(basis(a[i])) do ȧ
                ȧs = ntuple(j -> j == i ? ȧ : zerotangent(a[j]), length(a))
                Ω, Ω̇ = frule((NoTangent(), ȧs...), f, a...)
                @test Ω == y
                Ω̇ = unthunk(Ω̇)
                realcoords(y isa QuatVec ? quatvec(Ω̇) : quaternion(Ω̇))
            end
            @test relerr(reduce(hcat, columns) * projector(a[i]), Js[i], 0) < 1e-12
        end
        # The pullbacks of the basis cotangents of the output
        Ω, pullback = rrule(f, a...)
        @test Ω == y
        ȳs = y isa QuatVec ? basis(y) : basis(quaternion(y))
        for i in eachindex(a)
            rows = map(ȳs) do ȳ
                ∂ = unthunk(pullback(ȳ)[i+1])
                c = a[i] isa QuatVec ? collect(vec(∂)) : a[i] isa Real ? [∂] : collect(components(∂))
                permutedims(c)
            end
            @test relerr(reduce(vcat, rows) * projector(a[i]), Js[i], 0) < 1e-12
        end
    end

    # ChainRulesTestUtils, except where its Float64 finite differences are not accurate
    # enough: near the negative real axis, where they straddle the branch cut; at the norm
    # 1e-3, where the derivatives of `log` are of order 1e4 and the errors of the finite
    # differences exceed the tolerance; and at w = 0 for a `Rotor`, where they straddle the
    # branches of `log` and powers of an unnormalized `Rotor`.
    @testset "ChainRulesTestUtils for $name at $pname" for (name, f, args) in rules, (pname, x) in EXPLOG_POINTS
        pname ∈ ("near the negative real axis", "norm 1e-3") && continue
        a = args(x)
        a[1] isa Rotor && iszero(a[1][1]) && continue
        # ChainRulesTestUtils draws its random tangents from the global generator, which is
        # reseeded so that each check is reproducible on its own.
        Random.seed!(20261002)
        test_frule(f, a...; rtol=1e-8, atol=1e-8)
        test_rrule(f, a...; rtol=1e-8, atol=1e-8)
    end
end


@testitem "exp and log: Enzyme" tags=[:slow, :ad, :enzyme] setup=[ADTestUtils, ExpLogCases, EnzymeHelpers] begin
    import Enzyme
    modes = [
        "forward jacobian" => (:jacobian, (F, x) -> Enzyme.jacobian(Enzyme.Forward, F, x)[1]),
        "reverse jacobian" => (:jacobian, (F, x) -> Enzyme.jacobian(Enzyme.Reverse, F, x)[1]),
        "forward gradient" => (:gradient, (g, x) -> collect(Enzyme.gradient(Enzyme.Forward, g, x)[1])),
        "reverse gradient" => (:gradient, (g, x) -> Enzyme.gradient(Enzyme.Reverse, g, x)[1]),
        "forward over reverse hessian" => (:hessian, (g, x) -> enzyme_hessian(g, x)),
    ]
    checkexplog(modes)
    checkexplog(modes[3:4]; T=Float32)
end


@testitem "exp and log: Mooncake" tags=[:slow, :ad, :mooncake] setup=[ADTestUtils, ExpLogCases] begin
    import Mooncake
    import DifferentiationInterface as DI
    reverse, forward = DI.AutoMooncake(config=nothing), DI.AutoMooncakeForward(config=nothing)
    fwdoverrev = DI.SecondOrder(forward, reverse)
    # Mooncake builds its rules for the types of the arguments, not their values, so one
    # preparation for each function, operator, backend, and element type serves every point.
    preps = IdDict()
    function prepared(op, prepare, backend)
        (f, x) -> op(f, get!(() -> prepare(f, backend, x), preps, (f, op, backend, eltype(x))), backend, x)
    end
    modes = [
        "reverse jacobian" => (:jacobian, prepared(DI.jacobian, DI.prepare_jacobian, reverse)),
        "forward jacobian" => (:jacobian, prepared(DI.jacobian, DI.prepare_jacobian, forward)),
        "reverse gradient" => (:gradient, prepared(DI.gradient, DI.prepare_gradient, reverse)),
        "forward gradient" => (:gradient, prepared(DI.gradient, DI.prepare_gradient, forward)),
        "forward over reverse hessian" => (:hessian, prepared(DI.hessian, DI.prepare_hessian, fwdoverrev)),
    ]
    # Mooncake 0.5.60 has two bugs in forward mode over reverse mode, both reproduced in
    # `mooncake_fwd_over_rev_zero.jl` in the upstream reports.  First, the Hessian of a
    # product with a factor that is exactly zero, such as x[2] * exp(x[1]) at x[2] = 0, is
    # wrong, which affects `log` at a zero vector part.  Second, the Hessian of
    # `atan(y, x)` at x = 0 throws an `ArgumentError` about `bitcast`, which affects `log` and
    # powers of pure vectors.
    zerovector = ("identity", "vector part 0")
    purevector = ("pure unit vector", "pure vector of norm 2")
    function broken(case, point, mode)
        mode == "forward over reverse hessian" || return false
        startswith(case, "log") && point ∈ (zerovector..., purevector...) ||
            case == "Rotor^s" && point ∈ purevector
    end
    # Every Mooncake rule is compiled for each function and mode, which dominates the time
    # of this item, so the gradients (which the Jacobians already check) are computed only
    # with `Float32` components.
    checkexplog(modes[[1, 2, 5]]; broken=broken)
    checkexplog(modes[3:4]; T=Float32)
end


@testitem "exp and log: Symbolics and FastDifferentiation" tags=[:slow, :ad] setup=[ADTestUtils, ExpLogCases] begin
    import Symbolics
    import FastDifferentiation

    # Symbolic components cannot be compared with thresholds, so `exp`, `log`, and powers of
    # a `Rotor` return their general formulas, which have removable singularities at a zero
    # vector part.  Their derivatives are evaluated at points whose vector parts are not
    # small.  The symbolic normalization of `R(x)` makes the expressions for the power so
    # large that their Jacobian takes a minute to compute, so the power is taken of the
    # unnormalized `Ru(x)`, at points normalized beforehand.  (The general formula for the
    # power, like the source away from the identity, does not depend on the norm.)
    using .ADTestUtils: Ru
    using LinearAlgebra: normalize
    Symbolics.@variables w x y z s
    vars = [w, x, y, z, s]
    numeric = ["generic", "norm 1e3", "pure unit vector", "pure vector of norm 2"]
    cases = [EXPLOG_CASES[1:4]; "Rotor^s" => (x -> Ru(x)^x[5])]
    for (cname, f) in cases
        J = Symbolics.jacobian(realcoords(f(vars)), vars)
        Jf = Symbolics.build_function(J, vars; expression=Val(false))[1]
        @testset "Symbolics: $cname at $pname" for (pname, p) in EXPLOG_POINTS
            pname ∈ numeric || continue
            p = cname == "Rotor^s" ? [normalize(p[1:4]); p[5]] : p
            @test relerr(Jf(p), bigfd(f, p), 0) < 1e-12
        end
    end

    # FastDifferentiation cannot follow the branches of `exp`, `log`, and powers, so they
    # raise explanatory errors.
    u = FastDifferentiation.make_variables(:u, 5)
    for (cname, f) in EXPLOG_CASES
        @test_throws ErrorException f(u)
    end
end
