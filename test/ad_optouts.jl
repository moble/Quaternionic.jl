# Tests that ChainRules' generic `Number` rules do not capture quaternion arguments.
#
# `AbstractQuaternion <: Number`, so the scalar rules of ChainRules, ChainRulesCore, and
# Zygote apply to quaternion arguments unless `QuaternionicChainRulesCoreExt` supplies a
# rule of its own or opts out.  Those generic rules assume that multiplication commutes, so
# a new generic rule in a later release of any of these packages would silently give wrong
# derivatives.  The capture audit below finds every such rule that a call with quaternion
# arguments can reach, and fails for any that is neither displaced nor on an explicit list
# of rules that are correct for quaternions.
#
# An `@opt_out rrule(sig)` also defines `rrule(sig) = nothing`, so an opt-out with exactly
# the signature of a supplied rule would silently delete that rule.  The second item checks
# that every supplied rule still returns a rule, and that every opt-out still opts out.

@testitem "ChainRules: capture audit" tags=[:ad, :chainrules, :slow] begin
    using ChainRulesCore
    import ChainRulesCore: rrule, frule, no_rrule, no_frule, RuleConfig
    import ChainRules, Zygote
    using LinearAlgebra

    const ZC = Zygote.ZygoteRuleConfig{Zygote.Context{false}}
    const QTYPES = (QuaternionF64, RotorF64, QuatVecF64, Quaternion{ComplexF64}, QuatVec{Bool})
    const OTHERTYPES = (Float64, Int, ComplexF64, Bool)

    # Every pattern of 1 to 4 argument types from the lists above with at least one
    # quaternion, grouped by arity
    combos(n) = [c for c ∈ Iterators.product(ntuple(_ -> (QTYPES..., OTHERTYPES...), n)...)
                 if any(T -> T <: AbstractQuaternion, c)]
    const COMBOS = Dict(n => combos(n) for n ∈ 1:4)

    sample(::Type{QuaternionF64}) = quaternion(0.6, -0.31, 0.47, 0.23)
    sample(::Type{RotorF64}) = rotor(0.6, -0.31, 0.47, 0.23)
    sample(::Type{QuatVecF64}) = quatvec(-0.31, 0.47, 0.23)
    sample(::Type{Quaternion{ComplexF64}}) =
        Quaternion{ComplexF64}(0.3 + 0.1im, 0.5 - 0.2im, -0.4 + 0.3im, 0.2 + 0.6im)
    sample(::Type{QuatVec{Bool}}) = 𝐣
    sample(::Type{Float64}) = 0.7
    sample(::Type{Int}) = 2
    sample(::Type{ComplexF64}) = 0.3 + 0.4im
    sample(::Type{Bool}) = true

    # The package (root module) that defines a method
    function rootname(m::Method)
        M = parentmodule(m)
        while parentmodule(M) !== M
            M = parentmodule(M)
        end
        nameof(M)
    end
    isthirdparty(m::Method) = rootname(m) ∈ (:ChainRules, :ChainRulesCore, :Zygote, :ZygoteRules)
    isours(m::Method) = startswith(string(parentmodule(m)), "Quaternionic")
    location(m::Method) = "$(rootname(m)) $(basename(string(m.file))):$(m.line)"

    upper(t) = t isa TypeVar ? upper(t.ub) : t

    # The position of the function argument in the signature of a rule method, which follows
    # a `RuleConfig` for configured rules, and the tangent tuple for `frule`s
    function functionslot(m::Method, kind)
        p = Base.unwrap_unionall(m.sig).parameters
        i = kind === :rrule ? 2 : 3
        if length(p) ≥ 2 && !(p[2] isa Core.TypeofVararg) && upper(p[2]) isa Type && upper(p[2]) <: RuleConfig
            i += 1
        end
        (length(p) < i || p[i] isa Core.TypeofVararg) && return nothing, i
        return upper(p[i]), i
    end

    # The arities (numbers of primal arguments) that a rule method accepts, up to 4
    function arities(m::Method, slot)
        p = Base.unwrap_unionall(m.sig).parameters
        nfixed = length(p) - slot
        if !isempty(p) && p[end] isa Core.TypeofVararg
            return max(nfixed - 1, 1):4
        end
        return 1 ≤ nfixed ≤ 4 ? (nfixed:nfixed) : (1:0)
    end

    # Functions whose primal calls could have side effects are never called: those whose
    # names end in `!`, which mutate their arguments, those of `Base.Threads`, and those
    # named below.  The names are matched whole, because a substring such as `rm` or `eval`
    # also occurs in the names of pure functions such as `norm` and `evalpoly`, whose rules
    # the audit must see.
    const UNSAFE = r"^(show|print|println|printstyled|display|write|read|readline|readlines|readavailable|readuntil|readdir|readlink|rand|randn|randexp|randstring|randperm|randcycle|shuffle|task_local_storage|unsafe_\w+|eval|evalfile|include|include_string|exit|run|rm|cd|mkdir|mkpath|open|close|flush|wait|sleep|lock|unlock|trylock|notify|yield|schedule|ccall)$"
    function unsafe(f)
        f isa Function && parentmodule(f) === Base.Threads && return true
        name = f isa Function ? string(nameof(f)) : string(f)
        return endswith(name, "!") || occursin(UNSAFE, name)
    end

    # Whether the primal call with sample arguments of the types `c` returns.  Calls that
    # have no method, or that inference proves always throw, are not made, because thrown
    # exceptions are slow.
    primalcache = Dict{Any,Bool}()
    function callable(f, c)
        get!(primalcache, (f, c)) do
            unsafe(f) && return false
            hasmethod(f, Tuple{c...}) || return false
            rt = try Core.Compiler.return_type(f, Tuple{c...}) catch; Any end
            rt === Union{} && return false
            try
                f(map(sample, c)...)
                true
            catch
                false
            end
        end
    end

    safewhich(f, T) = try which(f, T) catch; nothing end

    # Whether a method is one of the third-party rules for specific functions (not the
    # fallbacks of ChainRulesCore, whose function argument is `Any`, and which return
    # `nothing`)
    function isgeneric(m::Method, kind)
        isthirdparty(m) || return false
        F, _ = functionslot(m, kind)
        return F !== nothing && F !== Any
    end

    # Sample tangents for the arguments of an `frule`
    sampletangent(::AbstractQuaternion) = quaternion(0.3, -0.2, 0.5, 0.1)
    sampletangent(::Complex) = 0.3 + 0.1im
    sampletangent(::AbstractFloat) = 0.3
    sampletangent(::Integer) = NoTangent()

    # Whether the generic rule that dispatch selects returns no derivative information at
    # all for the sample arguments, as the rules of non-differentiable functions (`typeof`,
    # `isdefined`, predicates, …) do.  Such rules cannot be wrong for quaternions.
    function nondifferentiable(kind, f, c, configured)
        args = map(sample, c)
        try
            if kind === :rrule
                res = configured ? rrule(Zygote.ZygoteRuleConfig(), f, args...) : rrule(f, args...)
                res === nothing && return true
                Ω, pb = res
                return all(∂ -> unthunk(∂) isa AbstractZero, pb(Ω))
            else
                res = frule((NoTangent(), map(sampletangent, args)...), f, args...)
                res === nothing && return true
                return unthunk(res[2]) isa AbstractZero
            end
        catch
            return false
        end
    end

    # The generic rules that are correct for quaternions (section 2 of the specification,
    # "leave alone"): structural functions, rounding, conversions, reductions, selections
    # (`ifelse` and `clamp`), and constructors of arrays, whose rules do not multiply.
    # `ChainRules._unsum` is the internal reshaping step of the pullback of `sum`.  Unary `+`
    # and `*` with one argument are identities.
    const ALLOWED = Any[
        adjoint, identity, float, zero, one, round, floor, ceil, transpose, LinearAlgebra.tr,
        sum, deg2rad, rad2deg, maximum, minimum, findmax, findmin, Base.vect, vcat, hcat,
        hvcat, fill, ifelse, clamp, copy ∘ Base.Broadcast.broadcasted, ChainRules._unsum,
    ]
    function allowed(f, c, configured)
        # Identity comparisons, because `==` of some callable objects is not a `Bool`
        isin(f, fs) = any(g -> g === f, fs)
        isin(f, ALLOWED) && return true
        isin(f, (+, *, Base.FastMath.add_fast, Base.FastMath.mul_fast)) && length(c) == 1 && return true
        # The frule of `dot(x, y)` is `dot(ẋ, y) + dot(x, ẏ)`, which keeps the order of the
        # factors of Base's `dot(x::Number, y::Number) = conj(x) * y`, so it is correct when
        # one argument is not a quaternion.  (The extension supplies the rules for two
        # quaternions.)
        isin(f, (LinearAlgebra.dot,)) && !all(T -> T <: AbstractQuaternion, c) && return true
        # The known gap of section 2: the n-ary sum of ChainRules is reached when the first
        # quaternion is the fourth argument or later (the ninth or later under Zygote, whose
        # extension covers positions 4 through 8).  Only the type of a real argument's
        # cotangent is then wrong.  The Zygote extension does not cover `add_fast`.
        k = findfirst(T -> T <: AbstractQuaternion, c)
        if (!configured && isin(f, (+,))) || isin(f, (Base.FastMath.add_fast,))
            k !== nothing && k ≥ 4 && return true
        end
        # ChainRules' rule for the product of four or more factors multiplies the first three
        # and calls `rrule(*, Ω₃, more...)` on the rest, which reaches the fold rules of the
        # extension when the first quaternion is the fourth factor or later; this is checked
        # against finite differences in `ad_rules.jl`.
        isin(f, (*, Base.FastMath.mul_fast)) && k !== nothing && k ≥ 4 && return true
        return false
    end

    # The functions that have third-party rules, with the arities that those rules accept
    targets = Dict{Tuple{Symbol,Any},Set{Int}}()
    for (kind, rulefn) ∈ ((:rrule, rrule), (:frule, frule)), m ∈ methods(rulefn)
        isthirdparty(m) || continue
        occursin("nondiff", string(m.file)) && continue
        F, slot = functionslot(m, kind)
        (F === nothing || !(F isa DataType) || !isdefined(F, :instance)) && continue
        union!(get!(targets, (kind, F), Set{Int}()), arities(m, slot))
    end

    # For each such function and each pattern of arguments, find the rule that dispatch
    # selects, both for plain ChainRules consumers (`rrule(f, args...)` and
    # `frule(ṫ, f, args...)`) and for Zygote (`rrule(::ZygoteRuleConfig, f, args...)`, which
    # falls back to the plain rule).  A rule of the extension, or one of its opt-outs, is
    # what should be selected.  A generic rule is a problem unless it is allowed, its primal
    # call fails anyway, or it returns no derivative information.
    problems = Set{String}()
    nprobed = Ref(0)
    for ((kind, F), ns) ∈ targets, n ∈ sort!(collect(ns)), c ∈ COMBOS[n]
        f = F.instance
        nprobed[] += 1
        selections = if kind === :rrule
            sel = safewhich(rrule, Tuple{F, c...})
            zsel = safewhich(rrule, Tuple{ZC, F, c...})
            if zsel !== nothing && Base.unwrap_unionall(zsel.sig).parameters[2] === RuleConfig
                zsel = sel
            end
            ((sel, false), (zsel, true))
        else
            ((safewhich(frule, Tuple{Tuple, F, c...}), false),)
        end
        for (s, configured) ∈ selections
            if s === nothing
                # An ambiguity, or no method at all, which is the case only if there is an
                # ambiguity with the fallback
                callable(f, c) || continue
                push!(problems, "$kind $f$(c)$(configured ? " under Zygote" : ""): ambiguous")
            elseif isgeneric(s, kind) && !isours(s) && !allowed(f, c, configured)
                callable(f, c) || continue
                nondifferentiable(kind, f, c, configured) && continue
                push!(problems, "$kind $f$(c)$(configured ? " under Zygote" : ""): generic rule at $(location(s))")
            end
        end
    end
    @test nprobed[] > 10000
    if !isempty(problems)
        println("Generic rules that capture quaternion arguments:")
        foreach(p -> println("  ", p), sort!(collect(problems)))
    end
    @test isempty(problems)
end


@testitem "ChainRules: supplied rules and opt-outs" tags=[:ad, :chainrules] begin
    using ChainRulesCore
    import ChainRulesCore: rrule, frule
    import ChainRules, Zygote
    using LinearAlgebra: LinearAlgebra, norm, normalize, dot
    import Base.FastMath: add_fast, sub_fast, mul_fast, div_fast, inv_fast, exp_fast,
        log_fast, sqrt_fast, pow_fast, angle_fast, sign_fast

    q, p = quaternion(1.2, -0.7, 0.5, 0.3), quaternion(0.3, 0.8, -0.4, 1.1)
    r, r2 = rotor(0.5, 0.3, -0.2, 0.4), rotor(-0.8, 0.6, 0.5, -0.7)
    v = quatvec(0.3, -0.6, 0.2)
    qc = Quaternion{ComplexF64}(0.3 + 0.1im, 0.5 - 0.2im, -0.4 + 0.3im, 0.2 + 0.6im)
    lc = Rotor{ComplexF64}(1.1 + 0.1im, 0.2 - 0.3im, 0.1im, -0.2)
    t, z = 0.7, 0.3 + 0.4im
    A = [3.0 0.1 0.2 0.0; 0.1 -1.0 0.0 0.3; 0.2 0.0 -1.0 0.1; 0.0 0.3 0.1 -1.0]

    # Every signature of the rule table of section 2 (and the rules of section 3) must give a
    # rule, so that no opt-out has overwritten it.
    supplied = Any[
        (+, q, p), (+, q, t), (+, t, q), (-, q, p), (-, q, t), (-, t, q), (*, q, p), (*, q, t),
        (*, t, q), (*, r, r2), (*, qc, z), (*, qc, qc), (/, q, p), (/, q, t), (/, t, q),
        (/, r, r2), (/, q, r), (/, t, r), (\, q, p), (\, q, t), (\, t, q),
        (add_fast, q, p), (sub_fast, q, p), (mul_fast, q, p), (div_fast, q, p),
        (add_fast, q, t), (mul_fast, t, q), (div_fast, q, t),
        (*, q, p, r), (*, t, q, p), (*, q, t, p), (*, q, p, t), (*, t, t, q), (*, t, q, t),
        (*, q, t, t), (*, t, t, q, p), (+, q, p, r), (+, t, q, p), (+, t, t, q, p),
        (mul_fast, q, p, r), (add_fast, q, p, r),
        (-, q), (sub_fast, q), (inv, q), (inv, r), (inv, qc), (inv_fast, q), (conj, q), (conj, qc),
        (exp, q), (exp, qc), (exp, v), (exp, r), (exp, lc), (exp_fast, q),
        (log, q), (log, r), (log_fast, q), (sqrt, q), (sqrt, r), (sqrt, v), (sqrt_fast, q),
        (^, q, 2), (^, q, -2), (^, qc, 2), (^, v, 3), (^, r, 2), (^, q, 0.3), (^, q, 2//3),
        (^, v, 0.3), (^, q, p), (^, q, r), (^, r, p), (^, r, 0.3), (^, r, 2//3), (^, t, q),
        (pow_fast, q, 0.3), (pow_fast, q, 2), (Base.literal_pow, ^, q, Val(2)),
        (Base.literal_pow, ^, r, Val(-1)),
        (muladd, q, p, r), (muladd, q, p, t), (muladd, q, t, p), (muladd, t, q, p),
        (muladd, q, t, t), (muladd, t, q, t), (muladd, t, t, q),
        (real, q), (real, qc), (imag, q), (vec, q), (components, q), (getindex, q, 2),
        (abs, q), (abs, qc), (abs, r), (abs2, q), (abs2, qc), (abs2, r), (abs2, lc),
        (abs2vec, q), (absvec, q), (norm, q), (norm, qc), (norm, r),
        (angle, q), (angle, r), (angle_fast, q), (dot, q, p), (dot, qc, qc),
        (r, v), (distance2, r, r2), (distance, r, r2),
        (normalize, q), (normalize, r), (rotor, q), (rotor, 1.0, 2.0, 3.0, 4.0),
        (rotor, 2.0, 3.0, 4.0), (rotor, 2.0), (rotor, [1.0, 2.0, 3.0, 4.0]), (Rotor, q),
        (quaternion, 1.0, 2.0, 3.0, 4.0), (Quaternion, 1.0), (Quaternion{Float64}, 1.0, 2.0, 3.0, 4.0),
        (Rotor{Float64}, 1.0, 2.0, 3.0, 4.0), (QuatVec{Float64}, 2.0, 3.0, 4.0), (quatvec, v),
        (quaternion, [1.0, 2.0, 3.0, 4.0]), (quatvec, [2.0, 3.0, 4.0]),
        (Quaternionic.csqrt, z), (Quaternionic.dominant_eigenvector_lapack, A, 'U'),
    ]
    # The frules of n-ary `*` and `+` are opted out.
    narywithoutfrule(f, args) = f ∈ (*, +, mul_fast, add_fast) && length(args) ≥ 3
    for (f, args...) ∈ supplied
        @testset "$f$(map(typeof, args))" begin
            @test rrule(f, args...) !== nothing
            ṫ = map(a -> a isa Number ? zero(a) : NoTangent(), args)
            ṫ = map(a -> a isa AbstractQuaternion ? quaternion(zero(a)) : a, ṫ)
            if narywithoutfrule(f, args)
                @test frule((NoTangent(), ṫ...), f, args...) === nothing
            elseif f isa Rotor
                @test frule((quaternion(zero(f)), ṫ...), f, args...) !== nothing
            else
                @test frule((NoTangent(), ṫ...), f, args...) !== nothing
            end
            # Zygote uses the rule (and does not skip it as an opt-out).
            @test first(Zygote.has_chain_rrule(Tuple{Zygote.ZygoteRuleConfig{Zygote.Context{false}}, Core.Typeof(f), map(Core.Typeof, args)...}, Base.get_world_counter()))
        end
    end

    # Every opt-out of section 2 must still opt out, so that the AD system differentiates
    # the source: `rrule` returns `nothing`, and Zygote's `no_rrule` check skips the rule.
    optedout = Any[
        (log, qc), (sqrt, qc), (log_fast, qc), (sqrt_fast, qc), (log, lc), (sqrt, lc),
        (^, qc, 0.3), (^, qc, p), (^, t, qc), (^, lc, 2), (^, lc, 0.3), (pow_fast, qc, 0.3),
        (*, lc, lc), (*, lc, r), (*, r, lc), (/, lc, lc), (/, lc, r), (/, r, lc), (\, lc, r),
        (mul_fast, lc, lc), (div_fast, lc, lc),
        (sign, q), (sign_fast, q), (angle, qc), (angle, lc), (angle_fast, qc), (norm, q, 2), (norm, qc, 1),
        (fma, q, p, r), (fma, q, t, t), (fma, t, t, q), (LinearAlgebra.det, q), (LinearAlgebra.logdet, q),
        (LinearAlgebra.pinv, q),
    ]
    for (f, args...) ∈ optedout
        @testset "opt-out $f$(map(typeof, args))" begin
            @test rrule(f, args...) === nothing
            ṫ = map(a -> a isa AbstractQuaternion ? quaternion(zero(a)) : zero(a), args)
            @test frule((NoTangent(), ṫ...), f, args...) === nothing
            @test !first(Zygote.has_chain_rrule(Tuple{Zygote.ZygoteRuleConfig{Zygote.Context{false}}, Core.Typeof(f), map(Core.Typeof, args)...}, Base.get_world_counter()))
        end
    end
    @test frule((NoTangent(), zero(q)), LinearAlgebra.norm2, q) === nothing
end
