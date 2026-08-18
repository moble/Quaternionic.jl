# Property-based tests using Supposition.jl.
#
# These complement — they do not replace — the hand-built input grids in
# `algebra.jl`, `conversion.jl`, etc.  Those grids are exhaustive over a small
# set of exactly representable components, which is a genuinely different
# guarantee.  What Supposition adds is:
#
#   * an adaptive search over a much larger input space,
#   * automatic shrinking of any counterexample to a minimal reproducer,
#   * `target!`, which hill-climbs toward the *worst* case of a stated error
#     measure rather than hoping a random draw lands near a branch cut, and
#   * a failure database, so a counterexample found once is replayed first on
#     every later run.
#
# Conventions used here:
#
#   * Every tolerance is `eps(T)`-relative; nothing is hard-coded to Float64.
#   * Generators are built from `Data.Integers`, not `Data.Floats`, because the
#     latter is restricted to `Base.IEEEFloat` and so cannot reach `Double64`
#     or `BigFloat`.  See `PropertyGens.floats`.
#   * Properties are written as *named predicate functions* in the
#     `PropertyGens` test module, and invoked with `@check`'s call syntax.
#     This is not stylistic: `@check function foo(...)` defines a global
#     method, so it cannot appear inside a `for T ∈ FloatTypes` loop.  The
#     call syntax can, which is what lets one property sweep every precision.

@testmodule PropertyGens begin
    using Quaternionic
    using Supposition, Supposition.Data
    using DoubleFloats: Double64
    using Random: Xoshiro, RandomDevice

    # ── Configuration ─────────────────────────────────────────────────────────

    "Float types swept by the property tests."
    const FloatTypes = (Float32, Float64, Double64)

    """
    Float types for properties that route through `LinearAlgebra.eigen`
    (`from_rotation_matrix`, `align`).

    `Double64` is excluded because neither `GenericLinearAlgebra` nor
    `GenericSchur` provides `eigen` for `Symmetric{Double64,<:SMatrix{4,4}}`, nor
    the `eigen!(::Symmetric{Double64,<:Matrix}, ::UnitRange)` that `align` needs;
    both raise `MethodError`.  `BigFloat` is excluded because it silently returns
    the *wrong* rotor — see the `:broken` test item at the bottom of this file.
    """
    const EigenFloatTypes = (Float32, Float64)

    """
        CFG

    Shared `@check` configuration, passed as `config=CFG`.  Individual `@check`
    invocations may still override any single option (`broken=true`, a larger
    `max_examples`, …).

    * **RNG.** Fixed by default, so CI is reproducible and a reported
      counterexample can be reproduced locally.  Set
      `QUATERNIONIC_PROPERTY_SEED=random` in the environment for an exploratory
      run that searches somewhere new, or to any integer to replay a seed.
    * **Example count.** Modest, so the whole file stays in the `:fast` tier
      (~10 s for ~90 properties).
    * **Failure database.** Supposition records each counterexample it finds and
      replays it first on the next run.  Its default directory is the one holding
      the active `Project.toml`; in this workspace layout that is `test/` itself,
      so the default scatters one record file per property directly into the
      source tree.  Naming a subdirectory explicitly keeps them together.
    """
    const CFG = Supposition.CheckConfig(;
        rng = let s = get(ENV, "QUATERNIONIC_PROPERTY_SEED", "42")
            s == "random" ? Xoshiro(rand(RandomDevice(), UInt)) : Xoshiro(parse(UInt, s))
        end,
        max_examples = 250,
        db = Supposition.DirectoryDB(mkpath(joinpath(@__DIR__, "SuppositionDB"))),
    )

    # ── Generators ────────────────────────────────────────────────────────────

    """
        floats(T, N=1024)

    Dyadic rationals `k/N` for `|k| ≤ N`, shrinking toward zero.  `N` is a power
    of two, so every value is exactly representable in any binary float type and
    the generated set is genuinely the same across precisions.

    Deliberately *not* `Data.Floats{T}`: that samples uniformly over bit
    patterns (so almost every draw is around `1e-300` — useless for geometry)
    and only supports `Float16`/`Float32`/`Float64`.  Use `wildfloats` where
    bit-pattern coverage is what you actually want.
    """
    floats(::Type{T}, N=1024) where {T} = map(n -> T(n) / N, Data.Integers(-N, N))

    """
        wildfloats(T)

    Uniform over *bit patterns* with `|x| ≤ 1`: denormals, `±0`, `1e-300`.  Only
    available for `Base.IEEEFloat`.  Appropriate for scale-free algebraic
    identities, and a poor choice for anything involving angles.
    """
    wildfloats(::Type{T}) where {T<:Base.IEEEFloat} =
        filter(x -> abs(x) ≤ 1, Data.Floats{T}(; nans=false, infs=false))

    quats(g) = map(v -> quaternion(v...), Data.Vectors(g; min_size=4, max_size=4))
    quats(::Type{T}) where {T} = quats(floats(T))

    "Rotors, rejecting inputs too close to zero to normalize meaningfully."
    function rotors(::Type{T}) where {T}
        map(rotor, filter(q -> abs2(q) > sqrt(eps(T)), quats(T)))
    end

    quatvecs(::Type{T}) where {T} =
        map(v -> quatvec(v[1], v[2], v[3]),
            Data.Vectors(floats(T); min_size=3, max_size=3))

    "Interpolation parameters τ ∈ [0, 1]."
    unitinterval(::Type{T}, N=64) where {T} = map(n -> T(n) / N, Data.Integers(0, N))

    # ── Tolerances ────────────────────────────────────────────────────────────

    "Absolute tolerance of `n` ulps at the precision of `q`."
    tol(q, n=10) = n * eps(basetype(q))

    # ── Predicates: normed division algebra ───────────────────────────────────

    associativity(p, q, r) = isapprox((p * q) * r, p * (q * r); atol=tol(p, 100))
    left_distributivity(p, q, r) = isapprox(r * (p + q), r * p + r * q; atol=tol(p, 100))
    right_distributivity(p, q, r) = isapprox((p + q) * r, p * r + q * r; atol=tol(p, 100))
    conj_involution(q) = conj(conj(q)) == q
    conj_antiautomorphism(p, q) = isapprox(conj(p * q), conj(q) * conj(p); atol=tol(p, 100))

    "‖pq‖² = ‖p‖²‖q‖² — the composition property that makes ℍ a normed algebra."
    norm_composition(p, q) =
        isapprox(abs2(p * q), abs2(p) * abs2(q); rtol=tol(p, 100), atol=tol(p, 100))

    function inverse_law(q)
        abs2(q) > sqrt(eps(basetype(q))) || return true   # not invertible in practice
        isapprox(q * inv(q), one(q); atol=tol(q, 100)) &&
            isapprox(inv(q) * q, one(q); atol=tol(q, 100))
    end

    # ── Predicates: round-trips ───────────────────────────────────────────────

    euler_roundtrip(R) =
        distance(rotor(from_euler_angles(to_euler_angles(R)...)), R) ≤ tol(R, 100)

    euler_phase_roundtrip(R) =
        distance(rotor(from_euler_phases(to_euler_phases(R))), R) ≤ tol(R, 100)

    matrix_roundtrip(R) =
        distance(from_rotation_matrix(to_rotation_matrix(R)), R) ≤ tol(R, 100)

    """
    `(θ, ϕ)` fixes only the image of ẑ, not the full rotation, so the round-trip
    is stated on that image rather than on the rotor.
    """
    function spherical_roundtrip(R)
        θ, ϕ = to_spherical_coordinates(R)
        S = from_spherical_coordinates(θ, ϕ)
        isapprox(to_rotation_matrix(S)[:, 3], to_rotation_matrix(R)[:, 3]; atol=tol(R, 100))
    end

    """
    `to_float_array` of a scalar quaternion returns a `Vector`, and
    `from_float_array` of a length-4 `Vector` returns a *zero-dimensional* array
    of quaternions — so the round-trip needs a `[]` to get back to a scalar.
    """
    floatarray_roundtrip(q) = from_float_array(to_float_array(q))[] === q

    # ── Predicates: exp, log, sqrt, powers ────────────────────────────────────

    function explog(q)
        abs2(q) > sqrt(eps(basetype(q))) || return true
        err = distance(exp(log(q)), q) / abs(q)
        target!(Float64(err))          # steer the search toward the worst case
        err ≤ tol(q, 100)
    end

    logexp_rotor(R) = distance(exp(log(R)), R) ≤ tol(R, 100)

    function sqrt_squares(q)
        abs2(q) > sqrt(eps(basetype(q))) || return true
        err = distance(sqrt(q)^2, q) / abs(q)
        target!(Float64(err))
        err ≤ tol(q, 100)
    end

    "Rⁿ Rᵐ = R⁽ⁿ⁺ᵐ⁾ for integer powers."
    pow_addition(R, m, n) = distance(R^m * R^n, R^(m + n)) ≤ tol(R, 100)

    # ── Predicates: rotation semantics ────────────────────────────────────────

    "The matrix and the sandwich product must implement the *same* rotation."
    sandwich_matches_matrix(R, v) =
        isapprox(to_rotation_matrix(R) * vec(v), vec(R * v * conj(R)); atol=tol(R, 100))

    rotation_is_isometry(R, v) =
        isapprox(absvec(R * v * conj(R)), absvec(v); atol=tol(R, 100))

    "R ↦ ℛ(R) is a homomorphism onto SO(3)."
    matrix_homomorphism(R, S) =
        isapprox(to_rotation_matrix(R * S), to_rotation_matrix(R) * to_rotation_matrix(S);
                 atol=tol(R, 100))

    "±R are the two preimages of one rotation."
    double_cover(R, v) = isapprox(R * v * conj(R), (-R) * v * conj(-R); atol=tol(R, 100))

    # ── Predicates: interpolation ─────────────────────────────────────────────

    slerp_endpoints(R, S) =
        distance(slerp(R, S, 0), R) ≤ tol(R, 100) && distance(slerp(R, S, 1), S) ≤ tol(R, 100)

    slerp_stays_unit(R, S, τ) = abs(abs(slerp(R, S, τ)) - 1) ≤ tol(R, 100)

    """
    slerp traverses the geodesic at constant angular speed:
    `distance(R, slerp(R, S, τ)) == τ · distance(R, S)`.

    This needs `unflip=true`.  Without it, slerp takes the *long* way around
    whenever `R⋅S < 0`, while `distance` always reports the short way — so the
    identity fails for half of all rotor pairs.  Supposition walks straight to
    that boundary: stated without `unflip`, it shrinks to a pair with `R⋅S`
    just barely negative (≈ -2e-4 when this was first written).
    """
    function slerp_constant_speed(R, S, τ)
        d = distance(R, S)
        d > sqrt(eps(basetype(R))) || return true
        err = abs(distance(R, slerp(R, S, τ; unflip=true)) - τ * d)
        target!(Float64(err))
        err ≤ tol(R, 100) * max(d, 1)
    end

    """
    The `unflip=false` default really does take the long way when R⋅S < 0.

    Written as `a - b ≤ ε` rather than the more natural `a ≤ b + ε`, because the
    latter is not reliable at `Double64`: adding a tolerance of order `100eps` to
    a value of order `0.05` renormalizes the limbs, and `DoubleFloats`' `≤` then
    compares the result as *smaller* than the original even though the difference
    is exactly `+ε`.  Comparing differences against tolerances sidesteps that —
    and is better conditioned regardless.  (`isapprox` is already difference-based,
    so the `isapprox` predicates above are unaffected.)
    """
    function slerp_unflip_is_the_short_way(R, S, τ)
        d = distance(R, S)
        d > sqrt(eps(basetype(R))) || return true
        distance(R, slerp(R, S, τ; unflip=true)) -
            distance(R, slerp(R, S, τ; unflip=false)) ≤ tol(R, 100)
    end

    # ── Predicates: alignment ─────────────────────────────────────────────────

    """
    `align(a, b)` minimizes `Σ‖aᵢ - ℛbᵢ‖²` — it rotates *b* onto *a*.  So
    rotating a well-conditioned frame by `R` and aligning the rotated frame back
    onto the original must recover `R` itself.
    """
    function align_recovers_rotation(R, a, b, c)
        ε = sqrt(sqrt(eps(basetype(R))))
        # Reject near-degenerate frames, where the problem is ill-conditioned
        # and no tolerance in ulps is meaningful.
        (absvec(a × b) > ε && absvec(b × c) > ε && absvec(a × c) > ε) || return true
        vs = [a, b, c]
        ws = [quatvec(R * v * conj(R)) for v ∈ vs]
        distance(align(ws, vs), R) ≤ sqrt(eps(basetype(R)))
    end

    # ── Predicates: distance ──────────────────────────────────────────────────

    distance_symmetric(R, S) = isapprox(distance(R, S), distance(S, R); atol=tol(R, 100))

    "The rotor metric is bi-invariant: left multiplication is an isometry."
    distance_left_invariant(R, S, T) =
        isapprox(distance(T * R, T * S), distance(R, S); atol=tol(R, 100))

    distance_right_invariant(R, S, T) =
        isapprox(distance(R * T, S * T), distance(R, S); atol=tol(R, 100))

    function distance_triangle(R, S, T)
        slack = distance(R, S) + distance(S, T) - distance(R, T)
        target!(-Float64(slack))       # hunt for the tightest (or violated) case
        slack ≥ -tol(R, 100)
    end
end


# ── Test items ────────────────────────────────────────────────────────────────

@testitem "properties: normed division algebra" tags=[:unit, :fast, :validation] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, CFG, quats, wildfloats,
        associativity, left_distributivity, right_distributivity,
        conj_involution, conj_antiautomorphism, norm_composition, inverse_law

    for T ∈ FloatTypes
        @testset "$T" begin
            Q = quats(T)
            @check config=CFG associativity(Q, Q, Q)
            @check config=CFG left_distributivity(Q, Q, Q)
            @check config=CFG right_distributivity(Q, Q, Q)
            @check config=CFG conj_involution(Q)
            @check config=CFG conj_antiautomorphism(Q, Q)
            @check config=CFG norm_composition(Q, Q)
            @check config=CFG inverse_law(Q)
        end
    end

    # The identities above are scale-free, so they should also survive inputs
    # drawn uniformly over bit patterns — denormals, ±0, and everything between.
    for T ∈ (Float32, Float64)
        @testset "$T (bit-pattern inputs)" begin
            W = quats(wildfloats(T))
            @check config=CFG conj_involution(W)
            @check config=CFG norm_composition(W, W)
        end
    end
end


@testitem "properties: conversion round-trips" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, EigenFloatTypes, CFG, quats, rotors,
        euler_roundtrip, euler_phase_roundtrip, spherical_roundtrip,
        matrix_roundtrip, floatarray_roundtrip

    for T ∈ FloatTypes
        @testset "$T" begin
            R, Q = rotors(T), quats(T)
            @check config=CFG euler_roundtrip(R)
            @check config=CFG euler_phase_roundtrip(R)
            @check config=CFG spherical_roundtrip(R)
            @check config=CFG floatarray_roundtrip(Q)
        end
    end

    # `from_rotation_matrix` needs `eigen`; see EigenFloatTypes.
    for T ∈ EigenFloatTypes
        @testset "$T (via eigen)" begin
            @check config=CFG matrix_roundtrip(rotors(T))
        end
    end
end


@testitem "properties: exp, log, sqrt, powers" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, CFG, quats, rotors,
        explog, logexp_rotor, sqrt_squares, pow_addition

    for T ∈ FloatTypes
        @testset "$T" begin
            Q, R = quats(T), rotors(T)
            @check config=CFG explog(Q)
            @check config=CFG logexp_rotor(R)
            @check config=CFG sqrt_squares(Q)
            @check config=CFG pow_addition(
                R, Data.Integers(-4, 4), Data.Integers(-4, 4))
        end
    end
end


@testitem "properties: rotation semantics" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, CFG, rotors, quatvecs,
        sandwich_matches_matrix, rotation_is_isometry, matrix_homomorphism, double_cover

    for T ∈ FloatTypes
        @testset "$T" begin
            R, V = rotors(T), quatvecs(T)
            @check config=CFG sandwich_matches_matrix(R, V)
            @check config=CFG rotation_is_isometry(R, V)
            @check config=CFG matrix_homomorphism(R, R)
            @check config=CFG double_cover(R, V)
        end
    end
end


@testitem "properties: slerp" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, CFG, rotors, unitinterval,
        slerp_endpoints, slerp_stays_unit, slerp_constant_speed,
        slerp_unflip_is_the_short_way

    for T ∈ FloatTypes
        @testset "$T" begin
            R, τ = rotors(T), unitinterval(T)
            @check config=CFG slerp_endpoints(R, R)
            @check config=CFG slerp_stays_unit(R, R, τ)
            @check config=CFG slerp_constant_speed(R, R, τ)
            @check config=CFG slerp_unflip_is_the_short_way(R, R, τ)
        end
    end
end


@testitem "properties: align recovers a rotation" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: EigenFloatTypes, CFG, rotors, quatvecs,
        align_recovers_rotation

    for T ∈ EigenFloatTypes
        @testset "$T" begin
            R, V = rotors(T), quatvecs(T)
            @check config=CFG align_recovers_rotation(R, V, V, V)
        end
    end
end


@testitem "properties: distance is a bi-invariant metric" tags=[:validation, :fast] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: FloatTypes, CFG, rotors,
        distance_symmetric, distance_left_invariant, distance_right_invariant,
        distance_triangle

    for T ∈ FloatTypes
        @testset "$T" begin
            R = rotors(T)
            @check config=CFG distance_symmetric(R, R)
            @check config=CFG distance_left_invariant(R, R, R)
            @check config=CFG distance_right_invariant(R, R, R)
            @check config=CFG distance_triangle(R, R, R)
        end
    end
end


@testitem "properties: from_rotation_matrix at extended precision" tags=[:validation, :slow] setup=[PropertyGens] begin
    using Supposition: @check, Data
    using .PropertyGens: CFG, rotors, matrix_roundtrip

    # `dominant_eigenvector` (src/conversion.jl) takes `eigen(M).vectors[:, 4]`,
    # which assumes `eigen` returns eigenvalues in ascending order.  LAPACK does,
    # and so does GenericLinearAlgebra on its own — but `GenericSchur` provides a
    # symmetric path that takes precedence and does *not* sort.  The eigenvalues
    # of K₃3 are (-1, -1, -1, 3), so column 4 then holds an eigenvector of the
    # degenerate eigenvalue and `from_rotation_matrix` returns a systematically
    # wrong rotor: `distance` comes out at exactly π/2 for essentially every
    # input, not merely a loss of accuracy.
    #
    # Note the trigger is *loading* GenericSchur, which arrives indirectly with
    # DoubleFloats.  So `from_rotation_matrix(::Matrix{BigFloat})` is correct in a
    # bare session and wrong once the user loads an unrelated package — which also
    # explains the pre-existing BigFloat failures in `conversion.jl`, since the
    # `@testitem`s run (and load DoubleFloats) before that file is included.
    #
    # Selecting by `argmax` of the eigenvalues instead of by position fixes it
    # regardless of load order; flip this to a plain `@check` once that lands.
    @check config=CFG broken=true matrix_roundtrip(rotors(BigFloat))
end
