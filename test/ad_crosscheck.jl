# Cross-backend derivative harness
#
# Every case below is a function from a real vector to a real vector, built from the public
# API of Quaternionic.  Because both the inputs and the outputs are plain real vectors, the
# harness does not depend on how any backend represents tangents of quaternions.  For each
# case and each of its points, the Jacobian computed by an AD backend (through
# DifferentiationInterface) is compared with a reference computed by central finite
# differences of the same function, evaluated in 256-bit `BigFloat` arithmetic with a tiny
# step.  The truncation and rounding errors of that reference are both far below Float64
# roundoff, and it uses no AD backend at all.
#
# The points include generic random points and the special points at which the
# implementations branch or have removable singularities: the identity, vector parts of
# size 1e-4, 1e-8, and 1e-12, an exactly zero vector part, a pure-vector rotor, points near
# -1, equal and antipodal arguments, the ends of interpolation intervals, and the knots of
# `squad`.  Only points at which the function is actually differentiable are included.
#
# Each backend runs in its own test item, so that a crash or a failure in one backend does
# not hide the results of the others.

@testmodule ADCrossCheck begin
    using Test
    using Random
    using LinearAlgebra: LinearAlgebra, norm, normalize, dot
    using Quaternionic
    import Quaternionic: RB, BR, Rv, vR, KAN, ℂreal, ℂimag, ℂconj
    import DifferentiationInterface as DI

    # ── Building blocks ──────────────────────────────────────────────────────────────────

    # Each argument type is built from consecutive entries of the input vector.  A `Rotor`
    # is built with `rotor`, which normalizes its input, so its derivative includes the
    # projection onto the unit sphere.  A `Lorentz` transformation is built as a rotor
    # times a boost with velocity given by three further entries.
    nargs(::Val{:Q}) = 4
    nargs(::Val{:R}) = 4
    nargs(::Val{:V}) = 3
    nargs(::Val{:L}) = 7
    nargs(T::Symbol) = nargs(Val(T))

    Qarg(x, i) = quaternion(x[i], x[i+1], x[i+2], x[i+3])
    Rarg(x, i) = rotor(x[i], x[i+1], x[i+2], x[i+3])
    Varg(x, i) = quatvec(x[i], x[i+1], x[i+2])
    Larg(x, i) = rotor(x[i], x[i+1], x[i+2], x[i+3]) * Boost(quatvec(x[i+4], x[i+5], x[i+6]))

    # Where the argument type is a variable, `make(Val(T), x, i)` selects the constructor.
    # Closures must capture `Val(T)` rather than the symbol `T` itself, and must not pass
    # even a literal symbol to `make`, so that the type of their result can be inferred.
    # Enzyme's reverse mode, in particular, fails on functions whose results are boxed.
    make(::Val{:Q}, x, i) = Qarg(x, i)
    make(::Val{:R}, x, i) = Rarg(x, i)
    make(::Val{:V}, x, i) = Varg(x, i)
    make(::Val{:L}, x, i) = Larg(x, i)

    # The input entries that produce a given (four-component, real) quaternion value.
    embed(::Val{:V}, g) = g[2:4]
    embed(::Val, g) = g
    embed(T::Symbol, g) = embed(Val(T), collect(Float64, g))

    # Flatten any output into a real vector.
    flat(x::Real) = [x]
    flat(z::Complex) = [real(z), imag(z)]
    flat(q::AbstractQuaternion{<:Real}) = collect(components(q))
    flat(q::AbstractQuaternion{<:Complex}) =
        vcat(collect(real.(components(q))), collect(imag.(components(q))))
    flat(a::AbstractArray{<:Real}) = collect(vec(a))
    # Zygote cannot differentiate `vec` of an `SMatrix` (it returns a 3×3 gradient for the
    # 9-vector), which is not the subject of these tests, so matrices are copied first.
    flat(a::AbstractMatrix{<:Real}) = vec(Matrix(a))
    flat(a::AbstractArray{<:Complex}) = vcat(collect(real.(vec(a))), collect(imag.(vec(a))))
    flat(a::AbstractArray) = reduce(vcat, map(flat, vec(a)))
    flat(t::Tuple) = reduce(vcat, map(flat, t))

    # An output that is defined only up to an overall sign (an eigenvector, for example) is
    # replaced by the products of all pairs of its components, which do not depend on that
    # sign and determine the original up to that sign.
    signfree(q::AbstractQuaternion) = (c = flat(q); vec(c * c'))

    # ── Points ───────────────────────────────────────────────────────────────────────────

    const rng = Random.MersenneTwister(20261002)
    const n̂ = normalize([0.36, -0.48, 0.8])
    const m̂ = normalize([-0.6, 0.0, 0.8])
    randq() = randn(rng, 4)
    randv() = randn(rng, 3)

    # Quaternion values (as four-component vectors) with interesting vector parts
    small(s, ε) = [s; ε .* n̂]
    const unitpurevector = [0.0; n̂]

    # Points for a function of one argument of type `T`, as `label => x` pairs, selected by
    # the kind of singularity the function has.  For a `QuatVec`, the kinds are `:smooth`
    # (including the zero vector), `:nonzero`, and `:large` (only points of order 1).
    function points(T::Symbol, kind::Symbol)
        s = T === :Q ? 1.3 : 1.0  # Quaternion points are deliberately not unit-norm
        if T === :V
            generic = ["generic" => randv(), "generic 2" => randv()]
            tiny = ["vec 1e-4" => 1e-4 .* n̂, "vec 1e-8" => 1e-8 .* n̂, "vec 1e-12" => 1e-12 .* n̂]
            zero = ["vec 0" => [0.0, 0.0, 0.0]]
            kind === :generic && return generic[1:1]
            kind === :smooth && return [generic; tiny; zero]
            kind === :nonzero && return [generic; tiny]
            kind === :large && return generic
            error("Unknown kind $kind for $T")
        end
        generic = ["generic" => randq(), "generic 2" => randq()]
        tiny = ["vec 1e-4" => small(s, 1e-4), "vec 1e-8" => small(s, 1e-8),
                "vec 1e-12" => small(s, 1e-12)]
        zero = [(T === :R ? "identity" : "vec 0") => small(s, 0.0)]
        pure = ["pure vector" => [0.0; 1.1 .* n̂]]
        nearm1 = ["near -1 (1e-4)" => small(-s, 1e-4)]
        minus1 = [(T === :R ? "-identity" : "negative real") => small(-s, 0.0)]
        if kind === :generic
            return generic[1:1]
        elseif kind === :smooth  # Differentiable everywhere except, perhaps, at zero
            return [generic; tiny; zero; pure; nearm1; minus1]
        elseif kind === :log     # Not differentiable on the negative real axis
            return [generic; tiny; zero; pure; nearm1]
        elseif kind === :absvec  # Not differentiable where the vector part vanishes
            return [generic; tiny; pure; nearm1]
        end
        error("Unknown kind $kind for $T")
    end

    # ── Cases ────────────────────────────────────────────────────────────────────────────

    # A case is a named function together with the points at which to differentiate it.
    # Second derivatives are not checked at the labels listed in `nohessian`, where the
    # function is only once differentiable.
    const CASES = NamedTuple{(:name, :f, :points, :nohessian),
                             Tuple{String, Function, Vector{Pair{String,Vector{Float64}}}, Vector{String}}}[]

    function addcase!(name, f, pts; nohessian=String[])
        pts = Pair{String,Vector{Float64}}[p.first => collect(Float64, p.second) for p in pts]
        push!(CASES, (name=name, f=f, points=pts, nohessian=collect(String, nohessian)))
        nothing
    end

    # Unary functions of one argument
    function unary!(name, T::Symbol, F, kind::Symbol)
        addcase!("$name($T)", let VT=Val(T); x -> flat(F(make(VT, x, 1))) end, points(T, kind))
    end

    # Concatenate points for several arguments, pairing the i-th points of each list
    function pairup(pts...)
        n = minimum(length, pts)
        [join((p[i].first for p in pts), " | ") => reduce(vcat, (p[i].second for p in pts))
         for i in 1:n]
    end

    # Points for a function of two quaternionic arguments of types `T1` and `T2`
    function binarypoints(T1, T2; equal=true, antipodal=true)
        pts = Pair{String,Vector{Float64}}[]
        g1 = T1 === :V ? randv() : randq()
        g2 = T2 === :V ? randv() : randq()
        push!(pts, "generic" => [g1; g2])
        # A unit pure vector is a valid value for every argument type, so equal and
        # antipodal values can be built for any pair of types.
        if equal
            push!(pts, "equal" => [embed(T1, unitpurevector); embed(T2, unitpurevector)])
            if T1 !== :V && T2 !== :V
                u = normalize(randq())
                push!(pts, "equal 2" => [embed(T1, u); embed(T2, u)])
            end
        end
        if antipodal
            push!(pts, "antipodal" => [embed(T1, unitpurevector); embed(T2, -unitpurevector)])
        end
        pts
    end

    # Constructors and accessors
    #
    # The linear constructors are grouped into a few cases, to save compilation time.
    addcase!(
        "quaternion, quatvec (4 reals)",
        x -> flat((
            quaternion(x[1], x[2], x[3], x[4]), quaternion(x),
            quatvec(x[1], x[2], x[3], x[4]), quatvec(x), quatvec(quaternion(x))
        )),
        points(:Q, :generic)
    )
    addcase!(
        "quaternion, quatvec (3 reals)",
        x -> flat((quaternion(x[1], x[2], x[3]), quatvec(x[1], x[2], x[3]), quatvec(x))),
        points(:V, :smooth)
    )
    addcase!("rotor(w,x,y,z)", x -> flat(rotor(x[1], x[2], x[3], x[4])), points(:R, :smooth))
    addcase!("rotor(x,y,z)", x -> flat(rotor(x[1], x[2], x[3])), points(:V, :large))
    addcase!("rotor(::Vector)", x -> flat(rotor(x)), points(:R, :generic))
    addcase!("rotor(::Quaternion), quaternion(::Rotor)",
             x -> flat(quaternion(rotor(Qarg(x, 1)))), points(:Q, :smooth))
    addcase!(
        "real, imag, vec, getindex, getproperty",
        x -> (q = Qarg(x, 1); flat((real(q), imag(q), vec(q), q[2], q.z))),
        points(:Q, :generic)
    )

    # Algebra of two quaternionic arguments
    const TYPES = (:Q, :R, :V)
    #
    # Multiplication and division are checked for every pair of types, but addition and
    # subtraction, which are simpler, only for some.
    const SUMTYPES = ((:Q, :Q), (:R, :R), (:V, :V), (:Q, :R), (:V, :R))
    for T1 ∈ TYPES, T2 ∈ TYPES, (opname, op) ∈ (("+", +), ("-", -), ("*", *), ("/", /))
        opname ∈ ("+", "-") && (T1, T2) ∉ SUMTYPES && continue
        n1 = nargs(T1)
        addcase!(
            "$T1 $opname $T2",
            let op=op, V1=Val(T1), V2=Val(T2), n1=n1
                x -> flat(op(make(V1, x, 1), make(V2, x, n1 + 1)))
            end,
            binarypoints(T1, T2)
        )
    end

    # Algebra mixing a real scalar `t` (the first input) with a quaternionic argument
    for T ∈ TYPES, (opname, op) ∈ (("+", +), ("-", -), ("*", *), ("/", /))
        pts = Pair{String,Vector{Float64}}[
            "generic" => [0.7; T === :V ? randv() : randq()],
            # The quaternion equals the scalar, where shortcuts based on equality would bite
            "equal" => [1.0; embed(T, [1.0, 0, 0, 0])],
        ]
        T === :V && pop!(pts)  # A QuatVec cannot equal a nonzero scalar
        addcase!("t $opname $T", let op=op, VT=Val(T)
            x -> flat(op(x[1], make(VT, x, 2)))
        end, pts)
        addcase!("$T $opname t", let op=op, VT=Val(T)
            x -> flat(op(make(VT, x, 2), x[1]))
        end, pts)
    end

    # Unary algebra
    for T ∈ TYPES
        unary!("-, conj", T, q -> (-q, conj(q)), :smooth)
        unary!("inv", T, inv, T === :V ? :nonzero : :smooth)
        unary!("abs", T, abs, T === :V ? :nonzero : :smooth)
        unary!("abs2", T, abs2, :smooth)
        unary!("normalize", T, normalize, T === :V ? :nonzero : :smooth)
    end
    for T ∈ (:Q, :R)
        unary!("absvec", T, absvec, :absvec)
        unary!("abs2vec", T, abs2vec, :smooth)
    end

    # Transcendental functions
    unary!("exp", :Q, exp, :smooth)
    unary!("exp", :V, exp, :smooth)
    unary!("log", :Q, log, :log)
    unary!("log", :R, log, :log)
    unary!("sqrt", :Q, sqrt, :log)
    unary!("sqrt", :R, sqrt, :log)
    unary!("sqrt", :V, sqrt, :nonzero)
    unary!("angle", :Q, angle, :absvec)
    unary!("angle", :R, angle, :absvec)

    # Powers with integer exponents
    const EXPONENTS = Dict(:Q => (-2, -1, 0, 1, 2, 3), :R => (-1, 2), :V => (-1, 2, 3))
    for T ∈ TYPES, n ∈ EXPONENTS[T]
        unary!("^$n", T, let n=n; q -> q^n end, T === :V ? (n > 0 ? :smooth : :nonzero) : :smooth)
    end
    # Powers with real exponents, both fixed and as the last input
    for T ∈ TYPES
        kind = T === :V ? :nonzero : :log
        unary!("^0.37", T, q -> q^0.37, kind)
        pts = [p.first => [p.second; 0.37] for p in points(T, kind)]
        addcase!("$T^s", let VT=Val(T); x -> flat(make(VT, x, 1)^x[end]) end, pts)
    end
    addcase!("R^s (large s)", x -> flat(Rarg(x, 1)^x[end]),
             [p.first => [p.second; 7.3] for p in points(:R, :log)])
    addcase!("t^Q", x -> flat(x[1]^Qarg(x, 2)), ["generic" => [1.7; randq()]])
    addcase!("Q^Q", x -> flat(Qarg(x, 1)^Qarg(x, 5)),
             ["generic" => [randq(); randq()], "vec 1e-8 | generic" => [small(1.3, 1e-8); randq()]])

    # Products of vectors
    addcase!("Q ⋅ Q", x -> flat(Qarg(x, 1) ⋅ Qarg(x, 5)), binarypoints(:Q, :Q))
    addcase!("R ⋅ R", x -> flat(Rarg(x, 1) ⋅ Rarg(x, 5)), binarypoints(:R, :R))
    addcase!("V ⋅ V", x -> flat(Varg(x, 1) ⋅ Varg(x, 4)), binarypoints(:V, :V))
    addcase!("Q ⋅ V", x -> flat(Qarg(x, 1) ⋅ Varg(x, 5)), binarypoints(:Q, :V))
    addcase!("V × V", x -> flat(Varg(x, 1) × Varg(x, 4)), binarypoints(:V, :V))
    # The normalized cross product is not differentiable for parallel arguments
    addcase!("V ×̂ V", x -> flat(Varg(x, 1) ×̂ Varg(x, 4)),
             binarypoints(:V, :V; equal=false, antipodal=false))

    # Distances.  `distance` is the square root of `distance2`, so it is not differentiable
    # where `distance2` vanishes: for equal arguments and, for two rotors, for antipodal
    # ones.  For two rotors, `distance2` has a kink where the rotors are orthogonal in ℝ⁴,
    # which no point here approaches.
    for (T1, T2) ∈ ((:Q, :Q), (:R, :R), (:V, :V), (:Q, :R), (:R, :V), (:V, :Q))
        n1 = nargs(T1)
        bothrotors = T1 === :R && T2 === :R
        addcase!("distance2($T1, $T2)", let V1=Val(T1), V2=Val(T2), n1=n1
            x -> flat(distance2(make(V1, x, 1), make(V2, x, n1 + 1)))
        end, binarypoints(T1, T2))
        addcase!("distance($T1, $T2)", let V1=Val(T1), V2=Val(T2), n1=n1
            x -> flat(distance(make(V1, x, 1), make(V2, x, n1 + 1)))
        end, binarypoints(T1, T2; equal=false, antipodal=!bothrotors))
    end
    # Rotor distances near equal and near antipodal arguments, where the series branch is
    # taken
    let
        u = normalize(randq())
        near(ε) = normalize(u .+ ε .* [0.0; n̂])
        pts = ["near equal ($ε)" => [u; near(ε)] for ε ∈ (1e-4, 1e-8, 1e-12)]
        apts = ["near antipodal ($ε)" => [u; -near(ε)] for ε ∈ (1e-4, 1e-8, 1e-12)]
        addcase!("distance2(R, R) near", x -> flat(distance2(Rarg(x, 1), Rarg(x, 5))), [pts; apts])
        addcase!("distance(R, R) near", x -> flat(distance(Rarg(x, 1), Rarg(x, 5))), [pts; apts])
    end

    # Interpolation
    let
        u1, u2 = normalize(randq()), normalize(randq())
        near(u, ε) = normalize(u .+ ε .* [0.0; n̂])
        pts = [
            "generic" => [u1; u2; 0.37],
            "τ = 0" => [u1; u2; 0.0],
            "τ = 1" => [u1; u2; 1.0],
            "q₁ == q₂" => [u1; u1; 0.37],
            "q₁ == q₂, τ = 0" => [u1; u1; 0.0],
            "nearly equal (1e-4)" => [u1; near(u1, 1e-4); 0.37],
            "nearly equal (1e-8)" => [u1; near(u1, 1e-8); 0.37],
            "nearly antipodal (1e-4)" => [u1; -near(u1, 1e-4); 0.37],
            "generic, τ = 1.6" => [u1; u2; 1.6],
            "unnormalized" => [2.1 .* u1; 0.4 .* u2; 0.37],
        ]
        addcase!("slerp", x -> flat(slerp(Rarg(x, 1), Rarg(x, 5), x[9])), pts)
        apts = [
            "generic" => [u1; u2; 0.37],
            "antipodal" => [u1; -u1; 0.37],
            "nearly antipodal (1e-8)" => [u1; -near(u1, 1e-8); 0.37],
        ]
        addcase!("slerp(unflip=true)",
                 x -> flat(slerp(Rarg(x, 1), Rarg(x, 5), x[9]; unflip=true)), apts)
    end
    let
        # Four rotors at the times `tin`; the output time is the last input.  The function
        # is continuously differentiable across the knots, but not twice differentiable.
        tin = [0.0, 0.7, 1.5, 2.6]
        ω⃗ = quatvec(0.3, -0.8, 0.5)
        Rs = [components(exp(t * ω⃗) * rotor(1.0, 0.1, -0.2, 0.05)) for t ∈ tin]
        Rs[3] = Rs[3] .+ 0.01 .* randq()  # Make it a little irregular
        x0 = reduce(vcat, collect.(Rs))
        pts = [
            "interior interval" => [x0; 1.1],
            "first interval" => [x0; 0.3],
            "last interval" => [x0; 2.2],
            "knot t₂" => [x0; 0.7],
            "knot t₃" => [x0; 1.5],
        ]
        addcase!("squad",
                 x -> flat(squad([Rarg(x, 1), Rarg(x, 5), Rarg(x, 9), Rarg(x, 13)],
                                 tin, x[17])),
                 pts; nohessian=["knot t₂", "knot t₃"])
        # The angular velocity and the time derivative are continuous across the knots, but
        # their derivatives with respect to time are not.
        addcase!("squad(compute_angular_velocity=true, compute_derivative=true)",
                 x -> flat(squad([Rarg(x, 1), Rarg(x, 5), Rarg(x, 9), Rarg(x, 13)],
                                 tin, x[17]; compute_angular_velocity=true, compute_derivative=true)),
                 pts[1:3])
    end
    let
        # A sequence of quaternions in which some signs need flipping
        u = [normalize(randq()) for _ ∈ 1:4]
        for i ∈ 2:4
            if dot(u[i], u[i-1]) * (iseven(i) ? 1 : -1) < 0
                u[i] = -u[i]
            end
        end
        pts = ["generic" => reduce(vcat, u)]
        addcase!("unflip(::Vector{Quaternion})",
                 x -> flat(unflip([Qarg(x, 1), Qarg(x, 5), Qarg(x, 9), Qarg(x, 13)])), pts)
        addcase!("unflip(::Vector{Rotor})",
                 x -> flat(unflip([Rarg(x, 1), Rarg(x, 5), Rarg(x, 9), Rarg(x, 13)])), pts)
    end

    # Action of a rotor or a quaternion on a vector
    addcase!("R(v)", x -> flat(Rarg(x, 1)(Varg(x, 5))),
             pairup(points(:R, :smooth), [points(:V, :smooth); points(:V, :smooth)]))
    addcase!("Q(v)", x -> flat(Qarg(x, 1)(Varg(x, 5))),
             pairup(points(:Q, :generic), points(:V, :generic)))

    # Conversions
    let
        # A generic rotor whose Euler angle β is far from 0 and π
        u = [0.6, -0.3, 0.5, 0.4]
        addcase!("to_euler_angles(R)", x -> flat(to_euler_angles(Rarg(x, 1))),
                 ["generic" => u, "β = 1e-4" => collect(components(from_euler_angles(0.3, 1e-4, -2.1)))])
        addcase!("to_euler_angles(Q)", x -> flat(to_euler_angles(Qarg(x, 1))), ["generic" => 1.7u])
        addcase!("from_euler_angles", x -> flat(from_euler_angles(x[1], x[2], x[3])),
                 ["generic" => [0.3, 1.1, -2.1], "β = 0" => [0.3, 0.0, -2.1], "β = π" => [0.3, π, -2.1]])
        addcase!("from_euler_angles(::Vector)", x -> flat(from_euler_angles(x)), ["generic" => [0.3, 1.1, -2.1]])
        addcase!("to_euler_phases", x -> flat(to_euler_phases(Rarg(x, 1))), ["generic" => u])
        addcase!("to_euler_phases!",
                 x -> (R = Rarg(x, 1); flat(to_euler_phases!(Vector{Complex{basetype(R)}}(undef, 3), R))),
                 ["generic" => u])
        z = to_euler_phases(rotor(u...))
        zx = reduce(vcat, [[real(zi), imag(zi)] for zi ∈ z])
        addcase!("from_euler_phases",
                 x -> flat(from_euler_phases(complex(x[1], x[2]), complex(x[3], x[4]), complex(x[5], x[6]))),
                 ["generic" => zx, "off the unit circle" => 1.1 .* zx])
        # The polar angle θ is not differentiable where it vanishes, but it is for small θ.
        sph(θ) = collect(components(from_spherical_coordinates(θ, 0.7)))
        addcase!("to_spherical_coordinates", x -> flat(to_spherical_coordinates(Rarg(x, 1))),
                 ["generic" => u, "θ = 1e-4" => sph(1e-4), "θ = 1e-8" => sph(1e-8)])
        addcase!("from_spherical_coordinates", x -> flat(from_spherical_coordinates(x[1], x[2])),
                 ["generic" => [1.1, -2.1], "θ = 0" => [0.0, -2.1]])
        addcase!("to_rotation_matrix(R)", x -> flat(to_rotation_matrix(Rarg(x, 1))), points(:R, :smooth))
        addcase!("to_rotation_matrix(Q)", x -> flat(to_rotation_matrix(Qarg(x, 1))), points(:Q, :generic))
        # The sign of the result of `from_rotation_matrix` is arbitrary, so only quantities
        # independent of that sign are compared.
        M = vec(to_rotation_matrix(rotor(u...)))
        addcase!("from_rotation_matrix", x -> signfree(from_rotation_matrix(reshape(x, 3, 3))),
                 ["generic" => M, "identity" => vec(Matrix{Float64}(LinearAlgebra.I, 3, 3)),
                  "perturbed" => M .+ 0.01 .* randn(rng, 9)])
        addcase!("to_float_array(::Vector)", x -> flat(to_float_array([Qarg(x, 1), Qarg(x, 5)])),
                 ["generic" => [randq(); randq()]])
        addcase!("to_float_array(q)", x -> flat(to_float_array(Qarg(x, 1))), points(:Q, :generic))
        addcase!("from_float_array", x -> flat(from_float_array(reshape(x, 4, 2))),
                 ["generic" => [randq(); randq()]])
    end

    # Alignment
    let
        R0 = rotor(0.6, -0.3, 0.5, 0.4)
        b = [randv() for _ ∈ 1:3]
        a = [collect(vec(R0(quatvec(bi...)))) .+ 0.05 .* randv() for bi ∈ b]
        w = [0.5, 1.3, 0.9]
        x0 = [reduce(vcat, a); reduce(vcat, b); w]
        V3(x, i) = [Varg(x, i), Varg(x, i + 3), Varg(x, i + 6)]
        # The sign of the dominant eigenvector is arbitrary
        addcase!("align(::Vector{QuatVec}, ::Vector{QuatVec}, w)",
                 x -> signfree(align(V3(x, 1), V3(x, 10), x[19:21])), ["generic" => x0])
        addcase!("align(::Vector{QuatVec}, ::Vector{QuatVec})",
                 x -> signfree(align(V3(x, 1), V3(x, 10))), ["generic" => x0[1:18]])
        B = [normalize(randq()) for _ ∈ 1:3]
        A = [collect(components(R0 * rotor(Bi...))) .+ 0.05 .* randq() for Bi ∈ B]
        y0 = [reduce(vcat, A); reduce(vcat, B); w]
        R3(x, i) = [Rarg(x, i), Rarg(x, i + 4), Rarg(x, i + 8)]
        addcase!("align(::Vector{Rotor}, ::Vector{Rotor}, w)",
                 x -> flat(align(R3(x, 1), R3(x, 13), x[25:27])), ["generic" => y0])
        addcase!("align(::Vector{Rotor}, ::Vector{Rotor})",
                 x -> flat(align(R3(x, 1), R3(x, 13))), ["generic" => y0[1:24]])
    end

    # Lorentz transformations
    let
        vpts = ["generic" => [0.3, -0.2, 0.5], "zero velocity" => [0.0, 0.0, 0.0],
                "velocity 1e-4" => 1e-4 .* n̂, "velocity 1e-8" => 1e-8 .* n̂]
        addcase!("Boost(η, ::Vector)", x -> flat(Boost(x[1], x[2:4])),
                 ["generic" => [0.8; n̂], "η = 0" => [0.0; n̂], "unnormalized" => [0.8; 1.3 .* m̂]])
        addcase!("Boost(constant η, ::Vector)", x -> flat(Boost(0.8, x)), ["generic" => n̂])
        addcase!("Boost(η, ::QuatVec)", x -> flat(Boost(x[1], Varg(x, 2))),
                 ["generic" => [0.8; n̂], "η = 0" => [0.0; n̂]])
        addcase!("Boost(::QuatVec)", x -> flat(Boost(Varg(x, 1))), vpts)
        addcase!("Boost(::Vector)", x -> flat(Boost(x[1:3])), vpts)
        lpts = [p.first => [normalize(randq()); p.second] for p in vpts]
        addcase!("ga_components(Λ)", x -> flat(ga_components(Larg(x, 1))), lpts)
        addcase!("Λ * Λ", x -> flat(Larg(x, 1) * Larg(x, 8)), pairup(lpts, reverse(lpts)))
        addcase!("inv(Λ)", x -> flat(inv(Larg(x, 1))), lpts)
        addcase!("conj(Λ)", x -> flat(conj(Larg(x, 1))), lpts)
        addcase!("Λ * R", x -> flat(Larg(x, 1) * Rarg(x, 8)), [p.first => [p.second; randq()] for p in lpts])
        addcase!("Λ(v)", x -> flat(Larg(x, 1)(x[8:11])), [p.first => [p.second; randq()] for p in lpts])
        addcase!("ℂconj, ℂreal, ℂimag", x -> (Λ = Larg(x, 1); flat((ℂconj(Λ), ℂreal(Λ), ℂimag(Λ)))), lpts)
        addcase!("RB", x -> flat(RB(Larg(x, 1))), lpts)
        addcase!("BR", x -> flat(BR(Larg(x, 1))), lpts)
        addcase!("Rv", x -> flat(Rv(Larg(x, 1))), lpts)
        addcase!("vR", x -> flat(vR(Larg(x, 1))), lpts)
        addcase!("KAN", x -> flat(KAN(Larg(x, 1))), lpts)
    end

    # ── References ───────────────────────────────────────────────────────────────────────

    const PRECISION = 256

    # Central differences in 256-bit arithmetic.  With h = 1e-30, the truncation error is
    # of order h² f‴ ≈ 1e-60 f‴, and the rounding error is of order 2⁻²⁵⁶/h ≈ 1e-47.
    function reference_jacobian(f, x::Vector{BigFloat})
        h = BigFloat(10)^-30
        columns = map(eachindex(x)) do j
            e = zeros(BigFloat, length(x))
            e[j] = h
            (f(x .+ e) .- f(x .- e)) ./ (2h)
        end
        reduce(hcat, columns)
    end

    # The Hessian of `s` by central second differences in 256-bit arithmetic.  With
    # h = 1e-25, the truncation error is of order h² s⁗ ≈ 1e-50 s⁗, and the rounding error
    # is of order 2⁻²⁵⁶/h² ≈ 1e-27.
    function reference_hessian(s, x::Vector{BigFloat})
        h = BigFloat(10)^-25
        n = length(x)
        e(i) = (v = zeros(BigFloat, n); v[i] = h; v)
        H = zeros(BigFloat, n, n)
        for i ∈ 1:n, j ∈ i:n
            H[i, j] = H[j, i] = (
                s(x .+ e(i) .+ e(j)) - s(x .+ e(i) .- e(j))
                - s(x .- e(i) .+ e(j)) + s(x .- e(i) .- e(j))
            ) / (4h^2)
        end
        H
    end

    # Even an implementation that is backward stable — one that returns the exact result
    # for inputs perturbed by a few ulps — cannot have derivatives more accurate than the
    # change in the exact derivative under such a perturbation.  Where the derivative is
    # ill-conditioned, for example in the direction of a difference of two nearly equal
    # rotors, that change can be far larger than a few hundred ulps of the derivative
    # itself.  So, alongside each reference, we compute the change in the reference when
    # every input is perturbed by at most one ulp in relative terms, and the tolerance
    # allows `CONDITIONING` times that change.
    const CONDITIONING = 30

    function reference(derivative, f, x::Vector{Float64}, seed)
        setprecision(BigFloat, PRECISION) do
            xb = BigFloat.(x)
            r = 2 .* rand(Random.MersenneTwister(seed), length(x)) .- 1
            x̃ = xb .* (1 .+ eps(Float64) .* BigFloat.(r))
            D = derivative(f, xb)
            D̃ = derivative(f, x̃)
            Float64.(D), Float64(norm(D̃ .- D))
        end
    end

    # The references depend only on the case and the point, so they are computed once per
    # process and shared by all backends.
    const JACOBIANS = Dict{Tuple{String,String},Tuple{Matrix{Float64},Float64}}()
    const HESSIANS = Dict{Tuple{String,String},Tuple{Matrix{Float64},Float64}}()
    reference_jacobian(c, p::Pair) = get!(JACOBIANS, (c.name, p.first)) do
        reference(reference_jacobian, c.f, p.second, hash((c.name, p.first)))
    end
    reference_hessian(c, p::Pair) = get!(HESSIANS, (c.name, p.first)) do
        reference(reference_hessian, scalarize(c), p.second, hash((c.name, p.first)))
    end

    # A scalar function for the second-derivative checks: a fixed random linear combination
    # of the outputs of the case
    const WEIGHTS = Dict{String,Vector{Float64}}()
    function scalarize(c)
        m = length(c.f(c.points[1].second))
        w = get!(() -> 0.5 .+ rand(Random.MersenneTwister(hash(c.name)), m), WEIGHTS, c.name)
        x -> dot(w, c.f(x))
    end

    # Errors are measured relative to the size of the reference, but never relative to a
    # size smaller than 1, so that a vanishing reference does not demand an exact zero.
    scale(Dref) = max(norm(Dref), 1)
    relerr(D, Dref) = norm(D .- Dref) / scale(Dref)

    # ── Running the checks ───────────────────────────────────────────────────────────────

    # Every outcome is also recorded here, so that a script can tabulate the results.  If
    # the environment variable `QUATERNIONIC_AD_CROSSCHECK_RESULTS` names a file, each
    # outcome is also appended to that file as a line of tab-separated values: the backend,
    # the case, the point, the status, and the details.
    const RESULTS = NamedTuple{(:backend, :case, :point, :status, :detail),
                               Tuple{String,String,String,Symbol,String}}[]

    function record!(backend, c, p, status, detail)
        push!(RESULTS, (backend=backend, case=c.name, point=p.first, status=status, detail=detail))
        file = get(ENV, "QUATERNIONIC_AD_CROSSCHECK_RESULTS", "")
        if !isempty(file)
            open(file, "a") do io
                fields = (backend, c.name, p.first, string(status), detail)
                println(io, join((replace(f, '\t' => ' ', '\n' => ' ') for f ∈ fields), '\t'))
            end
        end
        nothing
    end

    firstline(e) = first(split(sprint(showerror, e), '\n'))

    # Cases that a backend cannot differentiate, because of a bug or a limitation of the
    # backend itself rather than of Quaternionic, with the reason.  Their checks are
    # `@test_broken`, so that a fix in the backend shows up as an unexpected pass, which
    # should then be removed from this list.  Each reason names a reproduction that does
    # not need this harness; the reproductions are kept with the notes of the 4.4.5 patch,
    # in `scratchpad/upstream/`.
    const BROKEN = Dict{Tuple{String,String},String}()
    for case ∈ ("squad", "squad(compute_angular_velocity=true, compute_derivative=true)",
                "unflip(::Vector{Quaternion})", "unflip(::Vector{Rotor})", "to_euler_phases!")
        BROKEN[("Zygote", case)] =
            "Zygote does not support mutating arrays, which this function does " *
            "(reproduction: upstream/zygote_mutation.jl)"
    end
    BROKEN[("Zygote", "from_euler_angles(::Vector)")] =
        "Zygote's `jacobian` fails for a function that splats its vector argument " *
        "(reproduction: upstream/zygote_jacobian_splat.jl)"
    for case ∈ ("to_float_array(::Vector)", "from_float_array")
        BROKEN[("Zygote", case)] =
            "Zygote has no adjoint for the constructor of a `ReinterpretArray` " *
            "(reproduction: upstream/zygote_reinterpret.jl)"
    end
    for case ∈ ("squad", "squad(compute_angular_velocity=true, compute_derivative=true)",
                "unflip(::Vector{Quaternion})", "unflip(::Vector{Rotor})")
        BROKEN[("Enzyme reverse", case)] =
            "Enzyme's batched reverse mode crashes the compiler on this function " *
            "(reproductions: upstream/enzyme_3.jl and upstream/enzyme_4.jl)"
    end
    isbroken(label, c) = haskey(BROKEN, (label, c.name))

    # An optional hook, called before each evaluation, so that a script can find out which
    # case was running when a backend crashed the process
    const BEFORE = Ref{Any}(nothing)

    function check(label, derivative, reference, rtol, c, p)
        BEFORE[] === nothing || BEFORE[](label, c, p)
        if isbroken(label, c)
            ok = try
                D = derivative(c, p)
                Dref, κ = reference(c, p)
                err, tol = relerr(D, Dref), rtol + CONDITIONING * κ / scale(Dref)
                record!(label, c, p, err ≤ tol ? :pass : :broken, "relerr=$(err) tol=$(tol)")
                err ≤ tol
            catch e
                record!(label, c, p, :broken, firstline(e))
                false
            end
            @test_broken ok
            return
        end
        D = try
            derivative(c, p)
        catch e
            record!(label, c, p, :error, firstline(e))
            rethrow()
        end
        err, tol = try
            Dref, κ = reference(c, p)
            relerr(D, Dref), rtol + CONDITIONING * κ / scale(Dref)
        catch e
            record!(label, c, p, :error, "comparison failed: " * firstline(e))
            rethrow()
        end
        record!(label, c, p, err ≤ tol ? :pass : :fail, "relerr=$(err) tol=$(tol)")
        @test err ≤ tol
    end

    """
        crosscheck(label, backend; rtol=300eps(), cases=CASES)

    Compare the Jacobian of every case in `cases` at every point, computed by `backend`,
    with the reference.  The error relative to the size of the reference may be at most
    `rtol`, plus `CONDITIONING` times the change in the reference under a one-ulp
    perturbation of the input, also relative to the size of the reference.  The outcome of
    each comparison is also recorded in `RESULTS`.
    """
    function crosscheck(label, backend; rtol=300eps(), cases=CASES)
        derivative(c, p) = DI.jacobian(c.f, backend, p.second)
        @testset "$label" begin
            for c ∈ cases
                @testset "$(c.name)" begin
                    for p ∈ c.points
                        @testset "$(p.first)" begin
                            check(label, derivative, reference_jacobian, rtol, c, p)
                        end
                    end
                end
            end
        end
    end

    # Second derivatives are checked for every case with a scalar output, and for the cases
    # with vector outputs whose second derivatives are at particular risk: those with
    # series expansions near removable singularities, and those that go through an
    # eigendecomposition.  Checking every case would take much longer.
    const HESSIAN_NAMES = [
        "rotor(w,x,y,z)", "Q * Q", "R * R", "Q / Q", "R / R", "V / V", "inv(Q)", "normalize(V)",
        "exp(Q)", "exp(V)", "log(Q)", "log(R)", "sqrt(Q)", "sqrt(R)",
        "^-1(Q)", "^2(Q)", "^0.37(Q)", "Q^s", "R^s", "V^s", "Q^Q",
        "slerp", "squad", "R(v)", "to_euler_angles(R)", "from_euler_phases",
        "to_spherical_coordinates", "from_rotation_matrix",
        "align(::Vector{QuatVec}, ::Vector{QuatVec}, w)", "align(::Vector{Rotor}, ::Vector{Rotor}, w)",
        "Boost(::QuatVec)", "Λ * Λ", "Λ(v)", "RB", "Rv", "KAN",
    ]
    hessian_cases() =
        filter(c -> c.name ∈ HESSIAN_NAMES || length(c.f(c.points[1].second)) == 1, CASES)

    """
        crosscheck_hessian(label, backend; rtol=1e4eps(), cases=hessian_cases())

    Compare the Hessian of a fixed random linear combination of the outputs of each case at
    every point at which it is twice differentiable, computed by `backend`, with the
    reference.  The tolerance is formed as in `crosscheck`.
    """
    function crosscheck_hessian(label, backend; rtol=1e4eps(), cases=hessian_cases())
        derivative(c, p) = DI.hessian(scalarize(c), backend, p.second)
        @testset "$label" begin
            for c ∈ cases
                @testset "$(c.name)" begin
                    for p ∈ c.points
                        p.first ∈ c.nohessian && continue
                        @testset "$(p.first)" begin
                            check(label, derivative, reference_hessian, rtol, c, p)
                        end
                    end
                end
            end
        end
    end

    # ── Running in a separate process ────────────────────────────────────────────────────

    # Some backends can crash the whole Julia process (with a segmentation fault inside the
    # Enzyme compiler, for example) rather than throw an error.  When the test items all run
    # in one process, as under `runtests.jl`, such a crash would end the whole test run.  So
    # these backends run in a child process, which is restarted without the case that was
    # running when it crashed.  The crash is then reported as an error for that case.

    # A test set that records nothing, for use in the child process, which reports its
    # outcomes through a file instead
    struct QuietTestSet <: Test.AbstractTestSet end
    QuietTestSet(description; kwargs...) = QuietTestSet()
    Test.record(::QuietTestSet, result) = nothing
    Test.finish(::QuietTestSet) = nothing

    # The body of the child process: run the cases not listed in `skipfile`, appending each
    # outcome to `resultsfile` and the name of each case and point to `progressfile` just
    # before it is evaluated.
    function child(label, backend, skipfile, progressfile, resultsfile; rtol)
        skip = Set(filter(!isempty, readlines(skipfile)))
        written = Ref(0)
        function flushresults()
            open(resultsfile, "a") do io
                for r ∈ RESULTS[written[]+1:end]
                    detail = replace(r.detail, '\t' => ' ', '\n' => ' ')
                    println(io, join((r.case, r.point, r.status, detail), '\t'))
                end
            end
            written[] = length(RESULTS)
        end
        BEFORE[] = (label, c, p) -> begin
            flushresults()
            write(progressfile, string(c.name, '\t', p.first))
        end
        Test.@testset QuietTestSet "child" begin
            crosscheck(label, backend; rtol=rtol, cases=filter(c -> c.name ∉ skip, CASES))
        end
        flushresults()
        write(progressfile, "DONE")
        nothing
    end

    """
        crosscheck_isolated(label, imports, backend; rtol=300eps(), cases=CASES, maxrestarts=20)

    Run `crosscheck` in a child process, and report its outcomes as tests in this process.
    Here, `imports` is code that loads the backend in the child, and `backend` is code that
    constructs the backend there.  If the child crashes, the case that was running is
    reported as an error, and a new child continues with the remaining cases.
    """
    function crosscheck_isolated(label, imports, backend; rtol=300eps(), cases=CASES, maxrestarts=20)
        self = joinpath(pkgdir(Quaternionic), "test", "ad_crosscheck.jl")
        dir = mktempdir()
        skipfile, progressfile, resultsfile, logfile =
            joinpath.(dir, ("skip", "progress", "results", "log"))
        script = """
            using Test
            import DifferentiationInterface as DI
            $imports
            let ex = Meta.parseall(read($(repr(self)), String); filename=$(repr(self)))
                m = only(a for a ∈ ex.args
                         if a isa Expr && a.head === :macrocall && a.args[1] === Symbol("@testmodule"))
                Core.eval(Main, Expr(:module, true, m.args[3], m.args[4]))
            end
            Main.ADCrossCheck.child($(repr(label)), $backend, $(repr(skipfile)),
                                    $(repr(progressfile)), $(repr(resultsfile)); rtol=$rtol)
            """
        # The child reports its outcomes to this process, which records them.
        env = delete!(copy(ENV), "QUATERNIONIC_AD_CROSSCHECK_RESULTS")
        cmd = setenv(
            `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $script`,
            env
        )
        crashes = Dict{Tuple{String,String},String}()
        skip = Set{String}(c.name for c ∈ CASES if c ∉ cases)
        for _ ∈ 1:maxrestarts
            write(skipfile, join(skip, '\n'))
            rm(progressfile; force=true)
            process = run(pipeline(ignorestatus(cmd); stdout=logfile, stderr=logfile))
            progress = isfile(progressfile) ? read(progressfile, String) : ""
            progress == "DONE" && break
            log = read(logfile, String)
            if isempty(progress)
                # The child failed before reaching any case, so restarting would not help.
                error("The child process for $label failed before running any case:\n" *
                      join(last(split(log, '\n'), 30), '\n'))
            end
            case, point = split(progress, '\t')
            signal = match(r"signal \(?\d+\)?[^\n]*|LLVM ERROR[^\n]*|Assertion[^\n]*", log)
            how = process.termsignal == 0 ? "exit code $(process.exitcode)" : "signal $(process.termsignal)"
            crashes[(case, point)] = "The process crashed ($how)" *
                (signal === nothing ? "" : ": $(signal.match)")
            push!(skip, case)
            # Cases already completed are not run again.
            isfile(resultsfile) &&
                foreach(l -> push!(skip, first(split(l, '\t'))), readlines(resultsfile))
        end
        outcomes = Dict{Tuple{String,String},Tuple{String,String}}()
        if isfile(resultsfile)
            for l ∈ readlines(resultsfile)
                case, point, status, detail = split(l, '\t'; limit=4)
                outcomes[(case, point)] = (status, detail)
            end
        end
        @testset "$label" begin
            for c ∈ cases
                @testset "$(c.name)" begin
                    for p ∈ c.points
                        @testset "$(p.first)" begin
                            key = (c.name, p.first)
                            if isbroken(label, c)
                                status, detail = get(outcomes, key, ("broken",
                                    get(crashes, key, "Not reached by the child process")))
                                record!(label, c, p, Symbol(status), detail)
                                @test_broken status == "pass"
                            elseif haskey(crashes, key)
                                record!(label, c, p, :crash, crashes[key])
                                error(crashes[key])
                            elseif !haskey(outcomes, key)
                                message = any(k -> k[1] == c.name, keys(crashes)) ?
                                    "Not reached, because another point of this case crashed the process" :
                                    "Not reached by the child process"
                                record!(label, c, p, :error, message)
                                error(message)
                            else
                                status, detail = outcomes[key]
                                record!(label, c, p, Symbol(status), detail)
                                status == "error" && error(detail)
                                @test status == "pass"
                            end
                        end
                    end
                end
            end
        end
    end
end


@testitem "AD cross-check: ForwardDiff" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    import DifferentiationInterface as DI
    import ForwardDiff
    ADCrossCheck.crosscheck("ForwardDiff", DI.AutoForwardDiff())
end

@testitem "AD cross-check: ForwardDiff second derivatives" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    import DifferentiationInterface as DI
    import ForwardDiff
    # With the default chunk size, nested dual numbers with as many as 12 partials at each
    # level make compilation take a minute or more for some cases.  A chunk size of 1 computes the same
    # Hessian in a fraction of a second.
    ADCrossCheck.crosscheck_hessian("ForwardDiff (Hessian)", DI.AutoForwardDiff(chunksize=1))
end

@testitem "AD cross-check: Zygote" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    import DifferentiationInterface as DI
    import Zygote
    ADCrossCheck.crosscheck("Zygote", DI.AutoZygote())
end

@testitem "AD cross-check: ReverseDiff" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    import DifferentiationInterface as DI
    import ReverseDiff
    ADCrossCheck.crosscheck("ReverseDiff", DI.AutoReverseDiff())
end

@testitem "AD cross-check: Mooncake" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    import DifferentiationInterface as DI
    import Mooncake
    ADCrossCheck.crosscheck("Mooncake", DI.AutoMooncake(config=nothing))
end

@testitem "AD cross-check: Enzyme forward" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    # Enzyme can crash the process, so the cases run in a child process.
    ADCrossCheck.crosscheck_isolated(
        "Enzyme forward",
        "import Enzyme",
        "DI.AutoEnzyme(mode=Enzyme.Forward, function_annotation=Enzyme.Const)"
    )
end

@testitem "AD cross-check: Enzyme reverse" setup=[ADCrossCheck] tags=[:slow, :validation] begin
    # Enzyme can crash the process, so the cases run in a child process.
    ADCrossCheck.crosscheck_isolated(
        "Enzyme reverse",
        "import Enzyme",
        "DI.AutoEnzyme(mode=Enzyme.Reverse, function_annotation=Enzyme.Const)"
    )
end
