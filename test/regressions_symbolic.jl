# Regression tests for fixes to the Symbolics, FastDifferentiation, and Latexify extensions.
#
# Numerical results are compared with references computed independently of the extensions:
# explicit formulas evaluated in BigFloat, or central finite differences in BigFloat.

@testmodule SymbolicReference begin
    using Test

    # The Jacobian of the vector-valued function `f` at the point `x`, computed by central
    # finite differences in 256-bit BigFloat, and rounded to Float64
    function bigfloat_jacobian(f, x)
        setprecision(BigFloat, 256) do
            h = big"1e-30"
            xb = big.(x)
            n = length(xb)
            columns = [
                (f(xb .+ h .* ((1:n) .== j)) .- f(xb .- h .* ((1:n) .== j))) ./ (2h)
                for j ∈ 1:n
            ]
            Float64.(reduce(hcat, columns))
        end
    end

    # The ambiguities among the methods of `checked_modules` that involve a method defined
    # in any of `modules`.  The methods in `src/matrices.jl` are excluded: they only throw,
    # for operations such as `transpose(x) * A` that would be wrong for quaternions, and
    # they are necessarily ambiguous with the methods that other packages, Symbolics among
    # them, define for those operations with their own types of arrays.  An ambiguous call
    # throws a `MethodError` instead of the `ArgumentError`, and Aqua checks these methods
    # for ambiguities with Base, LinearAlgebra, and StaticArrays.
    function ambiguities_involving(modules, checked_modules...)
        ambiguities = Test.detect_ambiguities(checked_modules...; recursive=true)
        isguard(m) = endswith(String(m.file), joinpath("src", "matrices.jl"))
        filter(ambiguities) do (m1, m2)
            (parentmodule(m1) ∈ modules || parentmodule(m2) ∈ modules) &&
                !isguard(m1) && !isguard(m2)
        end
    end
end

@testitem "Symbolics: no method ambiguities involving the extension" tags=[:unit] setup=[SymbolicReference] begin
    import Symbolics, FastDifferentiation, Latexify
    extensions = [
        Base.get_extension(Quaternionic, name)
        for name ∈ (
            :QuaternionicSymbolicsExt,
            :QuaternionicFastDifferentiationExt,
            :QuaternionicLatexifyExt,
        )
    ]
    modules = (extensions..., Symbolics, FastDifferentiation, Latexify)
    ambiguities = SymbolicReference.ambiguities_involving(modules, Quaternionic, extensions...)
    @test isempty(ambiguities)
end

@testitem "Symbolics: mixed == with numeric values" tags=[:unit] begin
    import Symbolics
    Symbolics.@variables a b c d
    qs = quaternion(a, b, c, d)
    qf = quaternion(1.0, 2.0, 3.0, 4.0)

    # These comparisons used to throw ambiguity errors.
    @test (a == qf) == false
    @test (qf == a) == false
    @test (𝐢 == a) == false
    @test (a == 𝐢) == false
    @test (qf == qs) == false
    @test (quatvec(qf) == qs) == false
    @test (rotor(qf) == qs) == false
    @test (qs == qf) == false

    # Numerically equal values compare equal in either order.
    qn = quaternion(Symbolics.Num(1.0), Symbolics.Num(2), Symbolics.Num(3), Symbolics.Num(4))
    @test qf == qn
    @test qn == qf
    @test quatvec(qf) == quatvec(qn)
    @test quatvec(qn) == quatvec(qf)
    @test quatvec(qn) == quaternion(0.0, 2, 3, 4)
    @test quaternion(0.0, 2, 3, 4) == quatvec(qn)
    @test rotor(qf) == Rotor{Symbolics.Num}(components(rotor(qf)))
    @test Rotor{Symbolics.Num}(components(rotor(qf))) == rotor(qf)
    @test Symbolics.Num(1) == quaternion(1.0)
    @test quaternion(1.0) == Symbolics.Num(1)
    @test quaternion((a + b)^2) == a^2 + 2a*b + b^2
    @test a^2 + 2a*b + b^2 == quaternion((a + b)^2)
end

@testitem "Symbolics: QuatVec equals a scalar exactly when both are zero" tags=[:unit] begin
    import Symbolics
    Symbolics.@variables a b c
    v = quatvec(a, b, c)
    z = v - v
    @test z isa QuatVec{Symbolics.Num}
    @test z == 0
    @test 0 == z
    @test z == 0.0
    @test z == Symbolics.Num(0)
    @test Symbolics.Num(0) == z
    @test (z == 1) == false
    @test (1 == z) == false
    @test (v == 0) == false
    @test (0 == v) == false
    @test (v == a) == false
    @test (a == v) == false
    # A numeric QuatVec compared with a symbolic zero
    @test zero(QuatVecF64) == Symbolics.Num(0)
    @test Symbolics.Num(0) == zero(QuatVecF64)
    @test (𝐢 == Symbolics.Num(0)) == false
end

@testitem "Symbolics: promotion with QuatVec and Rotor gives Quaternion" tags=[:unit] begin
    import Symbolics
    Num = Symbolics.Num
    Symbolics.@variables a
    @test promote_type(QuatVecF64, Num) === Quaternion{Num}
    @test promote_type(Num, QuatVecF64) === Quaternion{Num}
    @test promote_type(QuatVec{Bool}, Num) === Quaternion{Num}
    @test promote_type(RotorF64, Num) === Quaternion{Num}
    @test promote_type(Num, RotorF64) === Quaternion{Num}
    @test promote_type(QuaternionF64, Num) === Quaternion{Num}
    @test promote_type(QuatVec{Num}, Num) === Quaternion{Num}
    @test promote_type(Rotor{Num}, Num) === Quaternion{Num}
    @test promote_type(Quaternion{Num}, Num) === Quaternion{Num}
    @test [imx, a] isa Vector{Quaternion{Num}}
    @test [quatvec(1.0, 2, 3), a] isa Vector{Quaternion{Num}}
    @test promote(quatvec(1.0, 2, 3), a) isa Tuple{Quaternion{Num}, Quaternion{Num}}
    @test promote(rotor(1.0, 2, 3, 4), a) isa Tuple{Quaternion{Num}, Quaternion{Num}}
end

@testitem "Symbolics: dividing a symbolic quaternion by a numeric one" tags=[:unit] begin
    import Symbolics
    Symbolics.@variables a b c d
    # These values have unit norm to Float64 precision, so that they can be substituted
    # into a symbolic `Rotor` built without symbolic normalization.
    values = Tuple(components(rotor(0.1, -0.2, 0.35, 0.9)))
    # Each component is evaluated numerically by a compiled function, which does not depend
    # on Symbolics folding constant subexpressions such as `sqrt(1.0)`.
    evaluate(x) = Symbolics.build_function(x, a, b, c, d; expression=Val(false))(values...)
    symbolic = (quaternion(a, b, c, d), Rotor{Symbolics.Num}(quaternion(a, b, c, d)))
    for p ∈ (rotor(0.3, -0.4, 0.5, 0.7), quaternion(0.3, -0.4, 0.5, 0.7))
        for qs ∈ symbolic
            r = qs / p
            @test r isa (p isa Rotor && qs isa Rotor ? Rotor : Quaternion){Symbolics.Num}
            # Reference: q / p = q p̄ / |p|², from explicit components in BigFloat
            q1, q2, q3, q4 = big.(values)
            p1, p2, p3, p4 = big.(components(p))
            n = p1^2 + p2^2 + p3^2 + p4^2
            expected = [
                (q1*p1 + q2*p2 + q3*p3 + q4*p4) / n,
                (-q1*p2 + q2*p1 - q3*p4 + q4*p3) / n,
                (-q1*p3 + q2*p4 + q3*p1 - q4*p2) / n,
                (-q1*p4 - q2*p3 + q3*p2 + q4*p1) / n,
            ]
            computed = [evaluate(x) for x ∈ components(r)]
            @test maximum(abs, computed .- expected) < 10eps(Float64)
        end
    end
end

@testitem "Symbolics: scalar constructors" tags=[:unit] begin
    import Symbolics
    Num = Symbolics.Num
    Symbolics.@variables a

    # Zero components print as `0`, not as `false`.
    @test !occursin("false", string(quaternion(a)))
    @test !occursin("false", string(quatvec(a)))
    @test !occursin("false", string(rotor(a)))
    @test quaternion(a) == quaternion(a, 0, 0, 0)
    @test quatvec(a) == 0

    # The rotor of a numerical constant has the same components as the rotor of that
    # constant as a real number: it has the sign of a finite nonzero constant, and it is a
    # NaN rotor for a zero, NaN, or infinite constant.
    @test rotor(Num(-2)) isa Rotor{Num}
    @test Symbolics.value(rotor(Num(-2))[1]) == -1
    @test Symbolics.value(rotor(Num(2.5))[1]) == 1
    @test rotor(Num(-2)) == rotor(-2.0)
    @test rotor(Num(2.5)) == rotor(2.5)
    @test string(rotor(Num(-2))) == string(rotor(-2))
    @test all(
        isequal(
            Float64.(Symbolics.value.(components(rotor(Num(x))))),
            Float64.(components(rotor(x)))
        )
        for x ∈ (0, false, 0.0, -0.0, NaN, Inf, -Inf, -2, 2.5, 3//2)
    )
    @test all(isnan, Symbolics.value.(components(rotor(Num(0)))))
    @test all(isnan, Symbolics.value.(components(rotor(Num(NaN)))))
    # For a free symbol, the sign is unknown, and the identity rotor is returned.
    @test rotor(a) == one(QuaternionF64)
    @test string(rotor(a)) == "rotor(1 + 0𝐢 + 0𝐣 + 0𝐤)"
end

@testitem "Symbolics: mixed operations used in the precompile workload" tags=[:unit] begin
    import Symbolics
    Symbolics.@variables w x y z a b c d
    values = (
        0.7, quatvec(0.1, -0.2, 0.3), rotor(0.3, -0.4, 0.5, 0.7), quaternion(1.0, 2.0, -3.0, 0.5),
        w, quatvec(x, y, z), rotor(a, b, c, d), quaternion(w, x, y, z)
    )
    succeeds(f) = try
        f()
        true
    catch
        false
    end
    nfailures = count(
        !succeeds(() -> (conj(p); p * q; p / q; p + q; p - q))
        for p ∈ values for q ∈ values
    )
    @test nfailures == 0
end

@testitem "FastDifferentiation: unsupported functions raise their messages" tags=[:unit] begin
    import FastDifferentiation
    import Quaternionic: to_euler_phases
    u = FastDifferentiation.make_variables(:u, 4)
    Q = quaternion(u...)
    R = rotor(u...)
    V = quatvec(u[2], u[3], u[4])
    for (name, f) ∈ (
        ("Base.exp", () -> exp(Q)),
        ("Base.exp", () -> exp(V)),
        ("Base.log", () -> log(Q)),
        ("Base.log", () -> log(R)),
        ("Base.sqrt", () -> sqrt(Q)),
        ("Base.sqrt", () -> sqrt(V)),
        ("Base.sqrt", () -> sqrt(R)),
        ("Base.:^", () -> Q^2),
        ("Base.:^", () -> V^2),
        ("Base.:^", () -> R^2),
        ("Base.:^", () -> R^0.5),
        ("Base.:^", () -> R^(1//2)),
        ("Base.log", () -> Q^1.5),
        ("Quaternionic.to_euler_phases", () -> to_euler_phases(R)),
    )
        err = try
            f()
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test err isa ErrorException && occursin("`$(name)` cannot yet be used", err.msg)
    end
end

@testitem "FastDifferentiation: scalar constructors and promotion" tags=[:unit] begin
    import FastDifferentiation
    Node = FastDifferentiation.Node
    u = FastDifferentiation.make_variables(:u, 1)
    @test !occursin("false", string(quaternion(u[1])))
    @test !occursin("false", string(quatvec(u[1])))
    @test !occursin("false", string(rotor(u[1])))
    @test string(rotor(Node(-2))) == string(rotor(-2))
    # As for Symbolics, the rotor of a constant matches the rotor of that real number.
    @test all(
        isequal(
            Float64.(FastDifferentiation.value.(components(rotor(Node(x))))),
            Float64.(components(rotor(x)))
        )
        for x ∈ (0, false, 0.0, -0.0, NaN, Inf, -Inf, -2, 2.5, 3//2)
    )
    @test all(isnan, FastDifferentiation.value.(components(rotor(Node(0)))))
    @test string(rotor(u[1])) == "rotor(1 + 0𝐢 + 0𝐣 + 0𝐤)"
    @test promote_type(QuatVecF64, Node) === Quaternion{Node}
    @test promote_type(Node, QuatVecF64) === Quaternion{Node}
    @test promote_type(RotorF64, Node) === Quaternion{Node}
    @test promote_type(QuaternionF64, Node) === Quaternion{Node}
end

@testitem "FastDifferentiation: Jacobian of quaternion division" tags=[:unit] setup=[SymbolicReference] begin
    import FastDifferentiation
    u = FastDifferentiation.make_variables(:u, 8)
    D = quaternion(u[1:4]...) / quaternion(u[5:8]...)
    J = FastDifferentiation.jacobian(collect(components(D)), u)
    point = [0.3, -0.4, 0.5, 0.7, 0.2, 0.9, -0.6, 0.1]
    computed = FastDifferentiation.make_function(J, u)(point)
    # Reference: q / p = q p̄ / |p|², from explicit components
    function divide(x)
        q1, q2, q3, q4, p1, p2, p3, p4 = x
        n = p1^2 + p2^2 + p3^2 + p4^2
        [
            (q1*p1 + q2*p2 + q3*p3 + q4*p4) / n,
            (-q1*p2 + q2*p1 - q3*p4 + q4*p3) / n,
            (-q1*p3 + q2*p4 + q3*p1 - q4*p2) / n,
            (-q1*p4 - q2*p3 + q3*p2 + q4*p1) / n,
        ]
    end
    expected = SymbolicReference.bigfloat_jacobian(divide, point)
    @test maximum(abs, computed .- expected) < 100eps(Float64)
end

@testitem "FastDifferentiation: Jacobian through rotor normalization" tags=[:unit] setup=[SymbolicReference] begin
    import FastDifferentiation
    u = FastDifferentiation.make_variables(:u, 4)
    point = [0.3, -0.4, 0.5, 0.7]

    # Normalization alone is differentiated correctly.
    J = FastDifferentiation.jacobian(collect(components(rotor(u...))), u)
    computed = FastDifferentiation.make_function(J, u)(point)
    expected = SymbolicReference.bigfloat_jacobian(x -> collect(components(rotor(x...))), point)
    @test maximum(abs, computed .- expected) < 100eps(Float64)

    # FastDifferentiation 0.4.5 returns wrong derivatives for expressions in which several
    # products share a common subexpression, which rotor normalization produces.  This is an
    # upstream bug; the function values are correct.  A reproduction is kept with the notes
    # of the 4.4.5 patch, as `scratchpad/upstream/fastdifferentiation_1.jl`.
    M = to_rotation_matrix(rotor(u...))
    J = FastDifferentiation.jacobian(vec(M), u)
    computed = FastDifferentiation.make_function(J, u)(point)
    expected = SymbolicReference.bigfloat_jacobian(
        x -> vec(to_rotation_matrix(rotor(x...))), point
    )
    @test_broken maximum(abs, computed .- expected) < 100eps(Float64)
end

@testitem "Latexify: latexify agrees with the text/latex form" tags=[:unit] begin
    import Latexify, Symbolics
    Symbolics.@variables a b c d
    for q ∈ (
        quaternion(1.0, 2, 3, 4), Quaternion{Float64}(1, 2, 3, 4e-9), QuatVec{Float64}(1, 2, 3, 4),
        rotor(1, 3, 3, 9), quaternion(a, b, c, d), quaternion(a - b, b*c, c/d, d + a)
    )
        @test string(Latexify.latexify(q)) == repr(MIME("text/latex"), q)
    end
    @test string(Latexify.latexify(quaternion(1.0, 2, 3, 4))) ==
        "\$1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}\$"

    # Arrays of quaternions use the same form for each element, without Unicode units.
    s = string(Latexify.latexify([quaternion(1.0, 2, 3, 4), quaternion(5.0, 6, 7, 8)]))
    @test occursin("1.0 + 2.0\\,\\mathbf{i} + 3.0\\,\\mathbf{j} + 4.0\\,\\mathbf{k}", s)
    @test occursin("5.0 + 6.0\\,\\mathbf{i} + 7.0\\,\\mathbf{j} + 8.0\\,\\mathbf{k}", s)
    @test !occursin("𝐢", s)

    # NaN is set upright, as text, rather than as a product of italic letters.
    @test repr(MIME("text/latex"), quaternion(NaN, Inf, -Inf, 0)) ==
        "\$\\mathrm{NaN} + \\infty\\,\\mathbf{i} - \\infty\\,\\mathbf{j} + 0.0\\,\\mathbf{k}\$"
end
