# Definitions shared by the test items that sweep over many element types
#
# The test items load this module with `setup=[TestUtils]`, which makes its exported names
# available in them.  The module also adds methods to `eps` and `≈` for exact and symbolic
# element types, so that a tolerance such as `rtol=eps(T)` can be written once for every
# element type.  These are methods of `Base` functions, so they remain in effect for the
# rest of the test process once any test item has loaded this module.

@testmodule TestUtils begin
    using Quaternionic
    import Symbolics

    export FloatTypes, IntTypes, SymbolicTypes, Types, PrimitiveTypes, QTypes
    export isapproxexpanded

    # `FloatTypes` and `IntTypes` must be in descending order of width.
    const FloatTypes = [BigFloat, Float64, Float32, Float16]
    const IntTypes = [BigInt, Int128, Int64, Int32, Int16, Int8]
    const SymbolicTypes = [Symbolics.Num]
    const Types = [FloatTypes...; IntTypes...; SymbolicTypes...]
    const PrimitiveTypes = [T for T in Types if isbitstype(T)]
    const QTypes = [Quaternion, Rotor, QuatVec]

    # Exact and symbolic arithmetic make no rounding errors, so their tolerance is zero.
    Base.eps(::Quaternion{T}) where {T} = eps(T)
    Base.eps(T::Type{<:Integer}) = zero(T)
    Base.eps(n::Symbolics.Num) = zero(n)

    # Whether the symbolic expression `diff` simplifies to zero
    function symbolic_iszero(diff)
        d = Symbolics.simplify(diff; expand=true)
        iszero(d) || iszero(Symbolics.simplify(d^2; expand=true))
    end

    Base.:≈(a::Symbolics.Num, b::Number; kwargs...) = symbolic_iszero(a - b)
    Base.:≈(a::Number, b::Symbolics.Num; kwargs...) = symbolic_iszero(a - b)
    Base.:≈(a::AbstractQuaternion{Symbolics.Num}, b::AbstractQuaternion{Symbolics.Num}; kwargs...) =
        all(iszero(Symbolics.simplify(x - y; expand=true)) for (x, y) in zip(components(a), components(b)))

    # Symbolics defines `isapprox(::Num, ::Num)` itself, and its method returns `false`
    # unless the difference of the two expressions is already a number.  Redefining that
    # method here would overwrite Symbolics' own, so tests that compare two symbolic
    # expressions that are equal only after expansion use this function instead, which
    # falls back to `isapprox` for every other type.
    isapproxexpanded(a, b; kwargs...) = isapprox(a, b; kwargs...)
    isapproxexpanded(a::Symbolics.Num, b::Symbolics.Num; kwargs...) = symbolic_iszero(a - b)
end
