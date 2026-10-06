module QuaternionicSymbolicsExt

using StaticArrays: SVector
import Quaternionic: absvec,
    AbstractQuaternion, Quaternion, Rotor, QuatVec,
    quaternion, rotor, quatvec,
    QuatVecF64, RotorF64, QuaternionF64,
    wrapper, components, basetype, _pm_ascii
using PrecompileTools
isdefined(Base, :get_extension) ? (using Symbolics) : (using ..Symbolics)


Base.abs(q::AbstractQuaternion{Symbolics.Num}) = √sum(x->x^2, components(q))
Base.abs(q::QuatVec{Symbolics.Num}) = √sum(x->x^2, vec(q))
absvec(q::AbstractQuaternion{Symbolics.Num}) = √sum(x->x^2, vec(q))

# The source of `log` and of non-integer powers of a `Rotor` chooses among branches by
# comparing the values of the components, which symbolic components do not have.  Instead,
# these methods return the general formula, log(q) = log|q| + atan(|v⃗|, w) v⃗/|v⃗|, which is
# valid everywhere except on the real axis, where it has a removable singularity (or, on
# the negative real axis, a branch cut).  As with `exp`, it should not be evaluated there
# by direct substitution.  The formula for a `Rotor` assumes unit norm.
function Base.log(q::Quaternion{Symbolics.Num})
    a = absvec(q)
    log(abs(q)) + atan(a, q[1]) * (quatvec(q) / a)
end
function Base.log(q::Rotor{Symbolics.Num})
    a = absvec(q)
    atan(a, q[1]) * (quatvec(q) / a)
end
# Integer powers of a `Rotor` have no branches, so only the other exponents are listed.
Base.:^(q::Rotor{Symbolics.Num}, s::AbstractFloat) = exp(s * log(q))
Base.:^(q::Rotor{Symbolics.Num}, s::Rational) = exp(s * log(q))
Base.:^(q::Rotor{Symbolics.Num}, s::Symbolics.Num) = exp(s * log(q))


### Functions that used to appear in quaternion.jl
# The other components are filled with `zero(w)` rather than `false`, because `Num(false)`
# prints as `false`.
quaternion(w::Symbolics.Num) = quaternion(SVector{4}(w, zero(w), zero(w), zero(w)))
function rotor(w::Symbolics.Num)
    # When `w` wraps a numerical constant, the components are those of src's `rotor(v)` for
    # that constant, so that the result is `±1` for a finite nonzero constant and a NaN rotor
    # for a zero, NaN, or infinite constant.  The sign of any other expression is unknown, so
    # the identity rotor is returned without symbolic normalization.
    v = Symbolics.value(w)
    c = if v isa Real
        SVector{4,Symbolics.Num}(components(rotor(v)))
    else
        SVector{4,Symbolics.Num}(one(w), zero(w), zero(w), zero(w))
    end
    Rotor{Symbolics.Num}(c)
end
quatvec(w::Symbolics.Num) = quatvec(SVector{4}(zero(w), zero(w), zero(w), zero(w)))
for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Symbolics.Num}) where {OP} = quaternion
        wrapper(::Type{<:Symbolics.Num}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = quaternion
    end
end
let NT = Symbolics.Num
    for QT ∈ (QuatVec,)
        for OP ∈ (Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = quatvec
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = quatvec
            end
        end
    end
    for QT ∈ (Rotor,)
        for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = quaternion
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = quaternion
            end
        end
    end
end
let T = Symbolics.Num
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval begin
            wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = quaternion
            wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = quaternion
        end
    end
end
# This method resolves an ambiguity with Symbolics' rule
# `promote_rule(::Type{<:Number}, ::Type{Num})`.
Base.promote_rule(::Type{Q}, ::Type{Symbolics.Num}) where {Q<:AbstractQuaternion} =
    wrapper(Q){promote_type(basetype(Q), Symbolics.Num)}
# As in src/quaternion.jl, a `QuatVec` or a `Rotor` promoted with a scalar becomes a
# `Quaternion`.  These methods also resolve ambiguities with the rules in src.
for QT ∈ (QuatVec, Rotor)
    @eval Base.promote_rule(::Type{$QT{T}}, ::Type{Symbolics.Num}) where {T<:Number} =
        Quaternion{promote_type(T, Symbolics.Num)}
end


### Functions that used to appear in base.jl

# A symbolic expression is treated as zero when it simplifies to zero.  Other numbers are
# tested directly.
simplifies_to_zero(x::Symbolics.Num) = iszero(Symbolics.simplify(x; expand=true))
simplifies_to_zero(z::Complex) = simplifies_to_zero(real(z)) && simplifies_to_zero(imag(z))
simplifies_to_zero(x::Number) = iszero(x)

# A `QuatVec` always stores a zero scalar part, so comparing all four components gives the
# right answer for every combination of quaternion types, and a `QuatVec` equals a scalar
# exactly when both are zero.
function symbolic_equal(q1::AbstractQuaternion, q2::AbstractQuaternion)
    (
        simplifies_to_zero(q1[1]-q2[1]) &&
        simplifies_to_zero(q1[2]-q2[2]) &&
        simplifies_to_zero(q1[3]-q2[3]) &&
        simplifies_to_zero(q1[4]-q2[4])
    )
end
function symbolic_equal(q::AbstractQuaternion, x::Number)
    (
        simplifies_to_zero(q[1]-x) &&
        simplifies_to_zero(q[2]) &&
        simplifies_to_zero(q[3]) &&
        simplifies_to_zero(q[4])
    )
end

# Each of these combinations of argument types is needed to avoid method ambiguities with
# the methods in src/base.jl and in Symbolics.
let NT = Symbolics.Num
    QTs = (AbstractQuaternion, QuatVec, AbstractQuaternion{NT}, QuatVec{NT})
    for QT1 ∈ QTs, QT2 ∈ QTs
        if QT1 <: AbstractQuaternion{NT} || QT2 <: AbstractQuaternion{NT}
            @eval Base.:(==)(q1::$QT1, q2::$QT2) = symbolic_equal(q1, q2)
        end
    end
    for QT ∈ QTs
        @eval begin
            Base.:(==)(q::$QT, x::$NT) = symbolic_equal(q, x)
            Base.:(==)(x::$NT, q::$QT) = symbolic_equal(q, x)
        end
    end
    for QT ∈ (AbstractQuaternion{NT}, QuatVec{NT})
        @eval begin
            Base.:(==)(q::$QT, x::Number) = symbolic_equal(q, x)
            Base.:(==)(x::Number, q::$QT) = symbolic_equal(q, x)
        end
    end
end

function _pm_ascii(x::Symbolics.Num)
    # Utility function to print a component of a quaternion
    s = "$x"
    if s[1] ∉ "+-"
        s = "+" * s
    end
    if occursin(r"[+^/-]", s[2:end])
        if s[1] == '+'
            s = " + " * "(" * s[2:end] * ")"
        else
            s = " + " * "(" * s * ")"
        end
    else
        s = " " * s[1] * " " * s[2:end]
    end
    s
end


# Broadcast-like operations from Symbolics
(d::Symbolics.Differential)(q::Quaternion) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
(d::Symbolics.Differential)(q::Rotor) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
(d::Symbolics.Differential)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))


### Functions that used to appear in algebra.jl
for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    let TB = Symbolics.Num
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end
let S = Symbolics.Num
    @eval begin
        # These need to be explicit to avoid method ambiguities between SVector and Num
        function Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion}
            wrapper(Q, Val(*), $S)(s*p[1], s*p[2], s*p[3], s*p[4])
        end
        function Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion}
            wrapper($S, Val(*), Q)(s*p[1], s*p[2], s*p[3], s*p[4])
        end
        function Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion}
            wrapper(Q, Val(/), $S)(p[1] / s, p[2] / s, p[3] / s, p[4] / s)
        end
        function Base.:/(s::$S, p::Q) where {Q<:AbstractQuaternion}
            f = s / abs2(p)
            wrapper($S, Val(/), Q)(p[1] * f, -p[2] * f, -p[3] * f, -p[4] * f)
        end
    end
end


# Pre-compilation

@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the
    # size of the precompile file and potentially make loading faster.
    Symbolics.@variables w x y z a b c d e
    s = randn(Float64)
    v = randn(QuatVecF64)
    r = randn(RotorF64)
    q = randn(QuaternionF64)
    𝓈 = w
    𝓋 = quatvec(x, y, z)
    𝓇 = rotor(a, b, c, d)
    𝓆 = quaternion(w, x, y, z)

    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether they belong to
        # this package or not (on Julia 1.8 and higher)
        r(v)
        Symbolics.simplify.(𝓇(𝓋))
        # Tuples, unlike vectors, keep the type of each element, so that mixed operations
        # are compiled too.
        for a ∈ (s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆)
            conj(a)
            for b ∈ (s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆)
                a * b
                a / b
                a + b
                a - b
            end
        end

    end
end


end # module
