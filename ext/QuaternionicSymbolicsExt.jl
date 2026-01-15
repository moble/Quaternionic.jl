module QuaternionicSymbolicsExt

using StaticArrays: SVector
import Quaternionic: normalize, absvec,
    AbstractQuaternion, Quaternion, Rotor, QuatVec,
    quaternion, rotor, quatvec,
    QuatVecF64, RotorF64, QuaternionF64,
    wrapper, components, basetype, _pm_ascii
using PrecompileTools
isdefined(Base, :get_extension) ? (using Symbolics) : (using ..Symbolics)


normalize(v::AbstractVector{Symbolics.Num}) = v ./ √sum(x->x^2, v)
Base.abs(q::AbstractQuaternion{Symbolics.Num}) = √sum(x->x^2, components(q))
Base.abs(q::QuatVec{Symbolics.Num}) = √sum(x->x^2, vec(q))
absvec(q::AbstractQuaternion{Symbolics.Num}) = √sum(x->x^2, vec(q))


### Functions that used to appear in quaternion.jl
quaternion(w::Symbolics.Num) = quaternion(SVector{4}(w, false, false, false))
rotor(w::Symbolics.Num) = rotor(SVector{4}(one(w), false, false, false))
quatvec(w::Symbolics.Num) = quatvec(SVector{4,typeof(w)}(false, false, false, false))
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
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = quaternion
        end
    end
end
Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Symbolics.Num} =
    wrapper(Q){promote_type(basetype(Q), S)}
Base.promote_rule(::Type{Q}, ::Type{Symbolics.Num}) where {Q<:AbstractQuaternion} =
    wrapper(Q){promote_type(basetype(Q), Symbolics.Num)}


### Functions that used to appear in base.jl
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1[1]-q2[1]; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::AbstractQuaternion)
    (
        iszero(Symbolics.simplify(q1[1]-q2[1]; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Number)
    (
        iszero(Symbolics.simplify(q1[1]-q2; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]; expand=true))
    )
end
function Base.:(==)(q1::Number, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1-q2[1]; expand=true)) &&
        iszero(Symbolics.simplify(q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q2[4]; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::Symbolics.Num)
    (
        iszero(Symbolics.simplify(q1[1]-q2; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]; expand=true))
    )
end
function Base.:(==)(q1::Symbolics.Num, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1-q2[1]; expand=true)) &&
        iszero(Symbolics.simplify(q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q2[4]; expand=true))
    )
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::AbstractQuaternion{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q2[1]; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::QuatVec{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1[1]; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::QuatVec{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::AbstractQuaternion)
    (
        iszero(q2[1]) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q2::AbstractQuaternion, q1::QuatVec{Symbolics.Num})
    (
        iszero(q2[1]) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::AbstractQuaternion{Symbolics.Num}, q2::QuatVec)
    (
        iszero(Symbolics.simplify(q1[1]; expand=true)) &&
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::QuatVec)
    (
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q2::QuatVec, q1::QuatVec{Symbolics.Num})
    (
        iszero(Symbolics.simplify(q1[2]-q2[2]; expand=true)) &&
        iszero(Symbolics.simplify(q1[3]-q2[3]; expand=true)) &&
        iszero(Symbolics.simplify(q1[4]-q2[4]; expand=true))
    )
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::Number)
    false
end
function Base.:(==)(q1::Number, q2::QuatVec{Symbolics.Num})
    false
end
function Base.:(==)(q1::QuatVec{Symbolics.Num}, q2::Symbolics.Num)
    false
end
function Base.:(==)(q1::Symbolics.Num, q2::QuatVec{Symbolics.Num})
    false
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
# (d::Symbolics.Operator)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::Symbolics.Operator)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))
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
        for a ∈ [s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆]
            conj(a)
            for b ∈ [s, v, r, q, 𝓈, 𝓋, 𝓇, 𝓆]
                a * b
                a / b
                a + b
                a - b
            end
        end

    end
end


end # module
