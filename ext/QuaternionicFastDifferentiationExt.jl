module QuaternionicFastDifferentiationExt

using StaticArrays: SVector
using Latexify: latexify
import Quaternionic: AbstractQuaternion, Quaternion, Rotor, QuatVec,
    quaternion, rotor, quatvec,
    QuatVecF64, RotorF64, QuaternionF64,
    wrapper, components, _pm_ascii, _pm_latex
using PrecompileTools
isdefined(Base, :get_extension) ? (using FastDifferentiation) : (using ..FastDifferentiation)


### TYPE PIRACY!!!
# The following Base functions should be added to FastDifferentiation itself; meanwhile
# we'll define them here.

# This is a workaround that should be fixed in FastDifferentiation.jl; Elaine will probably
# make a PR to fix this.  The problem is that FD only defines
# `Base.promote_rule(::Type{<:Real}, ::Type{Node})`.  But julia/src/bool.jl defines
# `promote_rule(::Type{Bool}, ::Type{T}) where T<:Number` and `Node <: Number`, so there is
# an ambiguity.
Base.promote_rule(::Type{Bool}, ::Type{FastDifferentiation.Node}) = FastDifferentiation.Node

# These are essentially copied from Symbolics.jl:
# https://github.com/JuliaSymbolics/Symbolics.jl/blob/e4c328103ece494eaaab2a265524a64bfbe43dbd/src/num.jl#L31-L34
Base.eps(::Type{FastDifferentiation.Node}) = FastDifferentiation.Node(0)
Base.typemin(::Type{FastDifferentiation.Node}) = FastDifferentiation.Node(-Inf)
Base.typemax(::Type{FastDifferentiation.Node}) = FastDifferentiation.Node(Inf)
Base.float(x::FastDifferentiation.Node) = x

# This one is needed because julia/base/float.jl only defines `isinf` for `Real`, but `Node
# <: Number`.  (See https://github.com/brianguenter/FastDifferentiation.jl/issues/73)
Base.isinf(x::FastDifferentiation.Node) = !isnan(x) & !isfinite(x)


normalize(v::AbstractVector{FastDifferentiation.Node}) = v ./ √sum(x->x^2, v)


### Functions that used to appear in quaternion.jl
quaternion(w::FastDifferentiation.Node) = quaternion(SVector{4}(w, false, false, false))
rotor(w::FastDifferentiation.Node) = rotor(SVector{4}(one(w), false, false, false))
quatvec(w::FastDifferentiation.Node) = quatvec(SVector{4,typeof(w)}(false, false, false, false))
for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:FastDifferentiation.Node}) where {OP} = quaternion
        wrapper(::Type{<:FastDifferentiation.Node}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = quaternion
    end
end
let NT = FastDifferentiation.Node
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
let T = FastDifferentiation.Node
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = quaternion
        end
    end
end
Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:FastDifferentiation.Node} =
    wrapper(Q){promote_type(eltype(Q), S)}


# function _pm_ascii(x::FastDifferentiation.Node)
#     # Utility function to print a component of a quaternion
#     s = "$x"
#     if s[1] ∉ "+-"
#         s = "+" * s
#     end
#     if occursin(r"[+^/-]", s[2:end])
#         if s[1] == '+'
#             s = " + " * "(" * s[2:end] * ")"
#         else
#             s = " + " * "(" * s * ")"
#         end
#     else
#         s = " " * s[1] * " " * s[2:end]
#     end
#     s
# end
# function _pm_latex(x::Num)
#     # Utility function to print a component of a quaternion in LaTeX
#     s = latexify(x, env=:raw, bracket=true)
#     if s[1] ∉ "+-"
#         s = "+" * s
#     end
#     if occursin(r"[+^/-]", s[2:end])
#         if s[1] == '+'
#             s = " + " * "\\left(" * s[2:end] * "\\right)"
#         else
#             s = " + " * "\\left(" * s * "\\right)"
#         end
#     else
#         s = " " * s[1] * " " * s[2:end]
#     end
#     s
# end


# # Broadcast-like operations from FastDifferentiation
# # (d::FastDifferentiation.Operator)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# # (d::FastDifferentiation.Operator)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::Quaternion) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::Rotor) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))


### Functions that used to appear in algebra.jl
for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    let TB = FastDifferentiation.Node
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end
let S = FastDifferentiation.Node
    @eval begin
        Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(*), $S)(s*components(p))
        Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion} = wrapper($S, Val(*), Q)(s*components(p))
        Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(/), $S)(components(p)/s)
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
    FastDifferentiation.@variables w x y z a b c d e
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
        𝓇(𝓋)
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
