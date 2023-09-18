module QuaternionicSymbolicsExt

using StaticArrays: SVector
import Quaternionic: AbstractQuaternion, Quaternion, Rotor, QuatVec,
    QuatVecF64, RotorF64, QuaternionF64, wrapper, components
using PrecompileTools
isdefined(Base, :get_extension) ? (using Symbolics) : (using ..Symbolics)


### Functions that used to appear in quaternion.jl
quaternion(w::Symbolics.Num) = quaternion(SVector{4}(w, false, false, false))
rotor(w::Symbolics.Num) = rotor(SVector{4}(one(w), false, false, false))
quatvec(w::Symbolics.Num) = quatvec(SVector{4,typeof(w)}(false, false, false, false))
for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Symbolics.Num}) where {OP} = Quaternion
        wrapper(::Type{<:Symbolics.Num}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = Quaternion
    end
end
let NT = Symbolics.Num
    for QT ∈ (AbstractQuaternion, QuatVec)
        for OP ∈ (Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = $QT
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = $QT
            end
        end
    end
    for QT ∈ (Rotor,)
        for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
            @eval begin
                wrapper(::Type{<:$QT}, ::$OP, ::Type{<:$NT}) = Quaternion
                wrapper(::Type{<:$NT}, ::$OP, ::Type{<:$QT}) = Quaternion
            end
        end
    end
end
let T = Symbolics.Num
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = Quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = Quaternion
        end
    end
end
Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Symbolics.Num} =
    wrapper(Q){promote_type(eltype(Q), S)}


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
# Broadcast-like operations from Symbolics
# (d::Symbolics.Operator)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::Symbolics.Operator)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))
(d::Symbolics.Differential)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
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
        # your package or not (on Julia 1.8 and higher)
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
