module QuaternionicFastDifferentiationExt

using StaticArrays: SVector
import Quaternionic: Quaternionic, normalize, absvec,
    AbstractQuaternion, Quaternion, Rotor, QuatVec,
    quaternion, rotor, quatvec,
    QuatVecF64, RotorF64, QuaternionF64,
    wrapper, components, basetype
using PrecompileTools
isdefined(Base, :get_extension) ?
    (using FastDifferentiation: FastDifferentiation, Node) :
    (using ..FastDifferentiation: FastDifferentiation, Node)


normalize(v::AbstractVector{Node}) = v ./ √sum(x->x^2, v)
Base.abs(q::AbstractQuaternion{Node}) = √sum(x->x^2, components(q))
Base.abs(q::QuatVec{Node}) = √sum(x->x^2, vec(q))
absvec(q::AbstractQuaternion{Node}) = √sum(x->x^2, vec(q))


### Functions that used to appear in quaternion.jl
quaternion(w::Node) = quaternion(SVector{4}(w, false, false, false))
rotor(w::Node) = rotor(SVector{4}(one(w), false, false, false))
quatvec(w::Node) = quatvec(SVector{4,typeof(w)}(false, false, false, false))
for QT1 ∈ (AbstractQuaternion, Quaternion, QuatVec, Rotor)
    @eval begin
        wrapper(::Type{<:$QT1}, ::Val{OP}, ::Type{<:Node}) where {OP} = quaternion
        wrapper(::Type{<:Node}, ::Val{OP}, ::Type{<:$QT1}) where {OP} = quaternion
    end
end
let NT = Node
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
let T = Node
    for OP ∈ (Val{+}, Val{-}, Val{*}, Val{/})
        @eval wrapper(::Type{<:Quaternion}, ::$OP, ::Type{<:$T}) = quaternion
        if T !== Quaternion
            @eval wrapper(::Type{<:$T}, ::$OP, ::Type{<:Quaternion}) = quaternion
        end
    end
end
Base.promote_rule(::Type{Q}, ::Type{S}) where {Q<:AbstractQuaternion,S<:Node} =
    wrapper(Q){promote_type(basetype(Q), S)}


# # Broadcast-like operations from FastDifferentiation
# # (d::FastDifferentiation.Operator)(q::QT) where {QT<:AbstractQuaternion} = QT(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# # (d::FastDifferentiation.Operator)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::Quaternion) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::Rotor) = quaternion(d(q[1]), d(q[2]), d(q[3]), d(q[4]))
# (d::FastDifferentiation.Differential)(q::QuatVec) = quatvec(d(q[2]), d(q[3]), d(q[4]))


### Functions that used to appear in algebra.jl
for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    let TB = Node
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end
let S = Node
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

# Here, we disable FastDifferentiation support for functions that rely on conditionals;
# otherwise, they will fail with errors like
#
#    TypeError: non-boolean (Node) used in boolean context
#
# This could be removed once FastDifferentiation supports conditionals, and these functions
# are tested.
for method ∈ (
    :(Quaternionic.to_euler_phases(::AbstractQuaternion{Node})),
    :(Quaternionic.from_euler_phases(::Complex{Node}, ::Complex{Node}, ::Complex{Node})),
    :(Base.log(::Quaternion{Node})),
    :(Base.log(::Rotor{Node})),
    :(Base.exp(::Quaternion{Node})),
    :(Base.exp(::QuatVec{Node})),
    :(Base.sqrt(::Quaternion{Node})),
    :(Base.sqrt(::QuatVec{Node})),
    :(Base.sqrt(::Rotor{Node})),
    :(Base.:^(::Rotor{Node}, ::Number)),
    :(Base.:^(::Quaternion{Node}, ::Integer)),
    :(Base.:^(::QuatVec{Node}, ::Integer)),
    :(Base.:^(::Rotor{Node}, ::Integer)),
)
    func = first(split(string(method), '('))
    @eval begin
        $(method) = error(
            """FastDifferentiation does not yet support conditionals involving FD variables,
            which are needed to implement its derivative of `$(func)` involving quaternions.
            """
        )
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
                # a / b  # TODO: Uncomment when https://github.com/brianguenter/FastDifferentiation.jl/issues/98 is fixed
                a + b
                a - b
            end
        end

    end
end


end # module
