module QuaternionicChainRulesCoreExt

using Quaternionic
import Quaternionic: _sinc, _cosc
using StaticArrays
isdefined(Base, :get_extension) ?
    (using ChainRulesCore; import ChainRulesCore: rrule, rrule_via_ad, RuleConfig, ProjectTo) :
    (using ..ChainRulesCore; import ...ChainRulesCore: rrule, rrule_via_ad, RuleConfig, ProjectTo)


## StaticArrays
# It's likely that StaticArrays will have its own ChainRulesCore extension someday, so we
# need to check if there is already a ProjectTo defined for SArray.  If so, we'll use that.
# If not, we'll define one here.
if !any(method->occursin("SArray", repr(method.sig)), methods(ProjectTo))
    # These are ripped from https://github.com/JuliaArrays/StaticArrays.jl/pull/1068
    function (project::ProjectTo{<:Tangent{<:Tuple}})(dx::SArray)
        dy = reshape(dx, axes(project.elements))  # allows for dx::OffsetArray
        dz = ntuple(i -> project.elements[i](dy[i]), length(project.elements))
        return ChainRulesCore.project_type(project)(dz...)
    end
    function ProjectTo(x::SArray{S,T}) where {S, T}
        return ProjectTo{SArray}(; element=ChainRulesCore._eltype_projectto(T), axes=S)
    end
    function (project::ProjectTo{SArray})(dx::AbstractArray{S,M}) where {S,M}
        return SArray{project.axes}(dx)
    end
    function rrule(::Type{T}, x::Tuple) where {T<:SArray}
        project_x = ProjectTo(x)
        Array_pullback(ȳ) = (NoTangent(), project_x(ȳ))
        return T(x), Array_pullback
    end
end


function rrule(::Type{QT}, arg::AbstractVector) where {QT<:AbstractQuaternion}
    AbstractQuaternion_pullback(Δquat) = (@show 1; (NoTangent(), components(unthunk(Δquat))))
    return QT(arg), AbstractQuaternion_pullback
end
function rrule(::Type{QT}, w::AbstractQuaternion) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (@show 2; (NoTangent(), unthunk(Δquat)))
    return QT(w), Quaternion_pullback
end
function rrule(::Type{QT}, w, x, y, z) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (@show 3; (NoTangent(), components(unthunk(Δquat))...))
    return QT(SVector{4}(w, x, y, z)), Quaternion_pullback
end
function rrule(::Type{QT}, x, y, z) where {QT<:AbstractQuaternion}
   Quaternion_pullback(Δquat) = (@show 4; (NoTangent(), vec(unthunk(Δquat))...))
    return QT(SVector{4}(false, x, y, z)), Quaternion_pullback
end
function rrule(::Type{QT}, w::Number) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (@show 5; (NoTangent(), real(unthunk(Δquat))))
    return QT(SVector{4}(w, false, false, false)), Quaternion_pullback
end

# function rrule(::Type{QuatVec{QT}}, arg::AbstractVector{VT}) where {QT, VT}
#     function QuatVec_pullback(Δquat)
#         c = components(unthunk(Δquat))
#         @info "6" QT typeof(Δquat) Δquat c typeof(arg) length(arg) arg[1] arg[2] arg[3] arg[4]
#         (NoTangent(), SVector{4}(zero(eltype(c)), c[2], c[3], c[4]))
#     end
#     v = SVector{4}(false, arg[begin+1], arg[begin+2], arg[begin+3])
#     return QuatVec{QT}(v), QuatVec_pullback
# end

# function rrule(::Type{QuatVec{T}}, w, x, y, z) where {T}
#     function QuatVec_pullback(Δquat)
#         c = components(unthunk(Δquat))
#         @info "7" T typeof(Δquat) Δquat c typeof(w) w x y z
#         (NoTangent(), zero(eltype(c)), c[2], c[3], c[4])
#     end
#     v = SVector{4}(false, x, y, z)
#     return QuatVec{eltype(v)}(v), QuatVec_pullback
# end

rrule(config::RuleConfig{>:HasReverseMode}, ::Type{Rotor}, args...) = rrule_via_ad(config, rotor, args...)
rrule(config::RuleConfig{>:HasReverseMode}, ::Type{QuatVec}, args...) = rrule_via_ad(config, quatvec, args...)


## Modified from `Complex` entries in ChainRulesCore.jl/src/projection.jl
ProjectTo(::T) where {T<:AbstractQuaternion} = ProjectTo{T}()
ProjectTo(x::AbstractQuaternion{<:Integer}) = ProjectTo(float(x))
for T in (
    QuaternionF16, QuaternionF32, QuaternionF64,
    RotorF16, RotorF32, RotorF64,
    QuatVecF16, QuatVecF32, QuatVecF64
)
    @eval ProjectTo(::$T) = ProjectTo{$T}()
end
function (::ProjectTo{QT})(dx::AbstractQuaternion{<:AbstractFloat}) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    #@info "ProjectTo{QT}(dx::AbstractQuaternion{<:AbstractFloat})" QT dx typeof(dx) convert(QT, dx)
    return convert(QT, dx)
end
function (::ProjectTo{QT})(dx::AbstractFloat) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    #@info "ProjectTo{QT}(dx::AbstractFloat)"
    return convert(QT, dx)
end
function (::ProjectTo{QT})(dx::AbstractQuaternion{<:Integer}) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    #@info "ProjectTo{QT}(dx::AbstractQuaternion{<:Integer})"
    return convert(QT, dx)
end
function (::ProjectTo{QT})(dx::Integer) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    #@info "ProjectTo{QT}(dx::Integer)"
    return convert(QT, dx)
end
function (project::ProjectTo{QT})(dx::Real) where {QT<:AbstractQuaternion}
    return project(QT(dx))
end
function (project::ProjectTo{<:Number})(dx::Tangent{QT}) where {QT<:AbstractQuaternion}
    project(QT(dx[:components]))
end


## Copied from `Complex` entries in ChainRulesCore.jl/src/tangent_types/abstract_zero.jl
for pattern ∈ 1:15
    T1 = iszero(pattern & 1) ? Number : AbstractZero
    T2 = iszero(pattern & 2) ? Number : AbstractZero
    T3 = iszero(pattern & 4) ? Number : AbstractZero
    T4 = iszero(pattern & 8) ? Number : AbstractZero
    @eval (::Type{QT})(w::$T1, x::$T2, y::$T3, z::$T4) where {QT<:AbstractQuaternion} = QT(w, x, y, z)
end


## Copied from `Complex` entries in ChainRulesCore.jl/src/tangent_types/thunks.jl
function (::Type{QT})(a::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a))
end
function (::Type{QT})(a::AbstractThunk, b::AbstractThunk, c::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a, b, c))
end
function (::Type{QT})(a::AbstractThunk, b::AbstractThunk, c::AbstractThunk, d::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a, b, c, d))
end


# Following ChainRules <https://juliadiff.org/ChainRulesCore.jl/stable/maths/complex.html>,
# we define derivatives of a function of a quaternion in terms of its components:
#
#    f(w + 𝐢*x + 𝐣*y + 𝐤*z) = s + 𝐢*t + 𝐣*u + 𝐤*v
#
# The `frule(Δw+𝐢*Δx+𝐣*Δy+𝐤*Δz)` should return
#
#    (∂s/∂w Δw + ∂s/∂x Δx + ∂s/∂y Δy + ∂s/∂z Δz)
#    + 𝐢 * (∂t/∂w Δw + ∂t/∂x Δx + ∂t/∂y Δy + ∂t/∂z Δz)
#    + 𝐣 * (∂u/∂w Δw + ∂u/∂x Δx + ∂u/∂y Δy + ∂u/∂z Δz)
#    + 𝐤 * (∂v/∂w Δw + ∂v/∂x Δx + ∂v/∂y Δy + ∂v/∂z Δz)
#
# while the `rrule(Δs+𝐢*Δt+𝐣*Δu+𝐤*Δv)` should return
#
#    (∂s/∂w Δs + ∂t/∂w Δt + ∂u/∂w Δu + ∂v/∂w Δv)
#    + 𝐢 * (∂s/∂x Δs + ∂t/∂x Δt + ∂u/∂x Δu + ∂v/∂x Δv)
#    + 𝐣 * (∂s/∂y Δs + ∂t/∂y Δt + ∂u/∂y Δu + ∂v/∂y Δv)
#    + 𝐤 * (∂s/∂z Δs + ∂t/∂z Δt + ∂u/∂z Δu + ∂v/∂z Δv)

function rrule(::typeof(exp), v::QuatVec{T}) where T
    x, y, z = vec(v)
    a2 = abs2vec(v)
    a = sqrt(a2)
    sinc = _sinc(a)
    cosc = _cosc(a)

    s = cos(a)
    t = x * sinc
    u = y * sinc
    v = z * sinc
    R = s + 𝐢*t + 𝐣*u + 𝐤*v

    ∂sinc∂x = cosc * x / a
    ∂sinc∂y = cosc * y / a
    ∂sinc∂z = cosc * z / a
    ∂s∂x = -x * sinc
    ∂s∂y = -y * sinc
    ∂s∂z = -z * sinc
    ∂t∂x = sinc + x * ∂sinc∂x
    ∂t∂y = x * ∂sinc∂y
    ∂t∂z = x * ∂sinc∂z
    ∂u∂x = y * ∂sinc∂x
    ∂u∂y = sinc + y * ∂sinc∂y
    ∂u∂z = y * ∂sinc∂z
    ∂v∂x = z * ∂sinc∂x
    ∂v∂y = z * ∂sinc∂y
    ∂v∂z = sinc + z * ∂sinc∂z

    function exp_pullback(ΔR)
        Δs, Δt, Δu, Δv = components(unthunk(ΔR))
        return (
            NoTangent(),
            𝐢 * (∂s∂x * Δs + ∂t∂x * Δt + ∂u∂x * Δu + ∂v∂x * Δv)
            + 𝐣 * (∂s∂y * Δs + ∂t∂y * Δt + ∂u∂y * Δu + ∂v∂y * Δv)
            + 𝐤 * (∂s∂z * Δs + ∂t∂z * Δt + ∂u∂z * Δu + ∂v∂z * Δv)
        )
    end

    return R, exp_pullback
end


end
