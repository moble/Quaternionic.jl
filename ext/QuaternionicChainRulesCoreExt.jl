module QuaternionicChainRulesCoreExt

using Quaternionic
using StaticArrays
isdefined(Base, :get_extension) ? (using ChainRulesCore) : (using ..ChainRulesCore)


function ChainRulesCore.rrule(::Type{QT}, arg::AbstractVector) where {QT<:AbstractQuaternion}
    AbstractQuaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
    return QT(arg), AbstractQuaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w::AbstractQuaternion) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
    return QT(w), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w, x, y, z) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
    return QT(@SVector[w, x, y, z]), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, x, y, z) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
    return QT(@SVector[false, x, y, z]), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w::Number) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
    return QT(@SVector[w, false, false, false]), Quaternion_pullback
end

ChainRulesCore.rrule(config::ChainRulesCore.RuleConfig, ::Type{Rotor}, args...) = ChainRulesCore.rrule_via_ad(config, rotor, args...)
# function ChainRulesCore.rrule(::Type{Rotor}, w, x, y, z)
#     Rotor_pullback(Δquat) = (NoTangent(), (false * components(unthunk(Δquat)))...)
#     return rotor(w, x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{Rotor}, x, y, z)
#     Rotor_pullback(Δquat) = (NoTangent(), (false * vec(unthunk(Δquat)))...)
#     return rotor(x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{Rotor}, w::AbstractVector)
#     Rotor_pullback(Δquat) = (NoTangent(), false * components(unthunk(Δquat)))
#     return rotor(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{Rotor}, w::AbstractQuaternion)
#     Rotor_pullback(Δquat) = (NoTangent(), false * unthunk(Δquat))
#     return rotor(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{Rotor}, w::Number)
#     Rotor_pullback(Δquat) = (NoTangent(), false * real(unthunk(Δquat)))
#     return rotor(w), Rotor_pullback
# end

function ChainRulesCore.rrule(::typeof(rotor), arg::AbstractVector)
    rotor_pullback(Δquat) = (NoTangent(), false * components(unthunk(Δquat)))
    return rotor(arg), rotor_pullback
end

function ChainRulesCore.rrule(::Type{QuatVec{QT}}, arg::AbstractVector{VT}) where {QT, VT}
    function QuatVec_pullback(Δquat)
        c = components(unthunk(Δquat))
        (NoTangent(), @SVector[zero(eltype(c)), c[2], c[3], c[4]])
    end
    v = @SVector[false, arg[begin+1], arg[begin+2], arg[begin+3]]
    return QuatVec{QT}(v), QuatVec_pullback
end
function ChainRulesCore.rrule(::Type{QuatVec{T}}, w, x, y, z) where {T}
    function QuatVec_pullback(Δquat)
        c = components(unthunk(Δquat))
        (NoTangent(), zero(eltype(c)), c[2], c[3], c[4])
    end
    v = @SVector[false, x, y, z]
    return QuatVec{eltype(v)}(v), QuatVec_pullback
end


## Modified from `Complex` entries in ChainRulesCore.jl/src/projection.jl
ChainRulesCore.ProjectTo(::T) where {T<:AbstractQuaternion} = ProjectTo{T}()
ChainRulesCore.ProjectTo(x::AbstractQuaternion{<:Integer}) = ProjectTo(float(x))
for T in (QuaternionF16, QuaternionF32, QuaternionF64)
    @eval ChainRulesCore.ProjectTo(::$T) = ChainRulesCore.ProjectTo{$T}()
end
function (::ChainRulesCore.ProjectTo{QT})(dx::AbstractQuaternion{<:AbstractFloat}) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    return convert(QT, dx)
end
function (::ChainRulesCore.ProjectTo{QT})(dx::AbstractFloat) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    return convert(QT, dx)
end
function (::ChainRulesCore.ProjectTo{QT})(dx::AbstractQuaternion{<:Integer}) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    return convert(QT, dx)
end
function (::ChainRulesCore.ProjectTo{QT})(dx::Integer) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    return convert(QT, dx)
end
function (project::ChainRulesCore.ProjectTo{QT})(dx::Real) where {QT<:AbstractQuaternion}
    return project(QT(dx))
end
function (project::ChainRulesCore.ProjectTo{<:Number})(dx::Tangent{QT}) where {QT<:AbstractQuaternion}
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


## StaticArrays
# It's likely that StaticArrays will have its own ChainRulesCore extension someday, so we
# need to check if there is already a ProjectTo defined for SArray.  If so, we'll use that.
# If not, we'll define one here.
if !any(method->occursin("SArray", repr(method.sig)), methods(ChainRulesCore.ProjectTo))
    # These are ripped from https://github.com/JuliaArrays/StaticArrays.jl/pull/1068
    function (project::ChainRulesCore.ProjectTo{<:Tangent{<:Tuple}})(dx::SArray)
        dy = reshape(dx, axes(project.elements))  # allows for dx::OffsetArray
        dz = ntuple(i -> project.elements[i](dy[i]), length(project.elements))
        return ChainRulesCore.project_type(project)(dz...)
    end
    function ChainRulesCore.ProjectTo(x::SArray{S,T}) where {S, T}
        return ChainRulesCore.ProjectTo{SArray}(; element=ChainRulesCore._eltype_projectto(T), axes=S)
    end
    function (project::ChainRulesCore.ProjectTo{SArray})(dx::AbstractArray{S,M}) where {S,M}
        return SArray{project.axes}(dx)
    end
    function ChainRulesCore.rrule(::Type{T}, x::Tuple) where {T<:SArray}
        project_x = ProjectTo(x)
        Array_pullback(ȳ) = (NoTangent(), project_x(ȳ))
        return T(x), Array_pullback
    end
end


end
