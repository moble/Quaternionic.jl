module QuaternionicChainRulesCoreExt

using Quaternionic
using StaticArrays
isdefined(Base, :get_extension) ? (using ChainRulesCore) : (using ..ChainRulesCore)


function ChainRulesCore.rrule(::Type{QT}, arg) where {QT<:AbstractQuaternion}
    AbstractQuaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
    return QT(arg), AbstractQuaternion_pullback
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
