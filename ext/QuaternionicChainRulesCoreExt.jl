module QuaternionicChainRulesCoreExt

using Quaternionic
using StaticArrays
isdefined(Base, :get_extension) ? (using ChainRulesCore) : (using ..ChainRulesCore)


function ChainRulesCore.rrule(::Type{QT}, arg::AbstractVector) where {QT<:AbstractQuaternion}
    @info 1
    AbstractQuaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
    @info QT(arg)
    println()
    return QT(arg), AbstractQuaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w::AbstractQuaternion) where {QT<:AbstractQuaternion}
    @info 2
    Quaternion_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
    return QT(w), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w, x, y, z) where {QT<:AbstractQuaternion}
    @info 3
    Quaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
    v = @SVector[w, x, y, z]
    return QT(v), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, x, y, z) where {QT<:AbstractQuaternion}
    @info 4
    Quaternion_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
    v = @SVector[false, x, y, z]
    return QT(v), Quaternion_pullback
end
function ChainRulesCore.rrule(::Type{QT}, w::Number) where {QT<:AbstractQuaternion}
    @info 5
    Quaternion_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
    v = @SVector[w, false, false, false]
    return QT(v), Quaternion_pullback
end

function ChainRulesCore.rrule(::Type{Rotor}, w, x, y, z)
    @info "rrule(::Type{Rotor}, w, x, y, z)"
    ChainRulesCore.rrule(rotor, w, x, y, z)
end
function ChainRulesCore.rrule(::Type{Rotor}, x, y, z)
    @info "rrule(::Type{Rotor}, x, y, z)"
    ChainRulesCore.rrule(rotor, x, y, z)
end
function ChainRulesCore.rrule(::Type{Rotor}, w::AbstractVector)
    @info "rrule(::Type{Rotor}, w::AbstractVector)"
    ChainRulesCore.rrule(rotor, w)
end
ChainRulesCore.rrule(::Type{Rotor}, w::AbstractQuaternion) = ChainRulesCore.rrule(rotor, w)
# function ChainRulesCore.rrule(T::Type{Rotor}, w::AbstractQuaternion)
#     Rotor_pullback(Δquat) = (NoTangent(), false * unthunk(Δquat))
#     return rotor(w), Rotor_pullback
# end
function ChainRulesCore.rrule(::Type{Rotor}, w::Number)
    @info "rrule(::Type{Rotor}, w::Number)"
    ChainRulesCore.rrule(rotor, w)
end

# function ChainRulesCore.rrule(::Type{Rotor{T}}, arg::AbstractVector) where {T}
#     @info 6
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
#     return Rotor{T}(arg), Rotor_pullback
# end
function ChainRulesCore.rrule(::typeof(rotor), arg::AbstractVector)
    @info 7 rotor(arg) arg
    rotor_pullback(Δquat) = (NoTangent(), false * components(unthunk(Δquat)))
    return rotor(arg), rotor_pullback
end
# function ChainRulesCore.rrule(::Type{QT}, arg::AbstractVector) where {QT<:Rotor}
#     @info 8
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
#     return QT(arg), Rotor_pullback
# end

# function ChainRulesCore.rrule(::Type{Rotor{T}}, w::AbstractQuaternion) where {T}
#     @info 9
#     Rotor_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
#     return Rotor{T}(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::typeof(rotor), w::AbstractQuaternion)
#     @info 10
#     Rotor_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
#     return rotor(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{QT}, w::AbstractQuaternion) where {QT<:Rotor}
#     @info 11
#     Rotor_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
#     return QT(w), Rotor_pullback
# end

# function ChainRulesCore.rrule(::Type{Rotor{T}}, w, x, y, z) where {T}
#     @info 12
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
#     return Rotor{T}(w, x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::typeof(rotor), w, x, y, z)
#     @info 13
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
#     @info "rrule(::typeof(rotor), w, x, y, z)" rotor(w, x, y, z)
#     return rotor(w, x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{QT}, w, x, y, z) where {QT<:Rotor}
#     @info 14
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
#     return QT(w, x, y, z), Rotor_pullback
# end

# function ChainRulesCore.rrule(::Type{Rotor{T}}, x, y, z) where {T}
#     @info 15
#     Rotor_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
#     return Rotor{T}(x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::typeof(rotor), x, y, z)
#     @info 16
#     Rotor_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
#     return rotor(x, y, z), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{QT}, x, y, z) where {QT<:Rotor}
#     @info 17
#     Rotor_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
#     return QT(x, y, z), Rotor_pullback
# end

# function ChainRulesCore.rrule(::Type{Rotor{T}}, w::Number) where {T}
#     @info 18
#     Rotor_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
#     return Rotor{T}(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::typeof(rotor), w::Number)
#     @info 19
#     Rotor_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
#     return rotor(w), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{QT}, w::Number) where {QT<:Rotor}
#     @info 20
#     Rotor_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
#     return QT(w), Rotor_pullback
# end


# function ChainRulesCore.rrule(::Type{Rotor{QT}}, arg::AbstractVector{VT}) where {QT, VT}
#     @info "X"
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
#     return Rotor{QT}(arg), Rotor_pullback
# end
# function ChainRulesCore.rrule(::Type{Rotor{T}}, w, x, y, z) where {T}
#     @info "Y"
#     Rotor_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
#     v = @SVector[w, x, y, z]
#     return Rotor{eltype(v)}(v), Rotor_pullback
# end

function ChainRulesCore.rrule(::Type{QuatVec{QT}}, arg::AbstractVector{VT}) where {QT, VT}
    @info "rrule(::Type{QuatVec{QT}}, arg::AbstractVector{VT})"
    function QuatVec_pullback(Δquat)
        c = components(unthunk(Δquat))
        (NoTangent(), @SVector[zero(eltype(c)), c[2], c[3], c[4]])
    end
    v = @SVector[false, arg[begin+1], arg[begin+2], arg[begin+3]]
    return QuatVec{QT}(v), QuatVec_pullback
end
function ChainRulesCore.rrule(::Type{QuatVec{T}}, w, x, y, z) where {T}
    @info "rrule(::Type{QuatVec{T}}, w, x, y, z)"
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
