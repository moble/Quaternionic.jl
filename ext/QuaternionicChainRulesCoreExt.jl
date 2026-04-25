module QuaternionicChainRulesCoreExt

using Quaternionic
import Quaternionic: _sincu, _cossu
using StaticArrays
isdefined(Base, :get_extension) ?
    (using ChainRulesCore; import ChainRulesCore: rrule, rrule_via_ad, RuleConfig, ProjectTo) :
    (using ..ChainRulesCore; import ...ChainRulesCore: rrule, rrule_via_ad, RuleConfig, ProjectTo)


function rrule(::typeof(components), q::AbstractQuaternion)
    c = components(q)
    Π = ProjectTo(q)
    function components_pullback(Δ)
        Δs, Δt, Δu, Δv = Tuple(unthunk(Δ))
        Δq = typeof(q)(Δs, Δt, Δu, Δv)
        return NoTangent(), Π(Δq)
    end
    return c, components_pullback
end

function rrule(::typeof(Base.getindex), q::AbstractQuaternion, i::Integer)
    s = q[i]
    function quaternion_getindex_pullback(Δs)
        Δw = Δx = Δy = Δz = zero(s)
        @inbounds if i == 1
            Δw = unthunk(Δs)
        elseif i == 2
            Δx = unthunk(Δs)
        elseif i == 3
            Δy = unthunk(Δs)
        elseif i == 4
            Δz = unthunk(Δs)
        else
            throw(BoundsError(q, i))
        end
        return NoTangent(), typeof(q)(Δw, Δx, Δy, Δz), NoTangent()
    end
    return s, quaternion_getindex_pullback
end


##################
## Constructors ##
##################

# Quaternion
function rrule(::Type{QT}, arg::AbstractVector) where {QT<:AbstractQuaternion}
    # AbstractQuaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat)))
    # return QT(arg), AbstractQuaternion_pullback
        # Return a same-shaped tangent for the vector argument
    function AbstractQuaternion_pullback(Δquat)
        Δw, Δx, Δy, Δz = components(unthunk(Δquat))
        Δarg = similar(arg)
        if length(Δarg) == 4
            Δarg[1] = Δw; Δarg[2] = Δx; Δarg[3] = Δy; Δarg[4] = Δz
        else
            # Fall back: try to fill linearly
            for (i, v) in enumerate((Δw, Δx, Δy, Δz))
                i <= length(Δarg) && (Δarg[i] = v)
            end
        end
        return NoTangent(), Δarg
    end
    return QT(arg), AbstractQuaternion_pullback
end
function rrule(::Type{QT}, w::AbstractQuaternion) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), unthunk(Δquat))
    return QT(w), Quaternion_pullback
end
function rrule(::Type{QT}, w, x, y, z) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), components(unthunk(Δquat))...)
    return QT(SVector{4}(w, x, y, z)), Quaternion_pullback
end
function rrule(::Type{QT}, x, y, z) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), vec(unthunk(Δquat))...)
    return QT(SVector{4}(false, x, y, z)), Quaternion_pullback
end
function rrule(::Type{QT}, w::Number) where {QT<:AbstractQuaternion}
    Quaternion_pullback(Δquat) = (NoTangent(), real(unthunk(Δquat)))
    return QT(SVector{4}(w, false, false, false)), Quaternion_pullback
end

rrule(::typeof(quaternion), args...) = rrule(Quaternion{promote_type(typeof.(args)...)}, args...)
function rrule(::typeof(quaternion), v::AbstractVector)
    if length(v) == 4
        Q, Quaternion_pullback1 = rrule(quaternion, v[begin], v[begin+1], v[begin+2], v[begin+3])
    elseif length(v) == 3
        Q, Quaternion_pullback1 = rrule(quaternion, v[begin], v[begin+1], v[begin+2])
    elseif length(v) == 1
        Q, Quaternion_pullback1 = rrule(quaternion, v[begin])
    else
        throw(DimensionMismatch("Input vector must have 1, 3, or 4 inputs"))
    end
    function Quaternion_pullback2(ΔQ)
        _, Q̄... = Quaternion_pullback1(ΔQ)
        Q̄′ = similar(v)
        copyto!(Q̄′, Q̄)
        (NoTangent(), Q̄′)
    end
    return Q, Quaternion_pullback2
end
rrule(::typeof(quaternion), w::AbstractQuaternion{T}) where {T} = rrule(Quaternion{T}, w)
rrule(::typeof(quaternion), w::T) where {T<:Number} = rrule(Quaternion{T}, w)

rrule(::Type{Quaternion}, args...) = rrule(quaternion, args...)
rrule(::Type{Quaternion}, w::AbstractVector) = rrule(quaternion, w)
rrule(::Type{Quaternion}, w::AbstractQuaternion) = rrule(quaternion, w)
rrule(::Type{Quaternion}, w::Number) = rrule(quaternion, w)

# Rotor
function rrule(::typeof(rotor), w, x, y, z)
    n = √(w^2 + x^2 + y^2 + z^2)
    function Rotor_pullback(ΔR)
        # s = w/n
        # t = x/n
        # u = y/n
        # v = z/n
        ∂s∂w = (n+w)*(n-w)/n^3  # 1/n - w^2/n^3
        ∂s∂x = -x*w/n^3
        ∂s∂y = -y*w/n^3
        ∂s∂z = -z*w/n^3
        ∂t∂w = -x*w/n^3
        ∂t∂x = (n+x)*(n-x)/n^3  # 1/n - x^2/n^3
        ∂t∂y = -y*x/n^3
        ∂t∂z = -z*x/n^3
        ∂u∂w = -y*w/n^3
        ∂u∂x = -x*y/n^3
        ∂u∂y = (n+y)*(n-y)/n^3  # 1/n - y^2/n^3
        ∂u∂z = -z*y/n^3
        ∂v∂w = -z*w/n^3
        ∂v∂x = -x*z/n^3
        ∂v∂y = -y*z/n^3
        ∂v∂z = (n+z)*(n-z)/n^3  # 1/n - z^2/n^3
        Δs,Δt,Δu,Δv = components(unthunk(ΔR))
        (
            NoTangent(),
            (∂s∂w*Δs + ∂t∂w*Δt + ∂u∂w*Δu + ∂v∂w*Δv),
            (∂s∂x*Δs + ∂t∂x*Δt + ∂u∂x*Δu + ∂v∂x*Δv),
            (∂s∂y*Δs + ∂t∂y*Δt + ∂u∂y*Δu + ∂v∂y*Δv),
            (∂s∂z*Δs + ∂t∂z*Δt + ∂u∂z*Δu + ∂v∂z*Δv)
            # (1/n - w^2/n^3) + 𝐢*w*x/n^3 + 𝐣*w*y/n^3 + 𝐤*w*z/n^3,
            # - x*w/n^3 + 𝐢*(-1/n + x*x/n^3) + 𝐣*x*y/n^3 + 𝐤*x*z/n^3,
            # - y*w/n^3 + 𝐢*y*x/n^3 + 𝐣*(-1/n + y*y/n^3) + 𝐤*y*z/n^3,
            # - z*w/n^3 + 𝐢*z*x/n^3 + 𝐣*z*y/n^3 + 𝐤*(-1/n + z*z/n^3),
        )
    end
    a = SVector{4}(w, x, y, z)
    v = a ./ abs(Quaternion{eltype(a)}(a))
    return Rotor{eltype(v)}(v), Rotor_pullback
end

function rrule(::typeof(rotor), x, y, z)
    n = √(x^2 + y^2 + z^2)
    function Rotor_pullback(ΔR)
        # s = 0
        # t = x/n
        # u = y/n
        # v = z/n
        ∂t∂x = (n+x)*(n-x)/n^3  # 1/n - x^2/n^3
        ∂t∂y = -y*x/n^3
        ∂t∂z = -z*x/n^3
        ∂u∂x = -x*y/n^3
        ∂u∂y = (n+y)*(n-y)/n^3  # 1/n - y^2/n^3
        ∂u∂z = -z*y/n^3
        ∂v∂x = -x*z/n^3
        ∂v∂y = -y*z/n^3
        ∂v∂z = (n+z)*(n-z)/n^3  # 1/n - z^2/n^3
        Δt,Δu,Δv = vec(unthunk(ΔR))
        (
            NoTangent(),
            (∂t∂x*Δt + ∂u∂x*Δu + ∂v∂x*Δv),
            (∂t∂y*Δt + ∂u∂y*Δu + ∂v∂y*Δv),
            (∂t∂z*Δt + ∂u∂z*Δu + ∂v∂z*Δv)
        )
    end
    a = SVector{4}(false, x, y, z)
    v = a ./ abs(Quaternion{eltype(a)}(a))
    return Rotor{eltype(v)}(v), Rotor_pullback
end

function rrule(::typeof(rotor), v::AbstractVector)
    if length(v) == 4
        R, Rotor_pullback1 = rrule(rotor, v[begin], v[begin+1], v[begin+2], v[begin+3])
    elseif length(v) == 3
        R, Rotor_pullback1 = rrule(rotor, v[begin], v[begin+1], v[begin+2])
    elseif length(v) == 1
        R, Rotor_pullback1 = rrule(rotor, v[begin])
    else
        throw(DimensionMismatch("Input vector must have 1, 3, or 4 inputs"))
    end
    function Rotor_pullback2(ΔR)
        _, R̄... = Rotor_pullback1(ΔR)
        R̄′ = similar(v)
        copyto!(R̄′, R̄)
        (NoTangent(), R̄′)
    end
    return R, Rotor_pullback2
end

function rrule(::typeof(rotor), q::AbstractQuaternion)
    R, Rotor_pullback1 = rrule(rotor, q[1], q[2], q[3], q[4])
    function Rotor_pullback2(ΔR)
        nt, R̄w, R̄x, R̄y, R̄z = Rotor_pullback1(ΔR)
        (nt, typeof(q)(R̄w, R̄x, R̄y, R̄z))
    end
    return R, Rotor_pullback2
end

function rrule(::typeof(rotor), w::Number)
    n = √(w^2)
    function Rotor_pullback(ΔR)
        # s = w/n
        # t = 0
        # u = 0
        # v = 0
        ∂s∂w = (n+w)*(n-w)/n^3  # 1/n - w^2/n^3
        Δs = real(unthunk(ΔR))
        (
            NoTangent(),
            ∂s∂w*Δs
        )
    end
    v = SVector{4}(copysign(one(w), w), false, false, false)
    return Rotor{eltype(v)}(v), Rotor_pullback
end

rrule(::Type{Rotor}, args...) = rrule(rotor, args...)
rrule(::Type{Rotor}, w::AbstractVector) = rrule(rotor, w)
rrule(::Type{Rotor}, w::AbstractQuaternion) = rrule(rotor, w)
rrule(::Type{Rotor}, w::Number) = rrule(rotor, w)

# QuatVec
function rrule(::typeof(quatvec), w, x, y, z)
    function QuatVec_pullback(ΔV)
        (NoTangent(), ZeroTangent(), vec(unthunk(ΔV))...)
    end
    v = SVector{4}(false, x, y, z)
    return QuatVec{eltype(v)}(v), QuatVec_pullback
end

function rrule(::typeof(quatvec), x, y, z)
    function QuatVec_pullback(ΔV)
        (NoTangent(), vec(unthunk(ΔV))...)
    end
    v = SVector{4}(false, x, y, z)
    return QuatVec{eltype(v)}(v), QuatVec_pullback
end

function rrule(::typeof(quatvec), v::AbstractVector)
    if length(v) == 4
        V, QuatVec_pullback1 = rrule(quatvec, v[begin], v[begin+1], v[begin+2], v[begin+3])
    elseif length(v) == 3
        V, QuatVec_pullback1 = rrule(quatvec, v[begin], v[begin+1], v[begin+2])
    elseif length(v) == 1
        V, QuatVec_pullback1 = rrule(quatvec, v[begin])
    else
        throw(DimensionMismatch("Input vector must have 1, 3, or 4 inputs"))
    end
    function QuatVec_pullback2(ΔV)
        _, Q̄... = QuatVec_pullback1(ΔV)
        Q̄′ = similar(v)
        copyto!(Q̄′, Q̄)
        (NoTangent(), Q̄′)
    end
    return V, QuatVec_pullback2
end

function rrule(::typeof(quatvec), q::AbstractQuaternion)
    R, QuatVec_pullback1 = rrule(quatvec, q[1], q[2], q[3], q[4])
    function QuatVec_pullback2(ΔV)
        nt, Q̄w, Q̄x, Q̄y, Q̄z = QuatVec_pullback1(ΔV)
        (nt, typeof(q)(Q̄w, Q̄x, Q̄y, Q̄z))
    end
    return R, QuatVec_pullback2
end

function rrule(::typeof(quatvec), w::Number)
    function QuatVec_pullback(ΔV)
        (NoTangent(), ZeroTangent())
    end
    v = SVector{4}(w, false, false, false)
    return QuatVec{eltype(v)}(v), QuatVec_pullback
end

rrule(::Type{QuatVec}, args...) = rrule(quatvec, args...)
rrule(::Type{QuatVec}, w::AbstractQuaternion) = rrule(quatvec, w)
rrule(::Type{QuatVec}, w::AbstractVector) = rrule(quatvec, w)
rrule(::Type{QuatVec}, w::Number) = rrule(quatvec, w)



# rrule(config::RuleConfig{>:HasReverseMode}, ::Type{Rotor}, args...) = rrule_via_ad(config, rotor, args...)
# rrule(config::RuleConfig{>:HasReverseMode}, ::Type{QuatVec}, args...) = rrule_via_ad(config, quatvec, args...)






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
# COV_EXCL_START
# function (::ProjectTo{QT})(dx::AbstractFloat) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
#     return convert(QT, dx)
# end
function (::ProjectTo{QT})(dx::AbstractQuaternion{<:Integer}) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
    #@info "ProjectTo{QT}(dx::AbstractQuaternion{<:Integer})"
    return convert(QT, dx)
end
# function (::ProjectTo{QT})(dx::Integer) where {T<:AbstractFloat, QT<:AbstractQuaternion{T}}
#     return convert(QT, dx)
# end
# function (project::ProjectTo{QT})(dx::Real) where {QT<:AbstractQuaternion}
#     return project(QT(dx))
# end
# COV_EXCL_STOP

# function (project::ProjectTo{<:Number})(dx::Tangent{QT}) where {QT<:AbstractQuaternion}
#     project(QT(dx[:components]))
# end

# # ProjectTo for Number must NOT manufacture a Quaternion from a Tangent{Quaternion}.
# # Instead, take the appropriate scalar (here: the scalar component).
# function (project::ProjectTo{T})(dx::Tangent{QT}) where {QT<:AbstractQuaternion, T<:Number}
#     # dx[:components] is a 4-tuple (Δw, Δx, Δy, Δz); a Number primal wants a Number tangent.
#     return project(dx[:components][1])
# end


## Copied from `Complex` entries in ChainRulesCore.jl/src/tangent_types/abstract_zero.jl
for pattern ∈ 1:15
    T1 = iszero(pattern & 1) ? Number : AbstractZero
    T2 = iszero(pattern & 2) ? Number : AbstractZero
    T3 = iszero(pattern & 4) ? Number : AbstractZero
    T4 = iszero(pattern & 8) ? Number : AbstractZero
    w = iszero(pattern & 1) ? :w : false
    x = iszero(pattern & 2) ? :x : false
    y = iszero(pattern & 4) ? :y : false
    z = iszero(pattern & 8) ? :z : false
    @eval (QT::Type{Quaternion})(w::$T1, x::$T2, y::$T3, z::$T4) = QT($w, $x, $y, $z)  # COV_EXCL_LINE
end


## Copied from `Complex` entries in ChainRulesCore.jl/src/tangent_types/thunks.jl
# COV_EXCL_START
function (::Type{QT})(a::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a))
end
function (::Type{QT})(a::AbstractThunk, b::AbstractThunk, c::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a, b, c))
end
function (::Type{QT})(a::AbstractThunk, b::AbstractThunk, c::AbstractThunk, d::AbstractThunk) where {QT<:AbstractQuaternion}
    QT(unthunk(a, b, c, d))
end
# COV_EXCL_STOP


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

function rrule(::typeof(+), t::Real, q::AbstractQuaternion)
    y = t + q
    function add_pullback_tq(Δy)
        Δw, Δx, Δu, Δv = components(unthunk(Δy))
        return NoTangent(), Δw, typeof(q)(Δw, Δx, Δu, Δv)
    end
    return y, add_pullback_tq
end
function rrule(::typeof(+), q::AbstractQuaternion, t::Real)
    y = q + t
    function add_pullback_qt(Δy)
        Δw, Δx, Δu, Δv = components(unthunk(Δy))
        return NoTangent(), typeof(q)(Δw, Δx, Δu, Δv), Δw
    end
    return y, add_pullback_qt
end

function rrule(::typeof(*), t::Real, q::AbstractQuaternion)
    function mul_pullback(Δq)
        ∂t = @thunk q ⋅ unthunk(Δq)
        ∂q = @thunk t * unthunk(Δq)
        return (NoTangent(), ∂t, ∂q)
    end
    return t * q, mul_pullback
end
function rrule(::typeof(*), q::AbstractQuaternion, t::Real)
    function mul_pullback(Δq)
        ∂t = @thunk q ⋅ unthunk(Δq)
        ∂q = @thunk t * unthunk(Δq)
        return (NoTangent(), ∂q, ∂t)
    end
    return q * t, mul_pullback
end

function rrule(::typeof(exp), q::Quaternion{T}) where T
    w, x, y, z = components(q)
    a = absvec(q)
    e = exp(w)
    sinc = _sincu(a)
    coss = _cossu(a)

    s = e * cos(a)
    t = e * x * sinc
    u = e * y * sinc
    v = e * z * sinc
    R = quaternion(s, t, u, v)

    ∂sinc∂x = coss * x
    ∂sinc∂y = coss * y
    ∂sinc∂z = coss * z
    ∂s∂w = s
    ∂t∂w = t
    ∂u∂w = u
    ∂v∂w = v
    ∂s∂x = -e * x * sinc
    ∂s∂y = -e * y * sinc
    ∂s∂z = -e * z * sinc
    ∂t∂x = e * sinc + e * x * ∂sinc∂x
    ∂t∂y = e * x * ∂sinc∂y
    ∂t∂z = e * x * ∂sinc∂z
    ∂u∂x = e * y * ∂sinc∂x
    ∂u∂y = e * sinc + e * y * ∂sinc∂y
    ∂u∂z = e * y * ∂sinc∂z
    ∂v∂x = e * z * ∂sinc∂x
    ∂v∂y = e * z * ∂sinc∂y
    ∂v∂z = e * sinc + e * z * ∂sinc∂z

    function exp_pullback(ΔR)
        Δs, Δt, Δu, Δv = components(unthunk(ΔR))
        return (
            NoTangent(),
            typeof(q)(
                (∂s∂w * Δs + ∂t∂w * Δt + ∂u∂w * Δu + ∂v∂w * Δv),
                (∂s∂x * Δs + ∂t∂x * Δt + ∂u∂x * Δu + ∂v∂x * Δv),
                (∂s∂y * Δs + ∂t∂y * Δt + ∂u∂y * Δu + ∂v∂y * Δv),
                (∂s∂z * Δs + ∂t∂z * Δt + ∂u∂z * Δu + ∂v∂z * Δv)
            )
        )
    end

    return R, exp_pullback
end

function rrule(::typeof(exp), v⃗::QuatVec{T}) where T
    x, y, z = vec(v⃗)
    a = absvec(v⃗)
    sinc = _sincu(a)
    coss = _cossu(a)

    s = cos(a)
    t = x * sinc
    u = y * sinc
    v = z * sinc
    R = rotor(s, t, u, v)

    ∂sinc∂x = coss * x
    ∂sinc∂y = coss * y
    ∂sinc∂z = coss * z
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
            typeof(v⃗)(
                (∂s∂x * Δs + ∂t∂x * Δt + ∂u∂x * Δu + ∂v∂x * Δv),
                (∂s∂y * Δs + ∂t∂y * Δt + ∂u∂y * Δu + ∂v∂y * Δv),
                (∂s∂z * Δs + ∂t∂z * Δt + ∂u∂z * Δu + ∂v∂z * Δv)
            )
        )
    end

    return R, exp_pullback
end


end
