module QuaternionicChainRulesCoreExt

import Quaternionic
isdefined(Base, :get_extension) ?
    (import ChainRulesCore; import Zygote; import ZygoteRules) :
    (import ...ChainRulesCore; import ...Zygote; import ...ZygoteRules)


# Zygote.jl/src/lib/number.jl
function ChainRulesCore.rrule(::ZygoteRuleConfig, T::Type{<:AbstractQuaternion}, w, x, y, z)
    Quaternion_pullback(c̄) = (NoTangent(), c̄[1], c̄[2], c̄[3], c̄[4])
    return T(w, x, y, z), Quaternion_pullback
end
# # for real x, ChainRules pulls back a zero real adjoint, whereas we treat x
# # as embedded in the complex numbers and pull back a pure imaginary adjoint
# function ChainRulesCore.rrule(::ZygoteRuleConfig, ::typeof(imag), x::Number)
#     imag_pullback(ī) = (NoTangent(), real(ī)*im)
#     return imag(x), imag_pullback
# end

# Zygote.jl/src/lib/broadcast.jl
Zygote._dual_safearg(x::Numeric{<:AbstractQuaternion}) = true
Zygote._dual_safearg(x::Ref{<:Numeric{<:AbstractQuaternion}}) = true
@inline function Zygote.dual(x::Q, i, ::Val{N}) where {T<:Number, Q<:AbstractQuaternion{T}, N}
    w_dual = Dual(x[1], ntuple(==(i), 4N))
    x_dual = Dual(x[2], ntuple(==(N+i), 4N))
    y_dual = Dual(x[3], ntuple(==(2N+i), 4N))
    z_dual = Dual(x[4], ntuple(==(3N+i), 4N))
    return Q(w_dual, x_dual, y_dual, z_dual)
end

# Zygote.jl/src/lib/array.jl
# @adjoint real(x::AbstractArray) = real(x), r̄ -> (real(r̄),)
# @adjoint conj(x::AbstractArray) = conj(x), r̄ -> (conj(r̄),)
# @adjoint imag(x::AbstractArray) = imag(x), ī -> (complex.(0, real.(ī)),)

# Zygote.jl/src/lib/grad.jl
Zygote._jvec(x::AbstractArray{<:AbstractQuaternion}) = throw(ArgumentError("jacobian does not accept quaternionic output"))

# Zygote.jl/src/compiler/interface.jl
Zygote.sensitivity(y::AbstractQuaternion) = error("Output is quaternionic, so the gradient is not defined.")

# Zygote.jl/src/forward/number.jl
DiffRules._abs_deriv(x::AbstractQuaternion) = x/abs(x)


end # module
