module QuaternionicForwardDiffExt

using Quaternionic
import Quaternionic: AbstractQuaternion, Quaternion, QuatVec, quaternion, quatvec, components,
    wrapper, basetype
# Under Requires (Julia < 1.9), ForwardDiff must be reached relative to the parent module;
# an absolute `using ForwardDiff` works there but warns that Quaternionic does not depend
# on ForwardDiff.
isdefined(Base, :get_extension) ? (using ForwardDiff) : (using ..ForwardDiff)

# Recurse so that nested duals, as in higher-order derivatives, are stripped to the
# innermost value.  Since ForwardDiff 1.0, `iszero(::Dual)` also checks the partials, so
# stripping only one level would make `iszerovalue` false whenever an inner partial is
# nonzero.  See issue #113.
Quaternionic.value(x::ForwardDiff.Dual) = Quaternionic.value(ForwardDiff.value(x))
# Then, this will be automatic:
# Quaternionic.iszerovalue(x::ForwardDiff.Dual) = iszero(Quaternionic.value(x))

###############################################
# Extract values and derivatives from quaternion-valued results.
#
# `ForwardDiff.derivative` calls `extract_derivative(T, y)` on the result `y`, and
# DifferentiationInterface calls `value(T, y)` and `partials(T, y, i)` on it, where `T` is
# the tag of the differentiation in progress.  For a `Number` that is not a `Dual`,
# ForwardDiff's fallbacks return `zero(y)`, `y`, and `zero(y)`, respectively, which would
# silently discard the derivatives stored in the components of a quaternion.  So we apply
# ForwardDiff's own methods to each component, which handle a `Dual` of any tag and a
# constant.  Complex components, as in `Lorentz` rotors, are split into their real and
# imaginary parts, because only recent versions of ForwardDiff provide `extract_derivative`
# for `Complex{<:Dual}`, and none provides `value` or `partials` for it.  Passing the tag
# `T` through keeps nested differentiation free of perturbation confusion.
#
# The derivative of a `QuatVec` is a `QuatVec`, and the derivative of a `Quaternion` or a
# `Rotor` is a `Quaternion`; derivatives of rotors are not rotors.  The value, on the other
# hand, keeps the type of the result: a `Rotor` stays a `Rotor`, and its components are
# stored as they are, without normalization.
#
# `ForwardDiff.can_dual` is deliberately not extended to quaternions: inputs to ForwardDiff
# must be real numbers or arrays of them.

# Apply the ForwardDiff extraction function `f` to the component `c`.
componentwise(f, ::Type{T}, c, args...) where {T} = f(T, c, args...)
function componentwise(f, ::Type{T}, c::Complex, args...) where {T}
    complex(f(T, real(c), args...), f(T, imag(c), args...))
end

# Apply `f` to each component of `y`, and return the results as a `Quaternion`, or as a
# `QuatVec` if `y` is a `QuatVec`.
function componentwise(f, ::Type{T}, y::AbstractQuaternion, args...) where {T}
    rewrap(y, map(c -> componentwise(f, T, c, args...), Tuple(components(y))))
end
rewrap(::AbstractQuaternion, c) = quaternion(c...)
rewrap(::QuatVec, c) = quatvec(c...)

function ForwardDiff.extract_derivative(::Type{T}, y::AbstractQuaternion) where {T}
    componentwise(ForwardDiff.extract_derivative, T, y)
end

function ForwardDiff.value(::Type{T}, y::AbstractQuaternion) where {T}
    c = map(cᵢ -> componentwise(ForwardDiff.value, T, cᵢ), Tuple(components(y)))
    wrapper(y){promote_type(map(typeof, c)...)}(c...)
end

function ForwardDiff.partials(::Type{T}, y::AbstractQuaternion, i) where {T}
    componentwise(ForwardDiff.partials, T, y, i)
end

# `ForwardDiff.jacobian` of a function that returns an array of quaternions allocates its
# result with the element type `valtype(T, eltype(y))`, and fills it with `partials` above.
# For a quaternion of duals, ForwardDiff's fallback would return the quaternion type itself,
# so the Jacobian would silently hold dual numbers instead of derivatives.  This method gives
# the type that `partials` returns, so that each entry of the Jacobian is the derivative of
# one quaternion of the output with respect to one real input.
function ForwardDiff.valtype(::Type{T}, ::Type{Q}) where {T,Q<:AbstractQuaternion}
    derivativewrapper(Q){partialstype(T, basetype(Q))}
end
partialstype(::Type{T}, ::Type{C}) where {T,C} = ForwardDiff.valtype(T, C)
partialstype(::Type{T}, ::Type{Complex{C}}) where {T,C} = Complex{ForwardDiff.valtype(T, C)}
derivativewrapper(::Type{<:AbstractQuaternion}) = Quaternion
derivativewrapper(::Type{<:QuatVec}) = QuatVec

end # module
