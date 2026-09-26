module QuaternionicForwardDiffExt

using Quaternionic
isdefined(Base, :get_extension) ? (using ForwardDiff) : (using ..ForwardDiff)

# Recurse so that nested duals, as in higher-order derivatives, are stripped to the
# innermost value.  Since ForwardDiff 1.0, `iszero(::Dual)` also checks the partials, so
# stripping only one level would make `iszerovalue` false whenever an inner partial is
# nonzero.  See issue #113.
Quaternionic.value(x::ForwardDiff.Dual) = Quaternionic.value(ForwardDiff.value(x))
# Then, this will be automatic:
# Quaternionic.iszerovalue(x::ForwardDiff.Dual) = iszero(Quaternionic.value(x))

###############################################
# Let ForwardDiff act naturally on quaternions.

# These next three definitions are cribbed from similar expressions
# enabling differentiation of complex-valued functions in
# https://github.com/JuliaDiff/ForwardDiff.jl/blob/78c73afd9a21593daf54f61c7d0db67130cf29e1/src/derivative.jl#L83-L88

@inline ForwardDiff.extract_derivative(::Type{T}, y::AbstractQuaternion) where {T} = zero(y)

# Both Quaternion and Rotor, when differentiated, result in a Quaternion...
@inline function ForwardDiff.extract_derivative(::Type{T}, y::AbstractQuaternion{TD}) where {T, TD <: ForwardDiff.Dual}
    quaternion(
        ForwardDiff.partials(T, y[1], 1),
        ForwardDiff.partials(T, y[2], 1),
        ForwardDiff.partials(T, y[3], 1),
        ForwardDiff.partials(T, y[4], 1)
    )
end

# ...but QuatVec results in a QuatVec
# COV_EXCL_START
@inline function ForwardDiff.extract_derivative(::Type{T}, y::QuatVec{TD}) where {T, TD <: ForwardDiff.Dual}
    quatvec(
        ForwardDiff.partials(T, y[2], 1),
        ForwardDiff.partials(T, y[3], 1),
        ForwardDiff.partials(T, y[4], 1)
    )
end
# COV_EXCL_STOP

end # module
