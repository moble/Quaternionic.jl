# Zygote-specific glue for quaternions.  Every derivative rule lives in
# `QuaternionicChainRulesCoreExt`, and Zygote uses those rules directly.  This module adds only
# what ChainRules cannot express.  Zygote prefers its own rules for `ZygoteRuleConfig` and its
# own `@adjoint`s over config-free ChainRules rules, and several of those assume commutative
# multiplication or complex-number semantics, because `AbstractQuaternion <: Number`.  They
# are overridden here for quaternion arguments.
#
# This extension is loaded when both Zygote and ChainRulesCore are loaded.
module QuaternionicZygoteExt

import Quaternionic: Quaternionic, AbstractQuaternion
import Zygote
import Zygote: ZygoteRuleConfig, @adjoint
import ChainRulesCore: rrule, rrule_via_ad, ProjectTo, NoTangent, unthunk

const AQ = AbstractQuaternion
const QArr = AbstractArray{<:AbstractQuaternion}
const NArr = AbstractArray{<:Number}

# Zygote defines rules for these functions on `(::ZygoteRuleConfig, ::Number)`, which outrank
# the config-free quaternion rules.  Their pullbacks treat the argument as a complex number:
# `imag`'s pullback returns `real(Δ) * im`, for example, which turns the cotangent of a
# `QuatVec` into a complex vector.  These methods are more specific, and they forward to the
# config-free rules.  (Defining the quaternion rules on `RuleConfig` instead would be
# ambiguous with Zygote's.)
for f ∈ (:real, :imag, :conj, :abs2, :abs)
    @eval rrule(::ZygoteRuleConfig, ::typeof($f), q::AQ) = rrule($f, q)
end

# Zygote's rule for `literal_pow` (which `q^2` with a literal exponent calls) differentiates
# `x^p` as if multiplication commuted.  This method forwards to the quaternion rule.  Where
# there is none, because the power is opted out (as for complex `Rotor`s), it differentiates
# `q^p` instead, which reaches the rules or the source of `^`.
function rrule(
    config::ZygoteRuleConfig, ::typeof(Base.literal_pow), ::typeof(^), q::AQ, v::Val{p}
) where {p}
    result = rrule(Base.literal_pow, ^, q, v)
    result === nothing || return result
    Ω, power_pullback = rrule_via_ad(config, ^, q, p)
    literal_pow_pullback(ΔΩ) = (NoTangent(), NoTangent(), power_pullback(ΔΩ)[2], NoTangent())
    return Ω, literal_pow_pullback
end

# Zygote's rule for `+(xs::Number...)` passes the cotangent of the sum to every argument
# without projecting it, which gives a real argument a quaternion gradient.  This function
# returns the sum and a pullback that first projects the cotangent onto the sum (which
# normalizes the forms that Zygote may pass, such as a `NamedTuple`) and then onto each
# argument.  A real argument thus receives the scalar part of the cotangent.
function projected_plus(args...)
    Ω = +(args...)
    project_sum = ProjectTo(Ω)
    projections = map(ProjectTo, args)
    function plus_pullback(ΔΩ)
        Δ = project_sum(unthunk(ΔΩ))
        return (NoTangent(), map(p -> p(Δ), projections)...)
    end
    return Ω, plus_pullback
end

# These methods cover sums of one or two arguments, at least one of which is a quaternion.
rrule(::ZygoteRuleConfig, ::typeof(+), a::AQ) = projected_plus(a)
for (A, B) ∈ ((:AQ, :AQ), (:AQ, :Number), (:Number, :AQ))
    @eval rrule(::ZygoteRuleConfig, ::typeof(+), a::$A, b::$B) = projected_plus(a, b)
end

# Sums of three or more arguments with a quaternion among the first three.  The seven
# patterns of quaternions and numbers in those positions cover each other's intersections, so
# that no two of these methods are ambiguous.
for (A, B, C) ∈ Iterators.product((:AQ, :Number), (:AQ, :Number), (:AQ, :Number))
    :AQ ∈ (A, B, C) || continue
    @eval rrule(::ZygoteRuleConfig, ::typeof(+), a::$A, b::$B, c::$C, xs::Number...) =
        projected_plus(a, b, c, xs...)
end

# Sums whose first quaternion is in position 4 through 8, after real or complex numbers.  A
# first quaternion further along still gets Zygote's unprojected rule, which gives the other
# arguments quaternion gradients with the correct scalar parts.  These methods cannot be
# ambiguous with those above, because a quaternion is not a `Union{Real,Complex}`.
for k ∈ 4:8
    reals = [Symbol(:r, i) for i ∈ 1:k-1]
    signature = [:($(r)::Union{Real,Complex}) for r ∈ reals]
    @eval rrule(::ZygoteRuleConfig, ::typeof(+), $(signature...), q::AQ, xs::Number...) =
        projected_plus($(reals...), q, xs...)
end

# Zygote's `@adjoint inv(::Union{Number,AbstractMatrix})` bypasses ChainRules, and its formula
# is wrong for quaternions with complex components.  This adjoint calls the quaternion rule.
@adjoint function Base.inv(q::AQ)
    Ω, back = rrule(inv, q)
    function inv_pullback(Δ)
        return (Zygote.wrap_chainrules_output(back(Zygote.wrap_chainrules_input(Δ))[2]),)
    end
    return Ω, inv_pullback
end

# Zygote's broadcast adjoints for `.*`, `./`, `.^p`, `abs2.`, and `imag.` assume commutative
# multiplication or complex arguments, and scalar times quaternion array overflows the stack.
# When any operand is a quaternion or an array of quaternions, these adjoints use Zygote's
# generic broadcast instead, which applies the scalar rules elementwise.  The pullback of
# `Zygote._broadcast_generic` returns a leading cotangent for the broadcast style, which is
# dropped here.  The signatures with `Bool` resolve ambiguities with Zygote's adjoints for
# `Bool` operands.  The broadcast adjoints for `.+`, `.-`, `real.`, and `conj.` are already
# correct for quaternions.
const BinarySignatures = (
    (:AQ, :QArr), (:QArr, :AQ), (:QArr, :QArr), (:Number, :QArr), (:QArr, :Number),
    (:Bool, :QArr), (:QArr, :Bool), (:AQ, :NArr), (:NArr, :AQ), (:QArr, :NArr), (:NArr, :QArr),
)
for op ∈ (:*, :/), (A, B) ∈ BinarySignatures
    @eval @adjoint function Base.Broadcast.broadcasted(::typeof($op), x::$A, y::$B)
        Ω, back = Zygote._broadcast_generic(__context__, $op, x, y)
        return Ω, Base.tail ∘ back
    end
end
@adjoint function Base.Broadcast.broadcasted(
    ::typeof(Base.literal_pow), ::typeof(^), x::QArr, p::Val
)
    Ω, back = Zygote._broadcast_generic(__context__, y -> Base.literal_pow(^, y, p), x)
    return Ω, Δ -> (nothing, nothing, back(Δ)[3], nothing)
end
for f ∈ (:abs2, :imag)
    @eval @adjoint function Base.Broadcast.broadcasted(::typeof($f), x::QArr)
        Ω, back = Zygote._broadcast_generic(__context__, $f, x)
        return Ω, Base.tail ∘ back
    end
end

# A broadcast whose operands are all scalars computes the same thing as the scalar call, but
# Zygote's adjoints for scalar operands of `.*`, `./`, `.^p`, `abs2.`, and `imag.` share the
# commutative and complex assumptions of the array adjoints.  Zygote's generic broadcast does
# not handle these zero-dimensional broadcasts, so when a quaternion is among the operands,
# these adjoints differentiate the scalar call instead, which reaches the quaternion rules.
# A broadcast of `*` over three or more scalar operands, as in `broadcast(*, a, b, c)`, is
# handled in the same way when a quaternion is among the first three operands; the seven
# patterns are those of the sums above.  (Chained operators, as in `a .* b .* c`, are nested
# broadcasts of two operands.)
for op ∈ (:*, :/), (A, B) ∈ ((:AQ, :AQ), (:AQ, :Number), (:Number, :AQ))
    @eval @adjoint Base.Broadcast.broadcasted(::typeof($op), x::$A, y::$B) =
        Zygote._pullback(__context__, $op, x, y)
end
for (A, B, C) ∈ Iterators.product((:AQ, :Number), (:AQ, :Number), (:AQ, :Number))
    :AQ ∈ (A, B, C) || continue
    @eval @adjoint Base.Broadcast.broadcasted(
        ::typeof(*), x::$A, y::$B, z::$C, zs::Number...
    ) = Zygote._pullback(__context__, *, x, y, z, zs...)
end
@adjoint Base.Broadcast.broadcasted(::typeof(Base.literal_pow), ::typeof(^), x::AQ, p::Val) =
    Zygote._pullback(__context__, Base.literal_pow, ^, x, p)
for f ∈ (:abs2, :imag)
    @eval @adjoint Base.Broadcast.broadcasted(::typeof($f), x::AQ) =
        Zygote._pullback(__context__, $f, x)
end

# `Zygote.gradient` seeds the output of the function with `sensitivity(y)`, which is
# `one(y)` for a `Number`, and `Zygote.jacobian` treats a `Number` output, or each element of
# an array output, as a single scalar.  For quaternion outputs, both would silently return
# the derivatives of the scalar parts only.  Zygote raises errors for complex outputs in the
# same places, and these methods do the same for quaternions.  Quaternions inside the
# function are unaffected; only the final output matters.
function Zygote.sensitivity(::AQ)
    throw(ArgumentError(
        "Output is a quaternion, so the gradient is not defined.  Differentiate a real " *
        "number computed from it, call the pullback from `Zygote.pullback` with an explicit " *
        "quaternion cotangent, or use `Zygote.jacobian` on `x -> collect(components(f(x)))`."
    ))
end
Zygote._jvec(::AQ) = throw(quaternion_jacobian_error())
Zygote._jvec(::QArr) = throw(quaternion_jacobian_error())
function quaternion_jacobian_error()
    ArgumentError(
        "jacobian does not accept quaternion output.  Return the components as a real " *
        "array instead, for example with `collect(components(q))` or `to_float_array`."
    )
end

end # module
