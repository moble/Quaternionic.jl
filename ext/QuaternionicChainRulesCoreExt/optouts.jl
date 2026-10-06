# Opt-outs of ChainRules' generic `Number` rules (section 2 of the AD specification).
#
# `AbstractQuaternion <: Number`, so the scalar rules of ChainRules (and the rules that
# ChainRulesCore's `@scalar_rule` generates) apply to quaternion arguments, and they assume
# that multiplication commutes.  Where this extension supplies no rule of its own, the
# opt-outs below make the AD system differentiate the package's source instead.
#
# These opt-outs live in a file of their own because `@opt_out rrule(sig)` defines both
# `no_rrule(sig)` and `rrule(sig) = nothing`.  An opt-out with exactly the signature of a
# rule supplied in another file would therefore silently delete that rule, so each signature
# here must differ from every supplied one.  An opt-out that is broader than a supplied rule
# (such as `log(::AbstractQuaternion)` here and `log(::Quaternion{<:Real})` in
# `elementary.jl`) coexists with it: dispatch selects the narrower rule, and the AD system
# skips a rule only when the selected `rrule` and `no_rrule` methods have the same signature.
#
# The `frule` opt-outs take `::Tuple` tangents, like the supplied `frule`s; with `::Any`,
# they would be ambiguous with them.

# The real-component rules for `log` and `sqrt` are in `elementary.jl`.  These catch
# complex components (where the generic rule has a relative error of 2.3 for `log`) and any
# other element type.
for f ∈ (:log, :sqrt, :(Base.FastMath.log_fast), :(Base.FastMath.sqrt_fast))
    @eval @opt_out rrule(::typeof($f), ::AQ)
    @eval @opt_out frule(::Tuple, ::typeof($f), ::AQ)
end

# The rules for real powers (and integer powers of most quaternions) are in `elementary.jl`;
# these catch complex components and complex rotors.
for f ∈ (:^, :(Base.FastMath.pow_fast))
    for (A, B) ∈ ((:AQ, :Number), (:Number, :AQ), (:AQ, :AQ), (:(Rotor{<:Complex}), :Integer))
        @eval @opt_out rrule(::typeof($f), ::$A, ::$B)
        @eval @opt_out frule(::Tuple, ::typeof($f), ::$A, ::$B)
    end
end

# The source renormalizes the product or quotient of two rotors, at least one of which has
# complex components, only when the sum of the squared magnitudes of the raw components is
# less than 16 (`Quaternionic.complex_rotor_product`).  Tracing the source follows that test
# exactly.  Base computes `a \ b` as `adjoint(adjoint(b) / adjoint(a))`, which is such a
# quotient for two rotors, and the `Base.FastMath` versions call the ordinary operators.
for f ∈ (:*, :/, :\, :(Base.FastMath.mul_fast), :(Base.FastMath.div_fast))
    for (A, B) ∈ ((:(Rotor{<:Complex}), :Rotor), (:Rotor, :(Rotor{<:Complex})),
                  (:(Rotor{<:Complex}), :(Rotor{<:Complex})))
        @eval @opt_out rrule(::typeof($f), ::$A, ::$B)
        @eval @opt_out frule(::Tuple, ::typeof($f), ::$A, ::$B)
    end
end

# Base's `sign(x) = x / abs(x)` runs for quaternions, whereas the generic rule uses `im`.
for f ∈ (:sign, :(Base.FastMath.sign_fast))
    @eval @opt_out rrule(::typeof($f), ::AQ)
    @eval @opt_out frule(::Tuple, ::typeof($f), ::AQ)
end

# The angle of a quaternion with complex components is `2absvec(log(q))`; the rules for real
# components are in `geometry.jl`.
for f ∈ (:angle, :(Base.FastMath.angle_fast))
    @eval @opt_out rrule(::typeof($f), ::Union{Quaternion{<:Complex},Rotor{<:Complex}})
    @eval @opt_out frule(::Tuple, ::typeof($f), ::Union{Quaternion{<:Complex},Rotor{<:Complex}})
end

# `norm(q, p)` and `LinearAlgebra.norm2(q)` reach Base's `Number` methods, which are not the
# quaternion norm for complex components.  The rules for `norm(q)` are in `geometry.jl`.
@opt_out rrule(::typeof(norm), ::AQ, ::Real)
@opt_out frule(::Tuple, ::typeof(norm), ::AQ, ::Real)
@opt_out frule(::Tuple, ::typeof(LinearAlgebra.norm2), ::AQ)

# Base's `Number` methods of `det`, `logdet`, and `pinv` run for quaternions, but the generic
# rules assume commutativity.
for f ∈ (:(LinearAlgebra.det), :(LinearAlgebra.logdet), :(LinearAlgebra.pinv))
    @eval @opt_out rrule(::typeof($f), ::AQ)
    @eval @opt_out frule(::Tuple, ::typeof($f), ::AQ)
end

# Every pattern of three arguments, each a quaternion or another number, with at least one
# quaternion.  `fma` has a `@scalar_rule`, which assumes commutativity.  The `frule`s of
# ChainRules for `*` and `+` with three or more arguments project nothing, and the one for `*`
# with three arguments ignores the normalization of rotor products, so these are opted out
# too; their `rrule`s are the folds in `arithmetic.jl`.
for (A, B, C) ∈ Iterators.product(ntuple(i -> (:AQ, :Number), 3)...)
    :AQ ∈ (A, B, C) || continue
    @eval @opt_out rrule(::typeof(fma), ::$A, ::$B, ::$C)
    @eval @opt_out frule(::Tuple, ::typeof(fma), ::$A, ::$B, ::$C)
    for f ∈ (:*, :+, :(Base.FastMath.mul_fast), :(Base.FastMath.add_fast))
        @eval @opt_out frule(::Tuple, ::typeof($f), ::$A, ::$B, ::$C, ::Vararg{Number})
    end
end
