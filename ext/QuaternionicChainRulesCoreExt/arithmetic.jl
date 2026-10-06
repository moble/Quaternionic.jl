# Rules for arithmetic (section 3.2 of the AD specification): the binary operations `+`, `-`,
# `*`, `/`, and `\`, with the `Base.FastMath` versions of the first four; unary `-`, `conj`,
# and `inv`, with `sub_fast` and `inv_fast`; the left folds of `*` and `+` (and of `mul_fast`
# and `add_fast`) with three or more arguments; and `muladd`.
#
# Each rule is the exact derivative of what the source computes from the stored components:
#
#   * Quaternion multiplication is bilinear (ℂ-bilinear for complex components), so the
#     pullback of `a * b` is `(Δ * qadjoint(b), qadjoint(a) * Δ)`.
#   * Division by a `Quaternion` or a `QuatVec` is multiplication by its inverse.  Division
#     by a `Rotor` divides by `abs2(::Rotor)`, which is identically one, so `a / R` is
#     `a * conj(R)`, even for a `Rotor` that does not have unit norm.
#   * The product and the quotient of two `Rotor`s with real components are normalized
#     (through `rotor`).  Those with complex components are normalized only when the sum of
#     the squared magnitudes of the raw components is less than 16, as in
#     `Quaternionic.complex_rotor_product`.  The two-argument rules for complex `Rotor` pairs
#     are opted out in `optouts.jl`, so that the AD system traces the source; but the folds
#     and `muladd`, which cannot hand a single step back to the AD system, reproduce that test
#     in `isrenormalized`.
#   * `a \ b` is Base's generic `adjoint(adjoint(b) / adjoint(a))`, so its rule is the
#     composition of the rules for those three steps.  This gives `inv(a) * b` for a
#     `Quaternion` `a`, but `R \ b == conj(R) * b` for a `Rotor` `R`, and `t \ q == q / conj(t)`
#     for a complex scalar `t`.
#   * `muladd(a, b, c)` promotes its three arguments to a common type and then computes
#     `a * b + c`.  Promotion converts without normalizing, so it has a unit derivative; but it
#     leaves three `Rotor`s as `Rotor`s, whose product is then normalized.
#
# The rules for two arguments forward to `binary_rrule` and `binary_frule`, which dispatch on
# the operation and the argument types through `binary_cotangents` and `binary_pushforward`.

########################################################################################
## Normalizing tangents and cotangents
########################################################################################

"""
    outputcot(Ω, ΔΩ)

Normalize the cotangent `ΔΩ` of the output `Ω` of an arithmetic operation: `cot(ΔΩ)` for a
quaternion `Ω`, with the scalar part set to zero for a `QuatVec` (whose scalar part is
identically zero, so that its cotangent is meaningless), and `scalarcot(ΔΩ)` for any other
number.
"""
outputcot(::AQ, ΔΩ) = cot(ΔΩ)
function outputcot(::QuatVec, ΔΩ)
    Δ = cot(ΔΩ)
    Δ isa AbstractZero && return Δ
    return quaternion(zero(Δ[1]), Δ[2], Δ[3], Δ[4])
end
outputcot(::Number, ΔΩ) = scalarcot(ΔΩ)

"""
    iszerotangent(ẋ)

Return `true` when the tangent `ẋ` is `nothing` or an `AbstractZero`, possibly inside a
thunk.  The result depends only on the type of `ẋ`.
"""
iszerotangent(ẋ) = unthunk(ẋ) isa Union{Nothing,AbstractZero}

"""
    inputtangent(x, ẋ)

Return the tangent `ẋ` of the argument `x` in the form that the pushforwards use: a
`Quaternion` for a quaternion `x` (with zero scalar part for a `QuatVec`), a number for any
other number, and `false` (a strong zero) for a zero tangent.  A number that is the tangent
of a quaternion is the quaternion with only that scalar part, so that the tangent of a
number promoted to a quaternion (as in `muladd`) is embedded correctly.
"""
inputtangent(::AQ, ẋ) = zerofill(cot(ẋ))
function inputtangent(::QuatVec, ẋ)
    Δ = cot(ẋ)
    Δ isa AbstractZero && return false
    return quaternion(zero(Δ[1]), Δ[2], Δ[3], Δ[4])
end
inputtangent(::Number, ẋ) = zerofill(scalarcot(ẋ))

########################################################################################
## Products and quotients of rotors, which may be normalized
########################################################################################

"""
    rawproduct(op, a, b)

Return the raw product (for `op` being `*`) or quotient (for `/`) of the `Rotor`s `a` and
`b` as a `Quaternion`, before the normalization that the source applies to `Rotor`
products.  The components are computed by the same expressions as in the source.
"""
rawproduct(::typeof(*), a::Rotor, b::Rotor) = quaternion(a) * b
rawproduct(::typeof(/), a::Rotor, b::Rotor) = quaternion(a) / b

"""
    isrenormalized(P)

Return `true` when the raw product or quotient `P` of two `Rotor`s, at least one of which
has complex components, is normalized by the source.  This is the test of
`Quaternionic.complex_rotor_product`, applied to values, so that dual numbers take the same
branch as the numbers they contain.
"""
isrenormalized(P::AQ) = (value(sum(abs2, components(P))) < 16) === true

########################################################################################
## Pullbacks of the binary operations
########################################################################################

"""
    product_cotangents(a, b, Δ)

Return the cotangents `(∂a, ∂b)` of the arguments of the raw product `a * b`, given the
`Quaternion` (or, for two numbers, scalar) cotangent `Δ` of the product.  A cotangent of a
number that is not a quaternion is a number.
"""
product_cotangents(a::AQ, b::AQ, Δ) = (Δ * qadjoint(b), qadjoint(a) * Δ)
product_cotangents(a::Number, b::AQ, Δ) = (inner(b, Δ), conj(a) * Δ)
product_cotangents(a::AQ, b::Number, Δ) = (Δ * conj(b), inner(a, Δ))
product_cotangents(a::Number, b::Number, Δ) = (Δ * conj(b), conj(a) * Δ)

"""
    quotient_cotangents(Ω, a, b, Δ)

Return the cotangents `(∂a, ∂b)` of the arguments of the raw quotient `Ω = a / b`, given
the `Quaternion` cotangent `Δ` of the quotient.  For a `Quaternion` or `QuatVec` divisor,
`Ω = a * inv(b)`; for a `Rotor` divisor, the source computes `Ω = a * conj(b)`; and for a
number that is not a quaternion, `Ω` is `a` with each component divided by `b`.
"""
function quotient_cotangents(Ω, a, b::AQ, Δ)
    ib = inv(q_(b))
    ∂a = first(product_cotangents(a, ib, Δ))
    ∂b = -(qadjoint(Ω) * Δ * qadjoint(ib))
    return ∂a, ∂b
end
quotient_cotangents(Ω, a, b::Number, Δ) = (Δ / conj(b), -inner(Ω, Δ) / conj(b))
function quotient_cotangents(Ω, a, b::Rotor, Δ)
    ∂a = first(product_cotangents(a, conj(quaternion(b)), Δ))
    ∂b = conj(qadjoint(a) * Δ)
    return ∂a, ∂b
end

"""
    binary_cotangents(op, Ω, a, b, Δ)

Return the cotangents `(∂a, ∂b)` of the arguments of `Ω = op(a, b)`, before projection,
given the cotangent `Δ` of `Ω`, normalized by `outputcot`.
"""
binary_cotangents(::typeof(+), Ω, a, b, Δ) = (Δ, Δ)
binary_cotangents(::typeof(-), Ω, a, b, Δ) = (Δ, -Δ)
binary_cotangents(::typeof(*), Ω, a, b, Δ) = product_cotangents(a, b, Δ)
binary_cotangents(::typeof(/), Ω, a, b, Δ) = quotient_cotangents(Ω, a, b, Δ)
function binary_cotangents(::typeof(\), Ω, a, b, Δ)
    # Ω = adjoint(Ω′) with Ω′ = adjoint(b) / adjoint(a).  The `adjoint` of a quaternion is
    # its quaternion conjugate, and that of a complex number is its complex conjugate; each
    # is its own adjoint as a linear map over the reals, so cotangents pass through it as
    # tangents do.
    ∂b′, ∂a′ = binary_cotangents(/, adjoint(Ω), adjoint(b), adjoint(a), adjoint(Δ))
    return adjoint(∂a′), adjoint(∂b′)
end

# The normalized products and quotients of two rotors.
for op ∈ (:*, :/)
    @eval begin
        function binary_cotangents(::typeof($op), Ω, a::Rotor{<:Real}, b::Rotor{<:Real}, Δ)
            P = rawproduct($op, a, b)
            return binary_cotangents($op, P, quaternion(a), b, normalize_pullback(P, Δ))
        end
        function binary_cotangents(::typeof($op), Ω, a::Rotor, b::Rotor, Δ)
            P = rawproduct($op, a, b)
            ΔP = isrenormalized(P) ? normalize_pullback(P, Δ) : Δ
            return binary_cotangents($op, P, quaternion(a), b, ΔP)
        end
    end
end

"""
    binary_rrule(op, a, b)

Return the primal `op(a, b)` and its pullback, for `op` being `+`, `-`, `*`, `/`, or `\\`,
and for any combination of quaternions and other numbers.
"""
function binary_rrule(op, a, b)
    Ω = op(a, b)
    proja, projb = ProjectTo(a), ProjectTo(b)
    function binary_pullback(ΔΩ)
        Δ = outputcot(Ω, ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ, Δ)
        ∂a, ∂b = binary_cotangents(op, Ω, a, b, Δ)
        return NoTangent(), proja(∂a), projb(∂b)
    end
    return Ω, binary_pullback
end

########################################################################################
## Pushforwards of the binary operations
########################################################################################

"""
    binary_pushforward(op, Ω, a, b, ȧ, ḃ)

Return the tangent of `Ω = op(a, b)`, before projection, given the tangents `ȧ` and `ḃ`
of the arguments, normalized by `inputtangent`.
"""
binary_pushforward(::typeof(+), Ω, a, b, ȧ, ḃ) = ȧ + ḃ
binary_pushforward(::typeof(-), Ω, a, b, ȧ, ḃ) = ȧ - ḃ
binary_pushforward(::typeof(*), Ω, a, b, ȧ, ḃ) = ȧ * q_(b) + q_(a) * ḃ
binary_pushforward(::typeof(/), Ω, a, b, ȧ, ḃ) = (ȧ - q_(Ω) * ḃ) * inv(q_(b))
binary_pushforward(::typeof(/), Ω, a, b::Rotor, ȧ, ḃ) =
    ȧ * conj(quaternion(b)) + q_(a) * conj(ḃ)
function binary_pushforward(::typeof(\), Ω, a, b, ȧ, ḃ)
    # As in `binary_cotangents`, Ω = adjoint(adjoint(b) / adjoint(a)).
    Ω̇′ = binary_pushforward(/, adjoint(Ω), adjoint(b), adjoint(a), adjoint(ḃ), adjoint(ȧ))
    return adjoint(Ω̇′)
end

# The normalized products and quotients of two rotors.
for op ∈ (:*, :/)
    @eval begin
        function binary_pushforward(::typeof($op), Ω, a::Rotor{<:Real}, b::Rotor{<:Real}, ȧ, ḃ)
            P = rawproduct($op, a, b)
            return normalize_pushforward(P, binary_pushforward($op, P, quaternion(a), b, ȧ, ḃ))
        end
        function binary_pushforward(::typeof($op), Ω, a::Rotor, b::Rotor, ȧ, ḃ)
            P = rawproduct($op, a, b)
            Ṗ = binary_pushforward($op, P, quaternion(a), b, ȧ, ḃ)
            return isrenormalized(P) ? normalize_pushforward(P, Ṗ) : Ṗ
        end
    end
end

"""
    binary_frule(op, a, b, ȧ, ḃ)

Return the primal `op(a, b)` and its tangent, given the tangents `ȧ` and `ḃ` of the
arguments, for `op` being `+`, `-`, `*`, `/`, or `\\`.
"""
function binary_frule(op, a, b, ȧ, ḃ)
    Ω = op(a, b)
    iszerotangent(ȧ) && iszerotangent(ḃ) && return Ω, ZeroTangent()
    Ω̇ = binary_pushforward(op, Ω, a, b, inputtangent(a, ȧ), inputtangent(b, ḃ))
    return Ω, ProjectTo(Ω)(Ω̇)
end

########################################################################################
## Rules for the binary operations
########################################################################################

for (op, opfast) ∈ ((:+, :add_fast), (:-, :sub_fast), (:*, :mul_fast), (:/, :div_fast), (:\, nothing))
    fs = opfast === nothing ? (op,) : (op, :(Base.FastMath.$opfast))
    for f ∈ fs, (A, B) ∈ ((:AQ, :AQ), (:AQ, :Number), (:Number, :AQ))
        @eval rrule(::typeof($f), a::$A, b::$B) = binary_rrule($op, a, b)
        @eval frule((_, ȧ, ḃ)::Tuple, ::typeof($f), a::$A, b::$B) = binary_frule($op, a, b, ȧ, ḃ)
    end
end

# These are covered by the methods for `(AQ, AQ)` above, and are listed separately to make
# explicit that the normalized products of real rotors have rules.  The products of complex
# rotors are opted out in `optouts.jl`.
for (op, opfast) ∈ ((:*, :mul_fast), (:/, :div_fast)), f ∈ (op, :(Base.FastMath.$opfast))
    @eval rrule(::typeof($f), a::Rotor{<:Real}, b::Rotor{<:Real}) = binary_rrule($op, a, b)
    @eval frule((_, ȧ, ḃ)::Tuple, ::typeof($f), a::Rotor{<:Real}, b::Rotor{<:Real}) =
        binary_frule($op, a, b, ȧ, ḃ)
end

########################################################################################
## Unary operations
########################################################################################

"""
    unary_cotangent(op, Ω, q, Δ)

Return the cotangent of the argument of `Ω = op(q)`, before projection, given the
`Quaternion` cotangent `Δ` of `Ω`, for `op` being `-`, `conj`, or `inv`.  The inverse of a
`Rotor` is its conjugate, whatever its norm.
"""
unary_cotangent(::typeof(-), Ω, q, Δ) = -Δ
unary_cotangent(::typeof(conj), Ω, q, Δ) = conj(Δ)
function unary_cotangent(::typeof(inv), Ω, q, Δ)
    Ωa = qadjoint(Ω)
    return -(Ωa * Δ * Ωa)
end
unary_cotangent(::typeof(inv), Ω, q::Rotor, Δ) = conj(Δ)

"""
    unary_pushforward(op, Ω, q, q̇)

Return the tangent of `Ω = op(q)`, before projection, given the `Quaternion` tangent `q̇`
of `q`, for `op` being `-`, `conj`, or `inv`.
"""
unary_pushforward(::typeof(-), Ω, q, q̇) = -q̇
unary_pushforward(::typeof(conj), Ω, q, q̇) = conj(q̇)
function unary_pushforward(::typeof(inv), Ω, q, q̇)
    Ωq = quaternion(Ω)
    return -(Ωq * q̇ * Ωq)
end
unary_pushforward(::typeof(inv), Ω, q::Rotor, q̇) = conj(q̇)

"""
    unary_rrule(op, q)

Return the primal `op(q)` and its pullback, for `op` being `-`, `conj`, or `inv`.
"""
function unary_rrule(op, q::AQ)
    Ω = op(q)
    proj = ProjectTo(q)
    function unary_pullback(ΔΩ)
        Δ = outputcot(Ω, ΔΩ)
        Δ isa AbstractZero && return (NoTangent(), Δ)
        return NoTangent(), proj(unary_cotangent(op, Ω, q, Δ))
    end
    return Ω, unary_pullback
end

"""
    unary_frule(op, q, q̇)

Return the primal `op(q)` and its tangent, given the tangent `q̇` of `q`, for `op` being
`-`, `conj`, or `inv`.
"""
function unary_frule(op, q::AQ, q̇)
    Ω = op(q)
    iszerotangent(q̇) && return Ω, ZeroTangent()
    return Ω, ProjectTo(Ω)(unary_pushforward(op, Ω, q, inputtangent(q, q̇)))
end

# `sub_fast(q)` is `-q`, and `inv_fast(q)` is `inv(q)`.
for (f, op) ∈ ((:-, :-), (:(Base.FastMath.sub_fast), :-), (:conj, :conj), (:inv, :inv),
               (:(Base.FastMath.inv_fast), :inv))
    @eval rrule(::typeof($f), q::AQ) = unary_rrule($op, q)
    @eval frule((_, q̇)::Tuple, ::typeof($f), q::AQ) = unary_frule($op, q, q̇)
end

########################################################################################
## Products and sums of three or more arguments, and `muladd`
########################################################################################

"""
    fold_rrule(op, a, b, xs...)

Return the primal and the pullback of the left fold `op(op(op(a, b), xs[1]), ...)`, which
is how Base evaluates `*` and `+` (and the source evaluates `mul_fast` and `add_fast`) with
three or more arguments.  Each step uses `binary_rrule`.

These are rules rather than opt-outs, because ChainRules' own rule for `*` with four or more
arguments calls `rrule(*, Ω3, more...)` on the product `Ω3` of the first three arguments and
the rest; that call reaches these rules when the first quaternion is the fifth argument or
later, and an opt-out there would return `nothing`.
"""
fold_rrule(op, a, b) = binary_rrule(op, a, b)
function fold_rrule(op, a, b, c, xs...)
    ab, ab_pullback = binary_rrule(op, a, b)
    Ω, rest_pullback = fold_rrule(op, ab, c, xs...)
    function fold_pullback(ΔΩ)
        Δrest = rest_pullback(ΔΩ)
        _, ∂a, ∂b = ab_pullback(Δrest[2])
        return (NoTangent(), ∂a, ∂b, Base.tail(Base.tail(Δrest))...)
    end
    return Ω, fold_pullback
end

"""
    muladd_rrule(a, b, c)

Return the primal `muladd(a, b, c)` and its pullback.  Base promotes the three arguments to
a common type and computes `a * b + c`; promotion converts without normalizing, so the
cotangent of each promoted argument is simply projected onto the original argument.
"""
function muladd_rrule(a, b, c)
    pa, pb, pc = promote(a, b, c)
    ab, product_pullback = binary_rrule(*, pa, pb)
    Ω, sum_pullback = binary_rrule(+, ab, pc)
    function muladd_pullback(ΔΩ)
        _, ∂ab, ∂c = sum_pullback(ΔΩ)
        _, ∂a, ∂b = product_pullback(∂ab)
        return (NoTangent(), projectargs((a, b, c), (∂a, ∂b, ∂c))...)
    end
    return Ω, muladd_pullback
end

"""
    muladd_frule(a, b, c, ȧ, ḃ, ċ)

Return the primal `muladd(a, b, c)` and its tangent, given the tangents of the arguments.
"""
function muladd_frule(a, b, c, ȧ, ḃ, ċ)
    pa, pb, pc = promote(a, b, c)
    ab, ab_dot = binary_frule(*, pa, pb, ȧ, ḃ)
    return binary_frule(+, ab, pc, ab_dot, ċ)
end

for (A, B, C) ∈ Iterators.product(ntuple(i -> (:AQ, :Number), 3)...)
    :AQ ∈ (A, B, C) || continue
    for (f, op) ∈ ((:*, :*), (:+, :+), (:(Base.FastMath.mul_fast), :*), (:(Base.FastMath.add_fast), :+))
        @eval rrule(::typeof($f), a::$A, b::$B, c::$C, xs::Vararg{Number}) =
            fold_rrule($op, a, b, c, xs...)
    end
    @eval rrule(::typeof(muladd), a::$A, b::$B, c::$C) = muladd_rrule(a, b, c)
    @eval frule((_, ȧ, ḃ, ċ)::Tuple, ::typeof(muladd), a::$A, b::$B, c::$C) =
        muladd_frule(a, b, c, ȧ, ḃ, ċ)
end
