# Essential elements of making quaternions into an algebra

"""
    conj(q)

Return the quaternion conjugate, which flips the sign of each "vector"
component.

# Examples
```jldoctest
julia> conj(quaternion(1,2,3,4))
1 - 2𝐢 - 3𝐣 - 4𝐤
```
"""
Base.conj(q::Q) where {Q<:AbstractQuaternion} = Q(q[1], -q[2], -q[3], -q[4])
Base.conj(q::Q) where {Q<:AbstractQuaternion{Bool}} = wrapper(Q)(q[1], -q[2], -q[3], -q[4])

Base.:-(q::Q) where {Q<:AbstractQuaternion} = Q(-components(q))
Base.:-(q::Q) where {Q<:AbstractQuaternion{Bool}} = wrapper(Q)(-Int.(components(q)))

# Note that, in the two definitions above, the `wrapper` function is used to
# ensure that the result is of the same type as the input quaternion, but with a different
# basetype.  This is necessary because `-true` is not a `Bool`, it's the `Int` -1.  This is
# the only type I can think of where `-` changes the type of the input, so that's the only
# special case we need.  (I'm assuming `Unsigned` types are not used here.)



for TA ∈ (AbstractQuaternion, Rotor, QuatVec)
    for TB ∈ (AbstractQuaternion, Rotor, QuatVec)
        @eval begin
            Base.:+(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(+), $TB)(components(q)+components(p))
            Base.:-(q::T1, p::T2) where {T1<:$TA, T2<:$TB} = wrapper($TA, Val(-), $TB)(components(q)-components(p))
        end
    end
    let TB = Number
        @eval begin
            Base.:+(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(+), $TB)(q[1]+p, q[2], q[3], q[4])
            Base.:-(q::QT, p::$TB) where {QT<:$TA} = wrapper($TA, Val(-), $TB)(q[1]-p, q[2], q[3], q[4])
            Base.:+(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(+), $TA)(p+q[1], q[2], q[3], q[4])
            Base.:-(p::$TB, q::QT) where {QT<:$TA} = wrapper($TB, Val(-), $TA)(p-q[1], -q[2], -q[3], -q[4])
        end
    end
end


function Base.:*(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    wrapper(Q1, Val(*), Q2)(
        q[1]*p[1] - q[2]*p[2] - q[3]*p[3] - q[4]*p[4],
        q[1]*p[2] + q[2]*p[1] + q[3]*p[4] - q[4]*p[3],
        q[1]*p[3] - q[2]*p[4] + q[3]*p[1] + q[4]*p[2],
        q[1]*p[4] + q[2]*p[3] - q[3]*p[2] + q[4]*p[1]
    )
end


# There is deliberately no shortcut returning `one` when `p == q`.  Such a branch would
# return a constant, so AD backends whose `==` compares only primal values (Enzyme,
# Mooncake, ReverseDiff, and ForwardDiff 0.10) would lose the derivative there.  It would
# also make the return type depend on the values, and it would fail for symbolic types
# whose `==` does not return a `Bool`.  Instead, each vector component is written as a sum
# of two differences, `(a*b - c*d) + (e*f - g*h)`, in which the terms of each difference
# are identical products when `p == q`.  Both differences then vanish exactly in
# floating-point arithmetic.  The scalar part sums the products in the same order as
# `abs2`, so `q / q` is exactly one for real floating-point components.  (For complex
# components, the final complex division may still round, by about one ulp.)
#
# For real floating-point components, `abs2(p)` overflows when `abs(p)` exceeds about
# `sqrt(floatmax)`, and it underflows to a subnormal number or zero when `abs(p)` falls below
# about `sqrt(floatmin)`.  In that case both `q` and `p` are multiplied by a power of two that
# brings the largest component of `p` near 1.  This scaling is exact, so it changes neither
# the quotient nor the exactness of `q / q`.  The power is clamped to the range of normal
# numbers, so that it stays finite when the components of `p` are subnormal; the recursive
# call then scales again if one step was not enough.  The test is applied to `value(den)`,
# so that dual numbers are rescaled in the same way as the floats they contain.  It begins
# with type checks that the compiler resolves, so it is skipped entirely for `Rotor`s, whose
# `abs2` is identically one, and for other element types, such as complex and symbolic
# types.  The `maximum` is evaluated, using values only, when the scaling is needed.
function Base.:/(q::Q1, p::Q2) where {Q1<:AbstractQuaternion, Q2<:AbstractQuaternion}
    den = abs2(p)
    if p isa Union{Quaternion, QuatVec} && value(den) isa AbstractFloat
        d = value(den)
        T = typeof(d)
        if !(floatmin(T) ≤ d < T(Inf))
            m = maximum(x -> abs(value(x)), components(p))
            if !iszero(m) && isfinite(m)
                k = clamp(-exponent(m), exponent(floatmin(T)), exponent(floatmax(T)))
                σ = ldexp(one(T), k)
                return wrapper(Q1, Val(/), Q2)(components((σ * q) / (σ * p)))
            end
        end
    end
    wrapper(Q1, Val(/), Q2)(
        (q[1]*p[1] + q[2]*p[2] + q[3]*p[3] + q[4]*p[4]) / den,
        ((q[2]*p[1] - q[1]*p[2]) + (q[4]*p[3] - q[3]*p[4])) / den,
        ((q[3]*p[1] - q[1]*p[3]) + (q[2]*p[4] - q[4]*p[2])) / den,
        ((q[4]*p[1] - q[1]*p[4]) + (q[3]*p[2] - q[2]*p[3])) / den
    )
end


let S = Number
    @eval begin
        Base.:*(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(*), $S)(s*components(p))
        Base.:*(s::$S, p::Q) where {Q<:AbstractQuaternion} = wrapper($S, Val(*), Q)(s*components(p))
        Base.:/(p::Q, s::$S) where {Q<:AbstractQuaternion} = wrapper(Q, Val(/), $S)(components(p)/s)
        function Base.:/(s::$S, p::Q) where {Q<:AbstractQuaternion}
            den = abs2(p)
            # As in quaternion division above, `s` and `p` are rescaled by a power of two
            # when `abs2(p)` overflows or underflows.
            if p isa Union{Quaternion, QuatVec} && value(den) isa AbstractFloat
                d = value(den)
                T = typeof(d)
                if !(floatmin(T) ≤ d < T(Inf))
                    m = maximum(x -> abs(value(x)), components(p))
                    if !iszero(m) && isfinite(m)
                        k = clamp(-exponent(m), exponent(floatmin(T)), exponent(floatmax(T)))
                        σ = ldexp(one(T), k)
                        return wrapper($S, Val(/), Q)(components((σ * s) / (σ * p)))
                    end
                end
            end
            f = s / den
            wrapper($S, Val(/), Q)(p[1] * f, -p[2] * f, -p[3] * f, -p[4] * f)
        end
    end
end


# `@fastmath` replaces each arithmetic operator with the corresponding function from
# `Base.FastMath`, whose fallbacks `promote` mixed arguments before applying the operator.
# That would lose the careful choice of return types made by `wrapper` above: for example,
# `@fastmath 1.0 * imz` would return a `Quaternion` rather than a `QuatVec`, and `@fastmath
# 2.0 * R` would convert `2.0` to a `Rotor` and return a `Rotor` equal to `R`.  (See issue
# #46.)  So we skip the promotion and apply the ordinary operators directly.
for (op, op_fast) ∈ ((:+, :add_fast), (:-, :sub_fast), (:*, :mul_fast), (:/, :div_fast))
    @eval begin
        Base.FastMath.$op_fast(q::AbstractQuaternion, p::AbstractQuaternion) = $op(q, p)
        Base.FastMath.$op_fast(q::AbstractQuaternion, s::Number) = $op(q, s)
        Base.FastMath.$op_fast(s::Number, q::AbstractQuaternion) = $op(s, q)
    end
end

# Expressions like `@fastmath a * b * c` become a single call to `mul_fast(a, b, c)`, whose
# fallback promotes all of its arguments at once.  We catch any such call with a quaternion
# among its first three arguments, and reduce it to calls with two arguments.  A call with
# plain numbers in all of the first three positions cannot be caught without type piracy, so
# it will still promote; for example, `@fastmath 1.0 * 2.0 * 3.0 * imz` is a `Quaternion`.
for op_fast ∈ (:add_fast, :mul_fast)
    for (A, B, C) ∈ Iterators.product(ntuple(_ -> (AbstractQuaternion, Number), 3)...)
        AbstractQuaternion ∈ (A, B, C) || continue
        @eval Base.FastMath.$op_fast(a::$A, b::$B, c::$C, xs::Number...) =
            Base.FastMath.$op_fast(Base.FastMath.$op_fast(a, b), c, xs...)
    end
end


@doc raw"""
    p ⋅ q

Evaluate the inner ("dot") product between two quaternions.  Equal to the
scalar part of `p * conj(q)`.

Note that this function is not very commonly used, except as a quick way to
determine whether the two quaternions are more anti-parallel than parallel, for
functions like [`unflip`](@ref).

This method extends `LinearAlgebra.dot`, of which `⋅` is an alias, and returns only the
componentwise sum `p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4]`, which is the Euclidean
inner product for real components.  Generic linear-algebra code written for
`Complex` numbers usually assumes instead that the `dot` of two scalars is the full product
`conj(p) * q`.  Consequently, functions that rely on that assumption, such as `x' * y` for
vectors, `qr`, `svd`, and least-squares solutions with `\`, do not give the quaternionic
results for arrays of quaternions.

This can be typed as `\cdot<tab>` in a Julia-aware editor.
"""
@inline function LinearAlgebra.:⋅(p::AbstractQuaternion, q::AbstractQuaternion)
    p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4]
end

@doc raw"""
    a × b
    cross(a, b)

Return the cross product of two pure-vector quaternions.  Equal to ½ of the
commutator product `a*b-b*a`.

This extends `LinearAlgebra.cross`, of which `×` is an alias, so the same function is
used whether `×` comes from `Quaternionic` or from `LinearAlgebra`.

This can be typed as `\times<tab>` in a Julia-aware editor.
"""
@inline function LinearAlgebra.:×(a::QuatVec, b::QuatVec)
    quatvec(
        a[3] * b[4] - a[4] * b[3],
        a[4] * b[2] - a[2] * b[4],
        a[2] * b[3] - a[3] * b[2]
    )
end

@doc raw"""
    a ×̂ b

Return the *direction* of the cross product between `a` and `b`; the normalized vector along
[`a×b`](@ref LinearAlgebra.:×) — unless the magnitude is zero, in which case the zero vector
is returned.  Both cases give the same element type, which is a floating-point type for
integer or rational inputs.

For complex components, the magnitude is the spinor norm computed by [`absvec`](@ref), which
can vanish for a nonzero null vector.  Such a vector cannot be normalized, so it is returned
unchanged.

This can be typed as `\times<tab>\hat<tab>` in a Julia-aware editor.
"""
@inline function ×̂(a::QuatVec, b::QuatVec)
    axb = a × b
    av = absvec(axb)
    if iszerovalue(av)
        # Dividing by one gives this branch the same element type as the other branch,
        # which matters when `a` and `b` have integer or rational components.
        return axb / one(av)
    else
        return axb / av
    end
end

"""
    normalize(q)

Return a copy of this quaternion, normalized.

Note that this returns the same type as the input quaternion.  If you want to
convert to a `Rotor`, just call `rotor(q)`, which includes a normalization
step.

This extends `LinearAlgebra.normalize` for quaternion types.
"""
@inline LinearAlgebra.normalize(q::AbstractQuaternion) = q / abs(q)
@inline LinearAlgebra.normalize(q::Rotor) = rotor(q)

# The function-call syntax `R(v)` computes the sandwich product `R * v * conj(R)`.  Both
# methods below use the same homogeneous formula for the vector part, which is valid for
# any `R`, whatever its norm.  The product equals the rotation `R * v * inv(R)` only when
# `abs(R) == 1`; otherwise it is that rotation scaled by `abs2(R)`.  The `Rotor` method
# accepts only a `QuatVec`, and returns a `QuatVec`.  The `Quaternion` method sets the
# scalar part to `abs2(R) * v[1]`, and wraps the result in the wrapper type of `v`, which
# renormalizes the result when `v` is a `Rotor`.
function (R::Rotor)(v::QuatVec)
    quatvec(SA[
        false,
        ((R[1]^2 + R[2]^2 - R[3]^2 - R[4]^2)*v[2]
            + (R[1]*R[3] + R[2]*R[4])*2v[4] + (R[2]*R[3] - R[1]*R[4])*2v[3]),
        ((R[1]^2 - R[2]^2 + R[3]^2 - R[4]^2)*v[3]
            + (R[2]*R[3] + R[1]*R[4])*2v[2] + (R[3]*R[4] - R[1]*R[2])*2v[4]),
        ((R[1]^2 + R[4]^2 - R[2]^2 - R[3]^2)*v[4]
            + (R[1]*R[2] + R[3]*R[4])*2v[3] + (R[2]*R[4] - R[1]*R[3])*2v[2])
    ])
end
function (R::Quaternion)(v::QT) where {QT<:AbstractQuaternion}
    wrapper(QT)(
        abs2(R) * v[1],
        ((R[1]^2 + R[2]^2 - R[3]^2 - R[4]^2)*v[2]
            + (R[1]*R[3] + R[2]*R[4])*2v[4] + (R[2]*R[3] - R[1]*R[4])*2v[3]),
        ((R[1]^2 - R[2]^2 + R[3]^2 - R[4]^2)*v[3]
            + (R[2]*R[3] + R[1]*R[4])*2v[2] + (R[3]*R[4] - R[1]*R[2])*2v[4]),
        ((R[1]^2 + R[4]^2 - R[2]^2 - R[3]^2)*v[4]
            + (R[1]*R[2] + R[3]*R[4])*2v[3] + (R[2]*R[4] - R[1]*R[3])*2v[2])
    )
end
