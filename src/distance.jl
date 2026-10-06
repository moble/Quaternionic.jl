"""
    distance(q₁, q₂)
    distance2(q₁, q₂)

Measure the "distance" between two quaternions, or the squared distance with `distance2`.

By default, these functions use the natural measure in the *additive* group of quaternions,
so that `distance2` returns
```julia
abs2(q₁ - q₂)
```
and `distance` returns its square root.  For quaternions with complex components, `abs2` is
the spinor norm, so these values may be complex, and they may be zero for distinct
quaternions.

If both arguments are `Rotor`s, the functions use the natural measure in the *rotation*
group instead, so that `distance2` returns
```julia
min(abs2(log(q₁ / q₂)), abs2(log(-q₁ / q₂)))
```
and `distance` returns its square root.  This is half the angle (in radians) of the rotation
that takes `q₂` to `q₁`, and it lies in the range ``[0, π/2]``.  For example, `rotor(imx)`
and `rotor(imy)` represent rotations through ``π`` about the ``x`` and ``y`` axes.  The
rotation that takes one to the other is a rotation through ``π`` about the ``z`` axis, so
their distance is ``π/2``.  These methods are (efficiently) independent of the scaling of
`q₁` and `q₂`, including factors of -1, as is appropriate for the rotation group.  The
`Rotor` methods do not support complex components (`Lorentz` rotors), because the Lorentz
group has no positive-definite bi-invariant metric; they throw an `ArgumentError` for such
input.

# Examples
```jldoctest example
julia> distance(imx, imy)
1.4142135623730951
julia> distance(rotor(imx), rotor(imy))
1.5707963267948966
julia> distance(imz, -imz)
2.0
julia> distance(rotor(imz), rotor(-imz))
0.0
```
"""
distance(q₁::AbstractQuaternion, q₂::AbstractQuaternion) = √distance2(q₁, q₂)
distance2(q₁::AbstractQuaternion, q₂::AbstractQuaternion) = abs2(q₁ - q₂)
distance(q₁::Rotor, q₂::Rotor) = √distance2(q₁, q₂)
distance2(q₁::Rotor, q₂::Rotor) = _abs2_small_vec_log(q₁ / q₂)

# Like `min(abs2(log(q)), abs2(log(-q)))`, but assumes that the norm of `q` is 1.  The
# threshold that selects the series is chosen by the type of the value of the scalar part,
# so that dual numbers take the same branches as the floats they wrap.
_abs2_small_vec_log(q::Rotor) = _abs2_small_vec_log(value(q[1]), q)

@inline function _abs2_small_vec_log(w, q::Rotor)
    v² = abs2vec(q)
    w² = q[1]^2
    # This compares `x = v²/w²` with the threshold, but only on the values, and without
    # dividing, so that nothing singular is evaluated when `q[1]` is zero.  For element
    # types whose values are not floats, such as `ReverseDiff.TrackedReal`, the threshold
    # for `Float64` is used.  For symbolic types, the comparison does not evaluate to a
    # `Bool`, so the closed form is used.
    T = w isa AbstractFloat ? typeof(w) : Float64
    small = value(v²) ≤ (sqrt(sqrt(eps(T))) / 2) * value(w²)
    if small === true
        # Near the identity, `absvec` is the square root of a tiny number, whose derivatives
        # are huge and lose all accuracy when they cancel, or are NaN when it is exactly 0.
        # This series for `atan(√x)^2` in `x` is smooth instead; its relative truncation
        # error is about `x⁴/3`, which is below `eps` for this range of `x`.
        x = v² / w²
        x * evalpoly(x, (1, -2//3, 23//45, -44//105))
    else
        atan(absvec(q), abs(q[1]))^2
    end
end

function _abs2_small_vec_log(::Complex, q::Rotor)
    throw(ArgumentError(
        "`distance` and `distance2` are not defined for `Rotor`s with complex components "
        * "(`Lorentz` rotors), because the Lorentz group has no positive-definite "
        * "bi-invariant metric.  Use `distance(quaternion(q₁), quaternion(q₂))` for the "
        * "additive measure, which is complex in general."
    ))
end
