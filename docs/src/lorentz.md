# Lorentz transformations

Quaternions with complex components represent proper orthochronous
Lorentz transformations — more precisely, elements of their double
cover ``\mathrm{Spin}^+(3,1)`` — through the isomorphism between the
complex quaternions and the even subalgebra of the spacetime algebra.
The theory, including the sign conventions, is described on the
[Quaternions and the Spacetime Algebra](@ref) page; this page
describes how to use it.

The [`Lorentz`](@ref) type is simply an alias for `Rotor{Complex{T}}`.
The [`Boost`](@ref) function constructs a pure boost from a rapidity
and a direction, or from a velocity (in units where ``c = 1``), and
`Lorentz(R)` converts a real `Rotor` `R` into the corresponding pure
rotation.  Lorentz transformations are composed with `*` and inverted
with `inv`, and calling one on a Minkowski 4-vector ``[t, x, y, z]``
(with signature ``{-}{+}{+}{+}``) transforms it:
```jldoctest lorentz
julia> using Quaternionic

julia> B = Boost([0.6, 0.0, 0.0]);  # A boost with velocity 0.6 along x

julia> B([1, 0, 0, 0]) ≈ [1.25, 0.75, 0, 0]
true

julia> Boost(atanh(0.6), [1, 0, 0]) ≈ B  # The same boost, given by its rapidity
true
```

The action is *active*: `Λ(v)` is the transformed vector, expressed in
the original frame.  Above, `B` maps the 4-velocity of a particle at
rest to the 4-velocity ``\gamma\,(1, \vec{v})`` of a particle moving
with velocity ``\vec{v}``.  The components of a fixed vector as seen
by an observer moving with velocity ``\vec{v}`` are given instead by
the passive transformation `inv(Boost(v⃗))`.  Composition acts from
right to left:
```jldoctest lorentz
julia> R = Lorentz(exp(π/8 * imz));  # A rotation by π/4 about z

julia> v = [1.0, 2.0, 3.0, 4.0];

julia> (R * B)(v) ≈ R(B(v))
true

julia> inv(B)(B(v)) ≈ v
true
```

Every Lorentz transformation can be factored into simpler pieces.  The
functions [`RB`](@ref Quaternionic.RB) and [`BR`](@ref
Quaternionic.BR) factor it into a rotation and a boost (in either
order), [`Rv`](@ref Quaternionic.Rv) and [`vR`](@ref Quaternionic.vR)
do the same but return the boost's velocity, and [`KAN`](@ref
Quaternionic.KAN) computes the Iwasawa decomposition described in
[Iwasawa's ``KA\,N`` decomposition](@ref iwasawa-kan):
```jldoctest lorentz
julia> K, A, N = Quaternionic.KAN(R * B);

julia> K * A * N ≈ R * B
true
```

The functions [`ℂreal`](@ref Quaternionic.ℂreal), [`ℂimag`](@ref
Quaternionic.ℂimag), [`ℂreim`](@ref Quaternionic.ℂreim), and
[`ℂconj`](@ref Quaternionic.ℂconj) operate on the real and imaginary
parts of the complex components.  These functions are public but not
exported, so they must be qualified with the module name, or imported
explicitly.

Many generic quaternion functions, such as `exp`, `log`, `sqrt`, and
integer powers, also accept complex quaternions, and `abs` and `abs2`
compute the complex spinor norm.  Functions that are specific to
rotations in three dimensions, such as `distance`, `to_euler_angles`,
and `from_rotation_matrix`, do not support `Lorentz` input.

## Reference

The entry `Lorentz(::AbstractVector)` below documents the action
`Λ(v)` of a Lorentz transformation `Λ` on a 4-vector `v`.

```@meta
CurrentModule = Quaternionic
```

```@docs
Lorentz
Lorentz(::AbstractVector)
Boost
ga_components
RB
BR
Rv
vR
KAN
ℂreal
ℂimag
ℂreim
ℂconj
```

```@meta
CurrentModule = nothing
```
