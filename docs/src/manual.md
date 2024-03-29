# Quaternionic functions

From `AbstractQuaternion{T}` we define three subtypes:

  * `Quaternion{T}`, which is an element of the general algebra of quaternions
    over any `T<:Real`.
  * `Rotor{T}`, which is an element of the multiplicative group of unit
    quaternions, and is interpreted as mapping to a rotation.  The magnitude is
    *assumed* to be 1 (though, for efficiency, this is not generally
    confirmed), and the sign may be freely changed in certain cases.
  * `QuatVec{T}`, which is an element of the additive group of quaternions with
    0 scalar part; a "pure vector" quaternion.

For simplicity, almost every function in this package is defined for general
`Quaternion`s, so you may not need any other type.  However, it can frequently
be more accurate *and* more efficient to use the other subtypes where relevant.

## Constructors, constants, and conversions

At the most basic level, `Quaternion{T}` mimics `Complex{T}` as closely as
possible, including the behavior of most functions in `Base`.  The `Rotor{T}`
and `QuatVec{T}` subtypes behave very similarly, except that most of their
constructors automatically impose the constraints that the norm is 1 and the
scalar component is 0, respectively.  Also note that when a certain operation
is not defined for either of those subtypes, the functions will usually convert
to a general `Quaternion` automatically.

To create new `Quaternion`s interactively, it is typically most convenient to
use the constants `imx`, `imy`, and `imz` — or equivalently `𝐢`, `𝐣`, and `𝐤` —
multiplied by appropriate factors and added together.  For programmatic work,
it is more common to use the [`Quaternion`](@ref) function — which takes all
four components, the three vector components, or just the one scalar component,
and creates a new `Quaternion` of the type implied by the arguments.  You can
also *specify* the type, as in `Quaternion{Float64}(...)`.  Type conversions
with `promote`, `widen`, `float`, etc., work as expected.

```@autodocs
Modules = [Quaternionic]
Pages   = ["quaternion.jl"]
Order   = [:module, :type, :constant, :function, :macro]
```

## Number functions from Base

The standard [`Number`
functions](https://docs.julialang.org/en/v1/base/numbers/#General-Number-Functions-and-Constants)
that work for `Complex`, such as `isfinite`, `iszero`, etc., should work
analogously for `Quaternion`.  The `hash`, `read`, and `write` functions are
also implemented.  As noted in the [Examples](@ref), broadcasting to each
component is also implemented via `broadcasted`.

```@autodocs
Modules = [Quaternionic]
Pages   = ["base.jl"]
```


## Algebra and mathematical functions

Along with the basic binary operators, the essential mathematical functions
like [`conj`](@ref), [`abs`](@ref), [`abs2`](@ref), [`exp`](@ref),
[`log`](@ref), etc., are implemented.  Most of these functions are found in the
`Base` module, and are simply overloaded methods of functions that should also
be familiar from `Complex` types.  Note that we use a slightly different
interpretation of [`angle`](@ref) for `Quaternion`, compared to `Complex`.  We
also have [`absvec`](@ref) and [`abs2vec`](@ref), which are not useful in a
`Complex` context, but compute the relevant quantities for the "vector"
component of a `Quaternion`.

```@autodocs
Modules = [Quaternionic]
Pages   = ["algebra.jl", "math.jl"]
```


## Random quaternions

It is frequently convenient to construct random `Quaternion` objects, which can
be done just as with other types by passing the desired output type to the
[`randn`](@ref) function.  The `rand` function is not overloaded, because there
would be no geometric significance to such a `Quaternion`; `randn` results are
independent of the orientation of the basis used to define the quaternions.
Note that it is possible to get random *rotors* and *vectors* by passing the
appropriate types to the `randn` function.

```@autodocs
Modules = [Quaternionic]
Pages   = ["random.jl"]
```


## Conversions

It can sometimes be useful to convert between quaternions and other
representations.  Most of these functions are named `to_<representation>` and
have a corresponding `from_<representation>` function.  Furthermore, most
convert to/from representations of rotations.  While rotations are not the only
useful application of quaternions, they are probably the most common.  The only
conversions that are not specifically related to rotations are
[`to_float_array`](@ref) and [`from_float_array`](@ref).

```@autodocs
Modules = [Quaternionic]
Pages   = ["conversion.jl"]
```


## Distances

There are several ways of measuring the "distance" between two quaternions:
``d(q_1, q_2)``.  Fundamentally, any comparison between two quaternions ``q_1``
and ``q_2`` must make use of a binary operation, for which there are two
obvious choices: addition or multiplication.  For either choice, we operate on
``q_1`` and the appropriate inverse (either additive or multiplicative) of
``q_2``.  That is, ``d`` should be a function of either ``q_1 - q_2`` or
``q_1/q_2``.[^1]

[^1]:
    For ``q_1/q_2``, we are dealing with the *multiplicative* group of
    quaternions, which does not include 0, so we will assume that no quaternion
    involved in such a function can be 0.

Now, we also have a number of criteria we would like any distance function to
satisfy.  For any quaternions ``q_1`` and ``q_2`` and any *unit* quaternion
``q_3``, we require

  * real-valued: ``d(q_1, q_2) \in \mathbb{R}``
  * symmetry: ``d(q_1, q_2) = d(q_2, q_1)``
  * invariance: ``d(q_3\, q_1, q_3\, q_2) = d(q_1, q_2) = d(q_1\, q_3, q_2\, q_3)``
  * identity: ``d(q_1, q_1) = 0``
  * positive-definiteness: ``d(q_1, q_2) > 0`` whenever ``q_1 ≠ q_2``

(Of course, it should be noted that these criteria all hold in the *exact* case;
when using floating-point numbers, they will likely be violated near edge
cases.)

It is not hard to see that these criteria can be satisfied by any of

  * `abs(q₁ - q₂)`
  * `abs2(q₁ - q₂)`
  * `abs(log(q₁ / q₂))`
  * `abs2(log(q₁ / q₂)`

If ``q_1`` and ``q_2`` are interpreted as rotations, we frequently don't care
about their signs, and just want the *smallest* distance between them, for any
choice of sign.  Furthermore, in the multiplicative case, the `log` functions
will involve calculation of the `log` of the magnitudes of the quaternions,
which should be 1.  In this case, we relax the "positive-definiteness" criterion
to allow ``d(q_1, q_2)`` to equal zero when ``q_1`` and ``q_2`` are related by a
nonzero scalar multiple.

For `Rotor` types, the latter two multiplicative options are most relevant,
while for other types the additive options are more relevant.  These are the
default behaviors of the `distance` and `distance2` functions.

```@autodocs
Modules = [Quaternionic]
Pages   = ["distance.jl"]
```


## Alignment

There are many ways to optimize alignment with rotations.  In particular, we can
seek the optimal rotation that takes one set of points onto a corresponding set
of points, or the optimal quaternion that takes one set of quaternions onto a
corresponding set of quaternions.  In both cases, the "optimal" value depends on
the metric being used.  The simplest and most robust results are obtained when
the metric is the standard Euclidean metric (in the case of points), or the
magnitude of the difference (in the case of quaternions).  Here, we assume that
`QuatVec`s represent points, and any other type of quaternion should be treated
as rotors.

```@autodocs
Modules = [Quaternionic]
Pages   = ["alignment.jl"]
```
