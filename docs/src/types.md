# Types and construction

From `AbstractQuaternion{T}` we define three subtypes:

  * `Quaternion{T}`, which is an element of the general algebra of
    quaternions over any `T<:Number`.  Most commonly, `T` is a real
    floating-point type, but integers, symbolic types, and complex
    numbers are also supported.  Complex components are used to
    represent [Lorentz transformations](@ref).
  * `Rotor{T}`, which is an element of the multiplicative group of
    unit quaternions, and is interpreted as mapping to a rotation.
    The magnitude is *assumed* to be 1 (though, for efficiency, this
    is not generally confirmed), and the sign may be freely changed in
    certain cases.
  * `QuatVec{T}`, which is an element of the additive group of
    quaternions with 0 scalar part; a "pure vector" quaternion.

For simplicity, almost every function in this package is defined for
general `Quaternion`s, so you may not need any other type.  However,
it can frequently be more accurate *and* more efficient to use the
other subtypes where relevant.

## Constructors, constants, and conversions

At the most basic level, `Quaternion{T}` mimics `Complex{T}` as
closely as possible, including the behavior of most functions in
`Base`.  The `Rotor{T}` and `QuatVec{T}` subtypes behave very
similarly, except that most of their constructors automatically impose
the constraints that the norm is 1 and the scalar component is 0,
respectively.  Also note that when a certain operation is not defined
for either of those subtypes, the functions will usually convert to a
general `Quaternion` automatically.

To create new `Quaternion`s interactively, it is typically most
convenient to use the constants `imx`, `imy`, and `imz` — or
equivalently `𝐢`, `𝐣`, and `𝐤` — multiplied by appropriate factors
and added together.  For programmatic work, it is more common to use
the [`quaternion`](@ref) function — which takes all four components,
the three vector components, or just the one scalar component, and
creates a new `Quaternion` of the type implied by the arguments.  The
[`rotor`](@ref) and [`quatvec`](@ref) functions do the same for the
other subtypes.  You can also *specify* the type, as in
`Quaternion{Float64}(...)`.  Type conversions with `promote`, `widen`,
`float`, etc., work as expected.

```@docs
AbstractQuaternion
Quaternion
Rotor
QuatVec
quaternion
rotor
quatvec
imx
imy
imz
𝐢
𝐣
𝐤
components
basetype
```

## Number functions from Base

The standard [`Number`
functions](https://docs.julialang.org/en/v1/base/numbers/#General-Number-Functions-and-Constants)
that work for `Complex`, such as `isfinite`, `iszero`, etc., should
work analogously for `Quaternion`.  The `hash`, `read`, and `write`
functions are also implemented.  Because quaternions are `Number`s,
broadcasting treats each quaternion as a single scalar: `abs.(q)` is
the same as `abs(q)`, for example.  To apply a function to each
component, apply it to [`components`](@ref)`(q)` and construct a new
quaternion from the result.

## Random quaternions

It is frequently convenient to construct random `Quaternion` objects,
which can be done just as with other types by passing the desired
output type to the [`randn`](@ref) function.  The `rand` function is
not overloaded, because there would be no geometric significance to
such a `Quaternion`; `randn` results are independent of the
orientation of the basis used to define the quaternions.  Note that it
is possible to get random *rotors* and *vectors* by passing the
appropriate types to the `randn` function.

```@autodocs
Modules = [Quaternionic]
Pages   = ["random.jl"]
```

## Symbolic components and LaTeX display

Quaternions with `Num` components from
[`Symbolics.jl`](https://docs.sciml.ai/Symbolics/stable/) support the
basic algebra and most of the mathematical functions.  A few
functions, such as `sqrt`, choose an algorithm by comparing the values
of components, and raise an error ("non-boolean (Num) used in boolean
context") for symbolic components.  To simplify each component of a
symbolic quaternion, use `quaternion(simplify.(components(q))...)`.

When [`Latexify.jl`](https://github.com/korsbo/Latexify.jl) is loaded,
`latexify(q)` returns the LaTeX form of a quaternion, and quaternions
are displayed as LaTeX in environments that support the `text/latex`
MIME type, such as Jupyter notebooks.
