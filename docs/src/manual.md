# Quaternionic functions

## Constructors and Number methods

At the most basic level, `Quaternion{T}` mimics `Complex{T}` as closely as possible, including the
behavior of most functions in `Base`.

To create new `Quaternion`s interactively, it is typically most convenient to use the constants
`imx`, `imy`, and `imz` — or equivalently `𝐢`, `𝐣`, and `𝐤` — multiplied by appropriate factors and
added together.  For programmatic work, it is more common to use the [`Quaternion`](@ref) function —
which takes all four components, the three vector components, or just the one scalar component, and
creates a new `Quaternion` of the type implied by the arguments.  You can also *specify* the type,
as in `Quaternion{Float64}(...)`.  Type conversions with `promote`, `widen`, `float`, etc., work as
expected.  The standard [`Number`
functions](https://docs.julialang.org/en/v1/base/numbers/#General-Number-Functions-and-Constants)
that work for `Complex`, such as `isfinite`, `iszero`, etc., should work analogously for
`Quaternion`.  The `hash`, `read`, and `write` functions are also implemented.  As noted in the
[Examples](@ref), broadcasting to each component is also implemented via `broadcasted`.

```@autodocs
Modules = [Quaternionic]
Pages   = ["quaternion.jl"]
```


## Algebra and mathematical functions

The essential mathematical features of quaternions are implemented as functions like [`conj`](@ref),
[`abs`](@ref), [`abs2`](@ref), [`exp`](@ref), [`log`](@ref), etc.  Most of these functions are
elements of the `Base` module, and are simply overloaded methods of functions that should also be
familiar from `Complex` types.  We also have [`absvec`](@ref) and [`abs2vec`](@ref), which are not
useful in a `Complex` context, but compute the relevant quantities for the "vector" component of a
`Quaternion`.

```@autodocs
Modules = [Quaternionic]
Pages   = ["algebra.jl", "math.jl"]
```


## Random quaternions

It is frequently convenient to construct random `Quaternion` objects, which can be done just as with
other types by passing the desired output type to the [`randn`](@ref) function.  The `rand` function
is not overloaded, because there would be no geometric significance to such a `Quaternion`; `randn`
results are independent of the orientation of the basis used to define the quaternions.  A simple
convenience function [`randn_rotor`](@ref) is also provided, to normalize each result.

```@autodocs
Modules = [Quaternionic]
Pages   = ["random.jl"]
```


## Conversions

It can sometimes be useful to convert between quaternions and other representations.  Most of these
functions are named `to_<representation>` and have a corresponding `from_<representation>` function.
Furthermore, most convert to/from representations of rotations.  While rotations are not the only
useful application of quaternions, they are probably the most common.  The only conversions that are
not specifically related to rotations are [`to_float_array`](@ref) and [`from_float_array`](@ref).

```@autodocs
Modules = [Quaternionic]
Pages   = ["conversion.jl"]
```
