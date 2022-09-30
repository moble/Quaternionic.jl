# Introduction

*Quaternions for Julia*

The goal of this package is to provide a simple but flexible and complete implementation of
quaternions, without restricting the interpretation of quaternions to being rotations — while also
providing extensive support for rotations — along with thorough testing, documentation, and
integration with the rest of Julia.  Wherever possible, standard functions that work with `Complex`
will also work with `Quaternion`.

In addition to a basic `Quaternion{T}` type, we also have [`Rotor{T}`](@ref) and
[`QuatVec{T}`](@ref) specializations, which can improve the accuracy and efficiency of certain
applications.  Each of these can be defined over any `T<:Real`; in addition to the standard
primitive types (`Float64`, etc.), `BigFloat` and `Symbolics.Num` are tested extensively.

## Examples

There are numerous ways to construct a [`Quaternion`](@ref) — the simplest being to just give the
components:
```jldoctest example
julia> using Quaternionic

julia> q = Quaternion(1.0, 2.0, 3.0, 4.0)
1.0 + 2.0𝐢 + 3.0𝐣 + 4.0𝐤
julia> p = Quaternion(4, 3, 2, 1)
4 + 3𝐢 + 2𝐣 + 1𝐤
```
Each quaternion type is parametrized by the types of its components (which are promoted to be all
the same type).  Any subtype of `Real` is allowed, and is detected automatically.  For example,
`q` has type `Quaternion{Float64}`, while `p` has type `Quaternion{Int64}`.[^1] The base type may
be given explicitly if desired, to override the detected type:
```jldoctest example
julia> r = Quaternion{Float64}(4, 3, 2, 1)
4.0 + 3.0𝐢 + 2.0𝐣 + 1.0𝐤
```
The various `Float` and `Int` types work well, as do `BigFloat`, and the [`Num` type from
`Symbolics.jl`](https://symbolics.juliasymbolics.org/v0.1/manual/variables/#A-note-about-functions-restricted-to-Numbers-1).
In particular, we can use symbolic expressions as components:
```jldoctest symbolics
julia> using Quaternionic, Symbolics

julia> @variables a b c d e;

julia> Quaternion(a-b, b*c, c/d, d+e)
a - b + b*c𝐢 + (c / d)𝐣 + (d + e)𝐤
```
It is also possible to construct random quaternions using [`randn`](@ref) with a `Quaternion` type.
In analogy with the complex types, the aliases `QuaternionF64`, `QuaternionF32`, and `QuaternionF16`
are provided, as well as the constants `imx`, `imy`, and `imz`, and (for copy-paste convenience) the
aliases 𝐢, 𝐣, and 𝐤 (as Unicode bold characters):
```jldoctest example
julia> QuaternionF64
QuaternionF64 (alias for Quaternion{Float64})
julia> 0.1 + 2.3imx + 4.5imz
0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤
julia> 0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤
0.1 + 2.3𝐢 + 0.0𝐣 + 4.5𝐤
```
As with the complex `im`, the result of multiplying `imx`, etc., with any real number will be a
quaternion with the type of the other number.

[^1]:
    Note that, mathematically speaking, quaternions can only be defined over a
    [field](https://en.wikipedia.org/wiki/Field_(mathematics)#Definition), which necessarily cannot
    be an integer type (because the multiplicative inverse of an integer is not generally an
    integer).  Nonetheless, it is possible to define a `Quaternion{<:Integer}`, which should behave
    as expected.  However, many functions (such as [`exp`](@ref), [`log`](@ref), etc.)  will then
    return a `Quaternion` of some different type, just as is the case for `Complex{<:Integer}`.

Components of the quaternion are stored as a four-element static array (even for `QuatVec`):
```jldoctest example
julia> components(q)
4-element StaticArraysCore.SVector{4, Float64} with indices SOneTo(4):
 1.0
 2.0
 3.0
 4.0
```
Those components can be indexed directly, just like an ordinary array:
```jldoctest example
julia> q[1], q[2], q[3], q[4]
(1.0, 2.0, 3.0, 4.0)
julia> q[2:4]
3-element Vector{Float64}:
 2.0
 3.0
 4.0
julia> q[[3, 2]]
2-element Vector{Float64}:
 3.0
 2.0
```
For convenience, the scalar and vector components can also be accessed in analogy with complex
numbers as
```jldoctest example
julia> real(q)
1.0
julia> imag(q)
3-element Vector{Float64}:
 2.0
 3.0
 4.0
```
Alternatively, *and slightly less efficiently*, various parts can be accessed as fields:
```jldoctest example
julia> q.w, q.x, q.y, q.z
(1.0, 2.0, 3.0, 4.0)
julia> q.vec
3-element Vector{Float64}:
 2.0
 3.0
 4.0
julia> q.re
1.0
julia> q.im
3-element Vector{Float64}:
 2.0
 3.0
 4.0
```
Again, however, these field accesses incur a slight overhead, so it's more efficient to treat the
quaternion as an array and use indexing.

Functions may also be broadcast to *each component* of a `Quaternion`.  For example, this can be
particularly helpful when simplifying `Symbolics` expressions:
```jldoctest symbolics
julia> @variables q[1:4];  # Defines q[1] through q[4] as symbolic variables

julia> Q = Quaternion(q...);

julia> simplify.(Q * imz * conj(Q))
0 + (2q[1]*q[3] + 2q[2]*q[4])𝐢 + (2q[3]*q[4] - 2q[1]*q[2])𝐣 + (q[1]^2 + q[4]^2 - (q[2]^2) - (q[3]^2))𝐤
```

The basic algebraic operations work as you would expect:
```jldoctest example
julia> p + q
5.0 + 5.0𝐢 + 5.0𝐣 + 5.0𝐤
julia> p - q
3.0 + 1.0𝐢 - 1.0𝐣 - 3.0𝐤
julia> p * q
-12.0 + 16.0𝐢 + 4.0𝐣 + 22.0𝐤
julia> q * p  # Note the non-commutativity
-12.0 + 6.0𝐢 + 24.0𝐣 + 12.0𝐤
julia> q / p
0.6666666666666666 + 0.3333333333333333𝐢 + 0.0𝐣 + 0.6666666666666666𝐤
```
Essential mathematical functions familiar from complex math, such as [`conj`](@ref), [`abs`](@ref),
[`abs2`](@ref), [`log`](@ref), [`exp`](@ref), etc., are also available.


## Contents

```@contents
Depth = 4
```

## Function list

The following list contains the public functions inside the `Quaternionic` module.  Note that there
are also many standard math functions defined for `Quaternion`s that live in the `Base` module, as
noted above.

```@index
Modules = [Quaternionic]
```
