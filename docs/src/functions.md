# Function list

Functions that operate on `Quaternion`, `Rotor`, and `QuatVec` types
are defined both in the `Quaternionic` module and in Julia's own
`Base` and `LinearAlgebra` modules.  The functions defined in the
`Base` module are generally standard mathematical functions that have
been extended to work with quaternions, while the functions defined in
the `Quaternionic` module are more specialized functions that are
specific to quaternion algebra.

## External functions

  * [`Base.conj`](@ref)
  * [`Base.abs`](@ref)
  * [`Base.abs2`](@ref)
  * [`Base.angle`](@ref)
  * [`Base.exp`](@ref)
  * [`Base.log`](@ref)
  * [`Base.sqrt`](@ref)
  * [`Base.:^`](@ref)
  * [`LinearAlgebra.:⋅`](@ref)

In addition, many other functions from `Base` that are not listed here
have also been extended to work with quaternionic types, so that
quaternions can generally function as numbers.  This includes
functions such as `+`, `-`, `*`, `/`, `inv`, `==`, `isequal`, `isnan`,
`isinf`, `iszero`, `isone`, `show`, `read`, `write`, `hash`,
`promote_rule`, and so on.  These are not separately documented, but
should behave analogously to their behavior with `Complex` numbers.

## Functions defined in `Quaternionic`

```@index
Modules = [Quaternionic]
```
