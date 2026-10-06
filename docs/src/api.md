# API index

This page lists every documented function, type, and constant that the package
exports or declares public, with links to their documentation on the other pages.
Internal helper functions are not listed.

Some of these functions are defined in `Quaternionic`, while others extend functions
from Julia's `Base` and `LinearAlgebra` modules, such as `conj`, `abs`, `exp`, `log`,
`randn`, `normalize`, `⋅` (also known as `dot`), and `×` (also known as `cross`).
Many other functions from `Base` that are not listed here have also been extended to
work with quaternionic types, so that quaternions can generally function as numbers.
These include `+`, `-`, `*`, `/`, `inv`, `==`, `isequal`, `isnan`, `isinf`, `iszero`,
`isone`, `show`, `read`, `write`, `hash`, `promote_rule`, and so on.  They are not
separately documented, but should behave analogously to their behavior with `Complex`
numbers.

A few exported names are documented together with closely related names: `distance2`
with [`distance`](@ref); the type aliases `QuaternionF64`, `QuaternionF32`, and
`QuaternionF16` with [`Quaternion`](@ref), and similarly for [`Rotor`](@ref) and
[`QuatVec`](@ref).  Names that are public but not
exported, such as [`Quaternionic.value`](@ref) and [`Quaternionic.KAN`](@ref), must be
qualified with the module name or imported explicitly.

```@index
Modules = [Quaternionic, Base, LinearAlgebra]
```
