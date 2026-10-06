# Rotations and conversions

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
