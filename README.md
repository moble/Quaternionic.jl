# Quaternionic

[![Test Status](https://github.com/moble/Quaternionic.jl/workflows/tests/badge.svg)](https://github.com/moble/Quaternionic.jl/actions)
[![Test Coverage](https://codecov.io/gh/moble/Quaternionic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/moble/Quaternionic.jl)
[![Documentation
Status](https://github.com/moble/Quaternionic.jl/workflows/docs/badge.svg)](https://moble.github.io/Quaternionic.jl/dev)

*Quaternions for Julia*

## Installation
```julia
using Pkg
Pkg.add("Quaternionic")
```

## Usage
See [the documentation](https://moble.github.io/Quaternionic.jl/dev) for details.

## Overview
This is a simple project to build a quaternion type in Julia.  I am unhappy with certain choices
made in [the most popular](https://github.com/JuliaGeometry/Quaternions.jl) quaternion package in
Julia, and I have a lot of experience programming quaternions (including [this popular python
package](https://github.com/moble/quaternion), and [this newer and fancier
package](https://github.com/moble/quaternionic)), but I need a project to help me learn Julia
itself, so this seems like a good fit.  In particular, I am interested in understanding how a
general geometric algebra would be coded, so I am experimenting with static arrays, even though I
can imagine that hard-coding the four components would be faster.
