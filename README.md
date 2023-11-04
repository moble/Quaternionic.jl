# Quaternionic

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/Quaternionic.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/Quaternionic.html

<a href="https://zenodo.org/badge/latestdoi/375490468"><img align="right" hspace="0" alt="Latest DOI" src="https://zenodo.org/badge/375490468.svg"></a>
<a href="https://moble.github.io/Quaternionic.jl/stable/"><img align="right" hspace="0" alt="Stable docs" src="https://img.shields.io/badge/docs-stable-blue.svg"></a>
[![Documentation
Status](https://github.com/moble/Quaternionic.jl/workflows/docs/badge.svg)](https://moble.github.io/Quaternionic.jl/dev)
[![Test Status](https://github.com/moble/Quaternionic.jl/workflows/tests/badge.svg)](https://github.com/moble/Quaternionic.jl/actions)
[![Test Coverage](https://codecov.io/gh/moble/Quaternionic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/moble/Quaternionic.jl)
[![PkgEval][pkgeval-img]][pkgeval-url]

*Quaternions for Julia*

## Installation
```julia
using Pkg
Pkg.add("Quaternionic")
```

## Usage
See [the documentation](https://moble.github.io/Quaternionic.jl/dev) for details.

## Motivation
The goal of this project is to build a complete quaternion type in Julia.  There are
several features that I consider to be fairly basic to a good quaternion package, including

  * Subtypes for rotors and pure-vector quaternions, with corresponding specialized methods
  * Numerous conversions to and from other representations of things quaternions can
    represent — especially rotations
  * Smooth interpolation, and differentation of the interpolant
  * Intelligent handling of distance measures
  * Construction of random quaternions of the various special types
  * Enabling efficient integration of angular velocity
  * Construction of a minimal rotation
  * Documentation
  * Thorough testing

While there is a wide variety of [implementations of quaternions in
Julia](https://juliahub.com/ui/Search?q=quaternion&type=packages), none of them tick all
these boxes for me, and none seem easy to extend to do so.  In particular, [the most
popular](https://github.com/JuliaGeometry/Quaternions.jl) of those has made certain design
choices that conflict with my needs.

I have plenty of experience programming quaternions (including [this popular python
package](https://github.com/moble/quaternion), and [this newer and fancier
package](https://github.com/moble/quaternionic)), though not much experience with Julia, so
this seems like a good first project.  In particular, I am interested in understanding how
a general geometric algebra would be coded, so I will be experimenting with subtypes and
static arrays, even though I can imagine that hard-coding the four components could be
better in some ways.  Once I'm done experimenting, I may just rewrite that part, though
hopefully in a way that will be invisible to users.
