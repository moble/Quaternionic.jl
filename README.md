# Quaternionic

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/Quaternionic.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/Q/Quaternionic.html

<a href="https://zenodo.org/badge/latestdoi/375490468"><img align="right" hspace="0" alt="Latest DOI" src="https://zenodo.org/badge/375490468.svg"></a>
<a href="https://moble.github.io/Quaternionic.jl/stable/"><img align="right" hspace="0" alt="Stable docs" src="https://img.shields.io/badge/docs-stable-blue.svg"></a>
[![Documentation
Status](https://github.com/moble/Quaternionic.jl/workflows/docs/badge.svg)](https://moble.github.io/Quaternionic.jl/dev)
[![Test Status](https://github.com/moble/Quaternionic.jl/workflows/tests/badge.svg)](https://github.com/moble/Quaternionic.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
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
package](https://github.com/moble/quaternionic)), so transferring that experience to
Julia feels like the right thing to do.
