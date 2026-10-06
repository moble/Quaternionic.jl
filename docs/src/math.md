# Algebra and mathematical functions

Along with the basic binary operators, the essential mathematical
functions like [`conj`](@ref), [`abs`](@ref), [`abs2`](@ref),
[`exp`](@ref), [`log`](@ref), etc., are implemented.  Most of these
functions are found in the `Base` module, and are simply overloaded
methods of functions that should also be familiar from `Complex`
types.  Note that we use a slightly different interpretation of
[`angle`](@ref) for `Quaternion`, compared to `Complex`.  We also have
[`absvec`](@ref) and [`abs2vec`](@ref), which are not useful in a
`Complex` context, but compute the relevant quantities for the
"vector" component of a `Quaternion`.

Quaternions with complex components, which represent [Lorentz
transformations](@ref), are also supported by most of these functions.
For them, `abs2` and `abs` compute the complex *spinor norm*, ``w^2 +
x^2 + y^2 + z^2``, and its square root, rather than the Euclidean
norm; see [Quaternions and the Spacetime Algebra](@ref) for the
reasons.

Calling a quaternion `Q` on a `QuatVec` `v`, as in `Q(v)`, computes
the "sandwich" `Q * v * conj(Q)` efficiently.  For a `Rotor`, this is
the rotation of `v`; for a general `Quaternion`, the result is also
scaled by `abs2(Q)`.

```@autodocs
Modules = [Quaternionic]
Pages   = ["algebra.jl", "math.jl"]
```

## `@fastmath`

The arithmetic operators `+`, `-`, `*`, and `/`, and the power `^`,
can be used inside `@fastmath` expressions, and they return the same
types as without `@fastmath`.  The one exception is a chained
expression such as `@fastmath 1.0 * 2.0 * 3.0 * imz`, in which none of
the first three factors is a quaternion; there, Julia's generic
fallback promotes all of the factors to a common type, so the result
is a `Quaternion` rather than a `QuatVec`.

## Matrices of quaternions

Arrays of quaternions work with most generic code, but some generic
linear-algebra and statistics routines assume that multiplication of
scalars commutes, and give silently wrong answers for quaternions.
Products of matrices and vectors, `inv`, solutions of square systems
with `\` and `/`, `lu`, `cholesky`, and operations on `Diagonal`,
`Bidiagonal`, and triangular matrices preserve the order of products,
and are correct.  On the other hand, `det` and `logabsdet` return a
quaternion whose phase is arbitrary (only `abs(det(A))` is
meaningful), solving a system with a `Tridiagonal` matrix gives wrong
results, and `Statistics.cov` and `Statistics.cor` are not meaningful
for quaternions.  Functions such as `exp`, `log`, `sqrt`, and `eigen`
of a quaternion matrix are not supported.
