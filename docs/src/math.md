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

Matrices and vectors of quaternions work with LinearAlgebra's generic
algorithms, as long as those algorithms keep the factors of every
product in order.  Because `x ⋅ y` is `conj(x) * y` for quaternions,
the conjugate transpose `x'` behaves as it does for complex numbers:
`x' * y` is `sum(conj.(x) .* y)` for vectors `x` and `y`.  Products of
matrices and vectors, `inv`, solutions with `\` and `/`, `lu`,
`cholesky`, `qr`, least-squares solutions, `svd`, `svdvals`, `opnorm`,
`cond`, `rank`, `pinv`, and `eigen` and `eigvals` of `Hermitian`
matrices all give the quaternionic results.  (The singular-value and
Hermitian eigenvalue routines come from GenericLinearAlgebra, which
this package loads.)

A few functions assume that the entries of a matrix commute, and would
silently give wrong results for quaternions, so they throw an
`ArgumentError` instead:

  * `det`, `logdet`, and `logabsdet`.  A determinant of quaternion
    matrices, such as the Study determinant, may be added later; see
    [issue #119](https://github.com/moble/Quaternionic.jl/issues/119).
  * `transpose(x) * A`, `transpose(x) / A`, `transpose(x) * A * y`,
    and `muladd(transpose(x), A, z)` for a vector `x`, and
    `inv(transpose(A))`.  For quaternion matrices, the transpose of a
    product is not the product of the transposes in reverse order, but
    LinearAlgebra computes these as if it were.  Use `permutedims(x)`
    in place of `transpose(x)`, or the conjugate transpose `x'`, which
    does reverse products.
  * `inv` of a `Symmetric` matrix, whose inverse need not be symmetric
    for quaternions.  Use `inv(Matrix(A))`.
  * Solving with the factorization `lu!(A)` of a `Tridiagonal` matrix,
    and `ldlt` of a `SymTridiagonal` matrix.  Use `lu(A)` or `A \ b`
    instead.  Systems with `Tridiagonal` and `SymTridiagonal` matrices
    are otherwise solved correctly, with a dense factorization, which
    costs O(n³) rather than the O(n) of LinearAlgebra's specialized
    solvers.

Functions such as `exp`, `log`, and `sqrt` of quaternion matrices, and
`eigen` and `eigvals` of non-`Hermitian` ones, are not supported;
issue #119 describes how they could be.  `Statistics.cov` and
`Statistics.cor` of two quaternion vectors are not meaningful either,
because they assume that multiplication commutes.
