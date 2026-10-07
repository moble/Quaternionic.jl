# Differentiating by quaternionic arguments

As with complex arguments, differentiation with respect to
quaternionic arguments treats the four components of each quaternion
as independent real arguments.  That is, a quaternion is regarded as a
vector in ``\mathbb{R}^4``, and derivatives are the usual Jacobians of
the resulting real functions.  The precise statement of this
convention, and its consequences for `Rotor` and `QuatVec`, are given
[below](@ref "Tangent conventions").

Several automatic-differentiation (AD) packages can differentiate code
that uses this package, but they obtain derivatives in different ways,
return gradients of different types, and have different limitations.
These are described in the following sections.  For most purposes,
[`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/) is the
simplest and most thoroughly tested choice for functions of a few real
variables, while the reverse-mode packages are preferable for
gradients of scalar functions of many variables.

As with complex differentiation, there are numerous notions of
quaternionic differentiation — including generalizations of the
holomorphic and Wirtinger derivatives, as well as left- and
right-multiplicative derivatives.  The goal here is to provide the
basic real derivatives from which any of these can be built.

!!! warning "Check the convention"
    If you need one of the specifically quaternionic derivatives
    mentioned above, check carefully how the real-component convention
    described here translates into your notion of derivative.  Getting
    this wrong will quietly give unexpected results.


## Supported AD packages

The following table lists the AD packages that are tested with this
package, how each of them obtains derivatives of quaternionic
functions, and the type of the gradient it returns for a scalar
function of a quaternionic argument.  All hand-written differentiation
rules live in a single extension that is loaded with
[`ChainRulesCore`](https://juliadiff.org/ChainRulesCore.jl/stable/);
the other packages differentiate the source code of this package
directly.

| Package | How derivatives are obtained | Gradient with respect to `Quaternion`, `Rotor`, and `QuatVec` arguments |
|:--------|:-----------------------------|:-----------------------------------------------------------------------|
| [ForwardDiff](https://juliadiff.org/ForwardDiff.jl/) | Differentiates the source with dual-number components.  An extension makes `ForwardDiff.derivative` return quaternions for quaternion-valued functions, including those with complex components. | Arguments must be real arrays (for example, `components(q)`), so gradients are arrays.  `ForwardDiff.derivative` of a quaternion-valued function returns a `Quaternion` (a `QuatVec` for a `QuatVec`-valued function, and a `Quaternion` for a `Rotor`-valued function). |
| [ReverseDiff](https://juliadiff.org/ReverseDiff.jl/) | Differentiates the source with tracked components.  It does not use the ChainRules rules. | Arguments must be real arrays, so gradients are arrays. |
| [Zygote](https://fluxml.ai/Zygote.jl/) and other ChainRules consumers | Uses the ChainRules rules directly; a small Zygote extension adds what ChainRules cannot express (Zygote's own rules for `Number`, `inv`, and broadcasting). | `Quaternion` for `Quaternion` and `Rotor` arguments; `QuatVec` for `QuatVec` arguments. |
| [Enzyme](https://enzymead.github.io/Enzyme.jl/stable/) | Differentiates the source natively, in forward and reverse mode, including forward-over-reverse Hessians.  An extension supplies `Enzyme.onehot` for quaternion arguments and rules for the one LAPACK call (used by `from_rotation_matrix` and `align`). | Reverse mode returns an object of the argument's own type: a `Rotor` holding the raw gradient (so `abs` reports 1), or a `QuatVec` whose hidden scalar slot may be nonzero.  Read these with `components(g)`.  `Enzyme.gradient(Forward, f, q)` returns a tuple of the partial derivatives. |
| [Mooncake](https://chalk-lab.github.io/Mooncake.jl/stable/) | Differentiates the source natively, in forward and reverse mode.  It imports the ChainRules rules for only two primitives that it cannot trace: Base's complex `sqrt` and the LAPACK eigenvector call. | Mooncake's structural `Tangent` type by default; a `Quaternion` (or a `QuatVec` for a `QuatVec` argument) with `Mooncake.Config(friendly_tangents=true)`. |
| [FastDifferentiation](https://github.com/brianguenter/FastDifferentiation.jl) | Symbolic differentiation of expression graphs built from quaternions with `Node` components. | Arrays of expressions.  Functions that branch on component values are not supported; see below. |
| [Symbolics](https://docs.sciml.ai/Symbolics/stable/) | Symbolic differentiation (for example, with `Symbolics.Differential`) of quaternions with `Num` components. | Symbolic expressions.  A few functions that branch on component values, such as `sqrt`, are not supported. |

Every package in this table except FastDifferentiation differentiates
`exp` and `log` correctly, including at the points where their
implementations switch to series expansions or special formulas: the
identity, quaternions with tiny or zero vector part, quaternions near
the negative real axis, and pure vectors.  (The exceptions for second
derivatives with Mooncake are listed under [Known limitations](@ref).)
For FastDifferentiation, these functions raise an explanatory error.


## Tangent conventions

For a `Quaternion` argument, the tangent (and cotangent) is a
`Quaternion` whose four components are the four partial derivatives,
as described in detail in the next section.  `Rotor` and `QuatVec` are
subsets of the quaternions, and their tangents follow from that fact:

| Primal type | Tangent and cotangent type in the ChainRules rules | Notes |
|:------------|:------------------------------------------|:------|
| `Quaternion{T}` | `Quaternion{float(T)}` | |
| `Rotor{T}` | `Quaternion{float(T)}`, never a `Rotor` | The cotangent is the ambient gradient in ``\mathbb{R}^4``; see below. |
| `QuatVec{T}` | `QuatVec{float(T)}` | The scalar part is zero by construction. |
| `AbstractQuaternion{Bool}` (the constants `𝐢`, `𝐣`, and `𝐤`) | `NoTangent()` | As for `Bool` in ChainRulesCore. |
| `AbstractQuaternion{<:Integer}` | `Quaternion{Float64}` or `QuatVec{Float64}` | As for `Integer` in ChainRulesCore. |
| A `Real` or `Complex` scalar receiving a quaternion cotangent | The scalar part of the cotangent | |

**`Rotor` arguments.**  A `Rotor` is a unit quaternion, so the true
tangent space at a rotor ``R`` is the three-dimensional space ``\{R\,
v : v \in \mathrm{QuatVec}\}``, which consists of the quaternions
orthogonal to ``R`` in ``\mathbb{R}^4``.  The rules in this package do
not restrict cotangents to that space.  Instead, a `Rotor` is
differentiated as the four real numbers it stores, and the gradient
with respect to a `Rotor` argument is the ambient gradient in
``\mathbb{R}^4`` of whatever the code computes from those four
numbers.  Its component along ``R`` therefore depends on how the
function is written: for example, `abs(R)` is identically 1 for a
`Rotor`, but `sqrt(sum(abs2, components(R)))` has a nonzero gradient
along ``R``.  When only the part tangent to the unit sphere is
meaningful, project the gradient ``g`` onto the tangent space with
``g-(g \cdot R)\, R``, where ``g \cdot R`` is the dot product in
``\mathbb{R}^4``; the projections of the gradients of different
implementations of the same function agree.  The gradient is returned
as a `Quaternion` (or, with Enzyme, a `Rotor` holding the raw
components), never as a renormalized `Rotor`.  Note that `abs` and
`abs2` of a `Rotor`-typed gradient return 1, whatever its components
are, so read Enzyme's `Rotor` gradients with `components(g)`.

**`QuatVec` arguments and outputs.**  A `QuatVec` has zero scalar
part, so tangents and cotangents of `QuatVec`s have zero scalar part,
and the scalar component of the input is treated as a constant.  The
ChainRules rules return `QuatVec` cotangents for `QuatVec` arguments.
Enzyme and Mooncake differentiate the stored components structurally,
so the scalar slot of their raw gradient may be nonzero when the code
reads it; that slot should be ignored (for example, by using
`vec(g)`).

**Complex components.**  Quaternions with complex components (such as
the [`Lorentz`](@ref) transformations) follow ChainRules' convention
for complex numbers: a pullback applies the conjugate transpose of the
Jacobian on ``\mathbb{C}^4``.  Quaternion multiplication is bilinear
over ``\mathbb{C}``, so the pullback of ``x \mapsto a\, x\, b`` is
``\Delta \mapsto a^\dagger\, \Delta\, b^\dagger``, where ``a^\dagger``
is the quaternion conjugate of the componentwise complex conjugate of
``a``.  For real components, ``a^\dagger = \bar{a}``, and this reduces
to the real convention.


## Simple generalization of complex differentiation

The [`ChainRulesCore`
docs](https://juliadiff.org/ChainRulesCore.jl/stable/maths/complex.html)
have this to say (and the [`Zygote`
docs](https://fluxml.ai/Zygote.jl/stable/complex/) essentially the
same thing) about differentiation with respect to complex arguments:

> `ChainRules` follows the convention that `frule` applied to a function ``f(x + i y) = u(x,y) + i v(x,y)`` with perturbation ``\Delta x + i \Delta y`` returns the value and
> ```math
> \tfrac{\partial u}{\partial x} \, \Delta x + \tfrac{\partial u}{\partial y} \, \Delta y + i \, \Bigl( \tfrac{\partial v}{\partial x} \, \Delta x + \tfrac{\partial v}{\partial y} \, \Delta y \Bigr).
> ```
> Similarly, `rrule` applied to the same function returns the value and a pullback function which, when applied to the adjoint ``\Delta u + i \Delta v``, returns
> ```math
> \Delta u \, \tfrac{\partial u}{\partial x} + \Delta v \, \tfrac{\partial v}{\partial x} + i \, \Bigl(\Delta u \, \tfrac{\partial u }{\partial y} + \Delta v \, \tfrac{\partial v}{\partial y} \Bigr).
> ```
> If we interpret complex numbers as vectors in ``\mathbb{R}^2``, then `frule` (`rrule`) corresponds to multiplication with the (transposed) Jacobian of ``f(z)``, i.e. `frule` corresponds to
> ```math
> \begin{pmatrix}
> \tfrac{\partial u}{\partial x} \, \Delta x + \tfrac{\partial u}{\partial y} \, \Delta y
> \\
> \tfrac{\partial v}{\partial x} \, \Delta x + \tfrac{\partial v}{\partial y} \, \Delta y
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{\partial u}{\partial x} & \tfrac{\partial u}{\partial y} \\
> \tfrac{\partial v}{\partial x} & \tfrac{\partial v}{\partial y} \\
> \end{pmatrix}
> \begin{pmatrix}
> \Delta x \\ \Delta y
> \end{pmatrix}
> ```
> and `rrule` corresponds to
> ```math
> \begin{pmatrix}
> \tfrac{\partial u}{\partial x} \, \Delta u + \tfrac{\partial v}{\partial x} \, \Delta v
> \\
> \tfrac{\partial u}{\partial y} \, \Delta u + \tfrac{\partial v}{\partial y} \, \Delta v
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{\partial u}{\partial x} & \tfrac{\partial u}{\partial y} \\
> \tfrac{\partial v}{\partial x} & \tfrac{\partial v}{\partial y} \\
> \end{pmatrix}^\mathsf{T}
> \begin{pmatrix}
> \Delta u \\ \Delta v.
> \end{pmatrix}
> ```

We can extend that naturally for differentiation with respect to
quaternionic arguments.  We start by working with `Quaternion`-valued
functions of a single `Quaternion` argument, and then explain how
`QuatVec` and `Rotor` relate to these rules.  Now, the statement for
quaternionic differentiation analogous to the above is:

> `Quaternionic` follows the convention that `frule` applied to a
> function
> ```math
> f(w + 𝐢 x + 𝐣 y + 𝐤 z) = s(w,x,y,z) + 𝐢 t(w,x,y,z) + 𝐣 u(w,x,y,z) + 𝐤 v(w,x,y,z)
> ```
> with perturbation ``\Delta w + 𝐢 \Delta x + 𝐣 \Delta y + 𝐤 \Delta
> z`` returns the value and
> ```math
> \begin{aligned}
> &\left(
>     \tfrac{\partial s}{\partial w} \, \Delta w + \tfrac{\partial s}{\partial x} \, \Delta x + \tfrac{\partial s}{\partial y} \, \Delta y + \tfrac{\partial s}{\partial z} \, \Delta z
> \right)
> +
> 𝐢 \left(
>     \tfrac{\partial t}{\partial w} \, \Delta w + \tfrac{\partial t}{\partial x} \, \Delta x + \tfrac{\partial t}{\partial y} \, \Delta y + \tfrac{\partial t}{\partial z} \, \Delta z
> \right) \\
> &+
> 𝐣 \left(
>     \tfrac{\partial u}{\partial w} \, \Delta w + \tfrac{\partial u}{\partial x} \, \Delta x + \tfrac{\partial u}{\partial y} \, \Delta y + \tfrac{\partial u}{\partial z} \, \Delta z
> \right)
> +
> 𝐤 \left(
>     \tfrac{\partial v}{\partial w} \, \Delta w + \tfrac{\partial v}{\partial x} \, \Delta x + \tfrac{\partial v}{\partial y} \, \Delta y + \tfrac{\partial v}{\partial z} \, \Delta z
> \right).
> \end{aligned}
> ```
> Similarly, `rrule` applied to the same function returns the value and
> a pullback function which, when applied to the adjoint ``\Delta s + 𝐢
> \Delta t + 𝐣 \Delta u + 𝐤 \Delta v``, returns
> ```math
> \begin{aligned}
> &\left(
>     \Delta s \, \tfrac{\partial s}{\partial w} + \Delta t \, \tfrac{\partial t}{\partial w} + \Delta u \, \tfrac{\partial u}{\partial w} + \Delta v \, \tfrac{\partial v}{\partial w}
> \right)
> +
> 𝐢 \left(
>     \Delta s \, \tfrac{\partial s}{\partial x} + \Delta t \, \tfrac{\partial t}{\partial x} + \Delta u \, \tfrac{\partial u}{\partial x} + \Delta v \, \tfrac{\partial v}{\partial x}
> \right) \\
> &+
> 𝐣 \left(
>     \Delta s \, \tfrac{\partial s}{\partial y} + \Delta t \, \tfrac{\partial t}{\partial y} + \Delta u \, \tfrac{\partial u}{\partial y} + \Delta v \, \tfrac{\partial v}{\partial y}
> \right)
> +
> 𝐤 \left(
>     \Delta s \, \tfrac{\partial s}{\partial z} + \Delta t \, \tfrac{\partial t}{\partial z} + \Delta u \, \tfrac{\partial u}{\partial z} + \Delta v \, \tfrac{\partial v}{\partial z}
> \right).
> \end{aligned}
> ```
> If we interpret quaternionic numbers as vectors in ``\mathbb{R}^4``,
> then `frule` (respectively, `rrule`) corresponds to multiplication
> with the Jacobian (respectively, transposed Jacobian) of ``f``.  That
> is, `frule` corresponds to
> ```math
> \begin{pmatrix}
> \tfrac{\partial s}{\partial w} \, \Delta w + \tfrac{\partial s}{\partial x} \, \Delta x + \tfrac{\partial s}{\partial y} \, \Delta y + \tfrac{\partial s}{\partial z} \, \Delta z
> \\
> \tfrac{\partial t}{\partial w} \, \Delta w + \tfrac{\partial t}{\partial x} \, \Delta x + \tfrac{\partial t}{\partial y} \, \Delta y + \tfrac{\partial t}{\partial z} \, \Delta z
> \\
> \tfrac{\partial u}{\partial w} \, \Delta w + \tfrac{\partial u}{\partial x} \, \Delta x + \tfrac{\partial u}{\partial y} \, \Delta y + \tfrac{\partial u}{\partial z} \, \Delta z
> \\
> \tfrac{\partial v}{\partial w} \, \Delta w + \tfrac{\partial v}{\partial x} \, \Delta x + \tfrac{\partial v}{\partial y} \, \Delta y + \tfrac{\partial v}{\partial z} \, \Delta z
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{\partial s}{\partial w} & \tfrac{\partial s}{\partial x} & \tfrac{\partial s}{\partial y} & \tfrac{\partial s}{\partial z}
> \\
> \tfrac{\partial t}{\partial w} & \tfrac{\partial t}{\partial x} & \tfrac{\partial t}{\partial y} & \tfrac{\partial t}{\partial z}
> \\
> \tfrac{\partial u}{\partial w} & \tfrac{\partial u}{\partial x} & \tfrac{\partial u}{\partial y} & \tfrac{\partial u}{\partial z}
> \\
> \tfrac{\partial v}{\partial w} & \tfrac{\partial v}{\partial x} & \tfrac{\partial v}{\partial y} & \tfrac{\partial v}{\partial z}
> \end{pmatrix}
> \begin{pmatrix}
> \Delta w \\ \Delta x \\ \Delta y \\ \Delta z
> \end{pmatrix}
> ```
> and `rrule` corresponds to
> ```math
> \begin{pmatrix}
> \tfrac{\partial s}{\partial w} \, \Delta s + \tfrac{\partial t}{\partial w} \, \Delta t + \tfrac{\partial u}{\partial w} \, \Delta u + \tfrac{\partial v}{\partial w} \, \Delta v
> \\
> \tfrac{\partial s}{\partial x} \, \Delta s + \tfrac{\partial t}{\partial x} \, \Delta t + \tfrac{\partial u}{\partial x} \, \Delta u + \tfrac{\partial v}{\partial x} \, \Delta v
> \\
> \tfrac{\partial s}{\partial y} \, \Delta s + \tfrac{\partial t}{\partial y} \, \Delta t + \tfrac{\partial u}{\partial y} \, \Delta u + \tfrac{\partial v}{\partial y} \, \Delta v
> \\
> \tfrac{\partial s}{\partial z} \, \Delta s + \tfrac{\partial t}{\partial z} \, \Delta t + \tfrac{\partial u}{\partial z} \, \Delta u + \tfrac{\partial v}{\partial z} \, \Delta v
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{\partial s}{\partial w} & \tfrac{\partial s}{\partial x} & \tfrac{\partial s}{\partial y} & \tfrac{\partial s}{\partial z}
> \\
> \tfrac{\partial t}{\partial w} & \tfrac{\partial t}{\partial x} & \tfrac{\partial t}{\partial y} & \tfrac{\partial t}{\partial z}
> \\
> \tfrac{\partial u}{\partial w} & \tfrac{\partial u}{\partial x} & \tfrac{\partial u}{\partial y} & \tfrac{\partial u}{\partial z}
> \\
> \tfrac{\partial v}{\partial w} & \tfrac{\partial v}{\partial x} & \tfrac{\partial v}{\partial y} & \tfrac{\partial v}{\partial z}
> \end{pmatrix}^\mathsf{T}
> \begin{pmatrix}
> \Delta s \\ \Delta t \\ \Delta u \\ \Delta v
> \end{pmatrix}.
> ```

We can easily restrict this definition to handle cases like
``\mathbb{R} \to \mathbb{H}`` and ``\mathbb{H} \to \mathbb{R}``; we
just map ``\mathbb{R}`` into ``\mathbb{H}`` by inclusion, with just
the scalar component being nonzero, so that various terms drop out of
these matrices.  For example, our function might be ``\mathbb{R} \to
\mathbb{H}``:
```math
f(w) = s(w) + 𝐢 t(w) + 𝐣 u(w) + 𝐤 v(w).
```
So the `rrule` would just look like
```math
\tfrac{\partial s}{\partial w} \, \Delta s + \tfrac{\partial t}{\partial w} \, \Delta t + \tfrac{\partial u}{\partial w} \, \Delta u + \tfrac{\partial v}{\partial w} \, \Delta v
=
\begin{pmatrix}
\tfrac{\partial s}{\partial w}
&
\tfrac{\partial t}{\partial w}
&
\tfrac{\partial u}{\partial w}
&
\tfrac{\partial v}{\partial w}
\end{pmatrix}
\begin{pmatrix}
\Delta s \\ \Delta t \\ \Delta u \\ \Delta v
\end{pmatrix}.
```
Or we might have a function ``\mathbb{H} \to \mathbb{R}``:
```math
f(w + 𝐢 x + 𝐣 y + 𝐤 z) = s(w, x, y, z),
```
so the pullback maps the scalar cotangent ``\Delta s`` to the
quaternion
```math
\tfrac{\partial s}{\partial w} \, \Delta s
+ 𝐢 \, \tfrac{\partial s}{\partial x} \, \Delta s
+ 𝐣 \, \tfrac{\partial s}{\partial y} \, \Delta s
+ 𝐤 \, \tfrac{\partial s}{\partial z} \, \Delta s,
\qquad\text{or, in matrix form,}\qquad
\begin{pmatrix}
\tfrac{\partial s}{\partial w}
\\
\tfrac{\partial s}{\partial x}
\\
\tfrac{\partial s}{\partial y}
\\
\tfrac{\partial s}{\partial z}
\end{pmatrix}
\Delta s.
```


Similarly, we can extend this with multiple arguments —
``\mathbb{R}``, ``\mathbb{H}``, or other — by appending those
arguments to the arguments of ``s``, ``t``, ``u``, and ``v``, and
similarly for multiple outputs.  For example, a function ``\mathbb{R}
\times \mathbb{H} \to \mathbb{H}`` would look like
```math
f(\sigma, w, x, y, z) =
s(\sigma, w, x, y, z)
+ 𝐢 t(\sigma, w, x, y, z)
+ 𝐣 u(\sigma, w, x, y, z)
+ 𝐤 v(\sigma, w, x, y, z).
```
The `rrule` in matrix form looks like
```math
\begin{pmatrix}
\tfrac{\partial s}{\partial \sigma} & \tfrac{\partial t}{\partial \sigma} & \tfrac{\partial u}{\partial \sigma} & \tfrac{\partial v}{\partial \sigma} \\[3pt]
\tfrac{\partial s}{\partial w}      & \tfrac{\partial t}{\partial w}      & \tfrac{\partial u}{\partial w}      & \tfrac{\partial v}{\partial w}      \\[3pt]
\tfrac{\partial s}{\partial x}      & \tfrac{\partial t}{\partial x}      & \tfrac{\partial u}{\partial x}      & \tfrac{\partial v}{\partial x}      \\[3pt]
\tfrac{\partial s}{\partial y}      & \tfrac{\partial t}{\partial y}      & \tfrac{\partial u}{\partial y}      & \tfrac{\partial v}{\partial y}      \\[3pt]
\tfrac{\partial s}{\partial z}      & \tfrac{\partial t}{\partial z}      & \tfrac{\partial u}{\partial z}      & \tfrac{\partial v}{\partial z}
\end{pmatrix}
\begin{pmatrix}
\Delta s \\ \Delta t \\ \Delta u \\ \Delta v
\end{pmatrix}.
```
If this function is just scalar multiplication, we have ``s = \sigma\,
w``, etc., so the above becomes
```math
\begin{pmatrix}
w & x & y & z \\[3pt]
\sigma & 0 & 0 & 0 \\[3pt]
0 & \sigma & 0 & 0 \\[3pt]
0 & 0 & \sigma & 0 \\[3pt]
0 & 0 & 0 & \sigma
\end{pmatrix}
\begin{pmatrix}
\Delta s \\ \Delta t \\ \Delta u \\ \Delta v
\end{pmatrix} =
\begin{pmatrix}
w \Delta s + x \Delta t + y \Delta u + z \Delta v \\
\sigma \Delta s \\ \sigma \Delta t \\ \sigma \Delta u \\ \sigma \Delta v
\end{pmatrix}.
```

Essentially, we imagine a wrapper where the quaternions on input and
output are expanded to arrays, the AD proceeds as usual, and then the
resulting arrays are reshaped back into quaternions as needed.



## Quaternions at the boundary of a differentiated function

Every package in the table above can differentiate functions that use
quaternions *internally*, in whatever combination of operations the
function needs.  The questions in this section concern only the
*boundary* of the function being differentiated: the arguments that
the AD package differentiates with respect to, and the final output.
Because `AbstractQuaternion <: Number`, a generic AD interface sees a
quaternion there as a single scalar, and would compute the derivative
of (or with respect to) its scalar part only.  This package therefore
adds methods to the AD packages so that each such case either gives
the full derivative or raises an `ArgumentError` explaining what to do
instead; none silently returns a partial result.

The pattern that works with every package is to write the function so
that its argument is a real number or a real array, from which it
builds any quaternions it needs, and its output is a real number or a
real array, such as `collect(components(q))` or `to_float_array(qs)`.
The tests of this package use that pattern throughout.

**Quaternion-valued outputs.**  For a function ``f`` of a real
variable ``t`` (or of a real vector ``x``) whose output is a
quaternion, or an array of quaternions:

| Call | Result |
|:-----|:-------|
| `ForwardDiff.derivative(f, t)` | The full derivative: a `Quaternion` (a `QuatVec` if the output is a `QuatVec`), or an array of them |
| `ForwardDiff.jacobian(f, x)` with an array of quaternions as output | A matrix of quaternions, whose entry ``(i, j)`` is the derivative of the ``i``th output with respect to ``x_j`` |
| `Enzyme.autodiff(Forward, …)` and `Enzyme.jacobian(Forward, f, x)` | The full derivatives, with quaternion entries |
| `DifferentiationInterface.derivative` and `pushforward` | The full derivative with `AutoForwardDiff`, `AutoZygote`, and `AutoEnzyme` in either mode, through an extension of this package that seeds reverse-mode pullbacks with each of the four units; see below for the exceptions |
| `DifferentiationInterface.jacobian` of an array of quaternions | A matrix of quaternions in forward mode (`AutoForwardDiff`, `AutoEnzyme(mode=Forward)`); an `ArgumentError` in reverse mode |
| `Zygote.gradient` and `Zygote.jacobian` | An `ArgumentError`, as Zygote raises for complex outputs.  An explicit pullback from `Zygote.pullback`, applied to the cotangents ``1``, ``𝐢``, ``𝐣``, and ``𝐤`` in turn, gives the four components of the derivative |
| `ReverseDiff.gradient` and `ReverseDiff.jacobian` | An `ArgumentError` |
| `Enzyme.gradient(Reverse, …)` and `Enzyme.jacobian(Reverse, …)` | An error from Enzyme, which requires real outputs in reverse mode |
| `Mooncake.value_and_gradient!!` | An error from Mooncake, which requires a floating-point output; `Mooncake.value_and_pullback!!` with an explicit Mooncake tangent works |

The exceptions with `DifferentiationInterface.derivative` are these:
with `AutoMooncakeForward` the derivative is returned as Mooncake's
structural `Tangent` rather than as a `Quaternion`; with reverse-mode
`AutoMooncake` it raises an error, because Mooncake requires the seed
of a struct-valued output to be its own tangent type; and with
reverse-mode `AutoEnzyme` a `QuatVec`-valued function raises an error,
because Enzyme requires a `QuatVec` seed while
DifferentiationInterface seeds with `oneunit(y)`, which is a
`Quaternion`.

**Quaternion arguments.**  For a real-valued function of a quaternion
argument ``q``, or of an array of quaternions:

* Reverse-mode gradients work: `Zygote.gradient`,
  `Enzyme.gradient(Reverse, …)`, and Mooncake (directly or through
  `DifferentiationInterface.gradient`) return gradients of the types
  listed in the table of [Supported AD packages](@ref).
* `Enzyme.gradient(Forward, f, q)` returns a tuple of the four partial
  derivatives for a single quaternion, and raises an `ArgumentError`
  for an array of quaternions, because Enzyme's forward-mode seeds for
  arrays would perturb only the scalar parts.
* ForwardDiff and ReverseDiff accept only real numbers and real arrays
  as arguments.
* `DifferentiationInterface.derivative` with respect to a single
  quaternion argument returns the derivative along `oneunit(q)` — that
  is, with respect to the scalar component only — just as it does for
  a `Complex` argument.  This is how DifferentiationInterface defines
  the derivative, so it is not changed here, but it is rarely what is
  wanted; use a gradient, or pass the components as a real vector.
* Forward-mode `DifferentiationInterface.gradient` and `jacobian` with
  respect to an array of quaternions raise an `ArgumentError`, for the
  same reason as with Enzyme.


## Known limitations

The following limitations are known.  Several of them are bugs in the
AD packages themselves, which have been reported or reproduced
upstream.  Each item says whether it applies *anywhere* in the
function being differentiated, or only at its *boundary* (its
arguments and final output, as described in the previous section).

**Zygote**
  * Anywhere: functions that mutate arrays cannot be differentiated:
    `squad`, `unflip`, `unflip!`, `slerp` with the `unflip=true`
    keyword, and `to_euler_phases!`.
  * Anywhere: `to_float_array` and `from_float_array` use
    `reinterpret`, for which Zygote has no adjoint.
  * Anywhere: products of matrices of quaternions are not supported,
    because ChainRules restricts its array rules to commutative number
    types.
  * Anywhere: `Zygote.hessian` through `from_rotation_matrix` is not
    tested.
  * Boundary only: `Zygote.jacobian` fails for functions that splat
    their vector argument, such as `x -> from_euler_angles(x...)`,
    although `Zygote.gradient` works.  This is a limitation of Zygote
    with any splatted argument, not only with quaternions.

**Enzyme**
  * Anywhere: `align` applied to vectors of `QuatVec`s or `Rotor`s,
    and the forward-over-reverse Hessian of `sqrt`, need
    `Enzyme.set_runtime_activity` on the mode.
  * Anywhere: batched reverse mode (batch width greater than 1, which
    DifferentiationInterface uses for Jacobians) crashes the compiler
    on `squad` and `unflip`.  Gradients of scalar functions (batch
    width 1) work.

**Mooncake**
  * Anywhere: forward-over-reverse Hessians are affected by two bugs
    in Mooncake.  The Hessian of a product with a factor that is
    exactly zero is silently wrong, which affects `log` at quaternions
    with zero vector part (including the identity).  The Hessian of
    `atan(y, x)` at `x = 0` throws an error about `bitcast`, which
    affects `log` and real powers of pure-vector quaternions.  First
    derivatives in both modes are correct.
  * Boundary only: `Mooncake.TestUtils.test_rule` cannot be used with
    `QuatVec` arguments, because its finite differences perturb the
    hidden scalar slot.

**ReverseDiff**
  * Anywhere: compiled or replayed tapes freeze value-dependent
    branches at the point where the tape was recorded.  The functions
    in this package switch between formulas (for example, to series
    expansions near the identity) at thresholds where the formulas
    agree to rounding error, so a tape recorded at a generic point is
    accurate everywhere except exactly at the special points: the
    identity, quaternions with zero vector part, and the negative real
    axis.  Uncompiled `ReverseDiff.gradient` re-records the tape on
    each call and has no such problem.

**FastDifferentiation**
  * Anywhere: FastDifferentiation cannot follow branches that depend
    on the values of variables.  Functions that choose an algorithm
    that way — `exp`, `log`, `sqrt`, non-integer powers, powers of
    `Rotor`s, `to_euler_phases`, and `from_euler_phases` — raise an
    error explaining this when called with FastDifferentiation
    variables.
  * Anywhere: integer powers are also disabled, because
    FastDifferentiation itself returns incorrect Jacobians for some of
    them.  The same upstream problem gives wrong Jacobians for some
    expressions in which several products share a common
    subexpression, such as the rotation matrix of a normalized rotor,
    so results from FastDifferentiation should be checked
    independently.

**Symbolics**
  * Anywhere: a few functions, such as `sqrt`, choose an algorithm by
    comparing the values of components, and raise an error
    ("non-boolean (Num) used in boolean context") for symbolic
    components.


## For AD and extension authors

The functions below help code that must work both with ordinary
numbers and with the number types of AD packages, such as
ForwardDiff's `Dual` and ReverseDiff's `TrackedReal`.  The extensions
for those packages add methods to `value`, so that thresholds and
special cases in this package are chosen by the values of the
components rather than by their derivative parts.

```@docs
Quaternionic.value
Quaternionic.iszerovalue
```

