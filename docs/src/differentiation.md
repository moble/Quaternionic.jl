# Differentiating by quaternionic arguments

As with complex arguments, differentiation with respect to quaternionic
arguments treats the components of the quaternionic argument as independent real
arguments.  These rules are implemented for this package in `ChainRulesCore`,
which means that they should work seamlessly with [any package that relies on
`ChainRulesCore`](https://juliadiff.org/ChainRulesCore.jl/stable/#ChainRules-roll-out-status),
such as [`Zygote`](https://github.com/FluxML/Zygote.jl).  Derivatives can also
be calculated automatically using
[`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/)

As with complex differentiation, there are numerous notions of
quaternionic differentiation — including generalizations of the
holomorphic and Wirtinger derivatives, as well as left- and
right-multiplicative derivatives.  The goal here is to provide the
basic differentiation rules upon which these derivatives can be
implemented, but not to implement those derivatives themselves. It is
recommended that you carefully check how the definitions of `frule`
and `rrule` translate into your specific notion of quaternionic
derivatives, since getting this wrong will quietly give you wrong
results.


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
Similarly, we can extend this with multiple arguments —
``\mathbb{R}``, ``\mathbb{H}``, or other — by appending those
arguments to ``s``, ``t``, ``u``, and ``v``, for example.  And
equivalently for the outputs.

Essentially, we imagine a wrapper where the quaternions on input and
output are expanded to arrays, the AD proceeds as usual, and then the
resulting arrays are reshaped back into quaternions as needed.


## Applications to `QuatVec` and `Rotor`

To understand how this works for `QuatVec` and `Rotor` inputs or
outputs, we just consider that these are submanifolds of the
`Quaternion` manifold.  The only subtlety is that — while the tangent
space to `Quaternion` and `QuatVec` are naturally identified with
`Quaternion` and `QuatVec` themselves — the tangent space of the
`Rotor` submanifold is naturally identified with `Quaternion`.

Thus, for a `QuatVec` input, ``w`` must always be 0, which means that
the tangent must always have ``\Delta w = 0``, and we always treat the
output functions ``(s,t,u,v)`` as independent of ``w`` so that
``\partial s / \partial w`` and so on are always 0.  Similarly, for
`QuatVec` outputs, ``s`` must always be 0, so that the tangent must
always have ``\Delta s = 0``, and ``\partial s / \partial w`` and so
on are always 0.  With these considerations in mind, it's not hard to
simplify the expressions above for `QuatVec` inputs and outputs.

On the other hand, because the tangent space to the `Rotor`
submanifold is naturally identified with `Quaternion`, while there is
a natural constraint on the norms of the input and output arguments,
there are no structural constraints on the tangent vectors (just that
they must be orthogonal to the arguments themselves, when viewed as
vectors in ``\mathbb{R}^4`` with the Euclidean metric).  Thus, the
expressions above for `Quaternion` inputs and outputs will look
formally identical for `Rotor` inputs or outputs.


## Older functions

In this vein, we also have some very explicit functions for computing "primals"
(values) and derivatives of functions of `log` and `exp`.  These are older, and
likely to be deprecated at some point in favor of `ChainRulesCore`-based AD.
Also, because of massive simplifications that result when using the right types,
these derivatives are more strict about input types than the main functions
themselves.  For example, the derivatives of `exp` are defined only for
`QuatVec` arguments; the derivatives of `log` are defined only for `Rotor`
arguments; etc.


```@autodocs
Modules = [Quaternionic]
Pages   = ["gradients_exp_log.jl"]
```
