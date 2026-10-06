# Distances and alignment

There are several ways of measuring the "distance" between two
quaternions: ``d(q_1, q_2)``.  Fundamentally, any comparison between
two quaternions ``q_1`` and ``q_2`` must make use of a binary
operation, for which there are two obvious choices: addition or
multiplication.  For either choice, we operate on ``q_1`` and the
appropriate inverse (either additive or multiplicative) of ``q_2``.
That is, ``d`` should be a function of either ``q_1 - q_2`` or
``q_1/q_2``.[^1]

[^1]: For ``q_1/q_2``, we are dealing with the *multiplicative* group
    of quaternions, which does not include 0, so we will assume that
    no quaternion involved in such a function can be 0.

Now, we also have a number of criteria we would like any distance
function to satisfy.  For any quaternions ``q_1`` and ``q_2`` and any
*unit* quaternion ``q_3``, we require

  * real-valued: ``d(q_1, q_2) \in \mathbb{R}``
  * symmetry: ``d(q_1, q_2) = d(q_2, q_1)``
  * invariance: ``d(q_3\, q_1, q_3\, q_2) = d(q_1, q_2) = d(q_1\, q_3, q_2\, q_3)``
  * identity: ``d(q_1, q_1) = 0``
  * positive-definiteness: ``d(q_1, q_2) > 0`` whenever ``q_1 ≠ q_2``

(Of course, it should be noted that these criteria all hold in the
*exact* case; when using floating-point numbers, they will likely be
violated near edge cases.)

It is not hard to see that these criteria can be satisfied by any of

  * `abs(q₁ - q₂)`
  * `abs2(q₁ - q₂)`
  * `abs(log(q₁ / q₂))`
  * `abs2(log(q₁ / q₂))`

If ``q_1`` and ``q_2`` are interpreted as rotations, we frequently
don't care about their signs, and just want the *smallest* distance
between them, for any choice of sign.  Furthermore, in the
multiplicative case, the `log` functions will involve calculation of
the `log` of the magnitudes of the quaternions, which should be 1.  In
this case, we relax the "positive-definiteness" criterion to allow
``d(q_1, q_2)`` to equal zero when ``q_1`` and ``q_2`` are related by
a nonzero scalar multiple.

For `Rotor` types, the latter two multiplicative options are most
relevant, while for other types the additive options are more
relevant.  These are the default behaviors of the `distance` and
`distance2` functions, which are documented together below.

```@autodocs
Modules = [Quaternionic]
Pages   = ["distance.jl"]
```


## Alignment

There are many ways to optimize alignment with rotations.  In
particular, we can seek the optimal rotation that takes one set of
points onto a corresponding set of points, or the optimal quaternion
that takes one set of quaternions onto a corresponding set of
quaternions.  In both cases, the "optimal" value depends on the metric
being used.  The simplest and most robust results are obtained when
the metric is the standard Euclidean metric (in the case of points),
or the magnitude of the difference (in the case of quaternions).
Here, we assume that `QuatVec`s represent points, and any other type
of quaternion should be treated as rotors.

```@autodocs
Modules = [Quaternionic]
Pages   = ["alignment.jl"]
```
