# Quaternions as Geometric Algebra

Quaternions arise naturally as the *even subalgebra* of the [geometric
algebra](https://en.wikipedia.org/wiki/Geometric_algebra) over
three-dimensional Euclidean space.  This page derives the conventions used
throughout the package from first principles — starting from the geometric
product — so that every sign choice, basis assignment, and rotation formula
has a clear geometric justification.

The treatment here focuses on *real* quaternions.  The extension to complex
quaternions (needed for Lorentz boosts via the spacetime algebra) is covered
on a separate page.


## The geometric algebra over ℝ³

We start with the standard right-handed Cartesian coordinate system and its
unit basis vectors ``(𝐱, 𝐲, 𝐳)``.  The *geometric product* of two
vectors ``𝐯`` and ``𝐰`` is defined by
```math
𝐯 𝐰 = 𝐯 ⋅ 𝐰 + 𝐯 ∧ 𝐰,
```
where the dot product is the usual scalar inner product and the wedge product
is the antisymmetric outer product (the [exterior
product](https://en.wikipedia.org/wiki/Exterior_algebra)).  The geometric
product is linear, associative, and distributive, and satisfies
```math
𝐯 𝐯 = 𝐯 ⋅ 𝐯.
```
Two key consequences follow immediately.  Parallel vectors commute:
``𝐯 𝐰 = 𝐰 𝐯`` when ``𝐯 ∥ 𝐰``.  Orthogonal vectors anticommute:
``𝐯 𝐰 = -𝐰 𝐯`` when ``𝐯 ⟂ 𝐰``.

The full algebra over ℝ³ has dimension ``2^3 = 8``.  A basis is provided by
products of the Cartesian basis vectors, grouped by *grade* (number of
vector factors):
```math
\begin{array}{ll}
\text{grade 0 (scalar):}       & \boldsymbol{1}, \\[4pt]
\text{grade 1 (vectors):}      & 𝐱,\; 𝐲,\; 𝐳, \\[4pt]
\text{grade 2 (bivectors):}    & 𝐱𝐲,\; 𝐱𝐳,\; 𝐲𝐳, \\[4pt]
\text{grade 3 (pseudoscalar):} & 𝐈 = 𝐱𝐲𝐳.
\end{array}
```
Each bivector squares to ``-1``.  For example,
```math
𝐱𝐲𝐱𝐲 = -𝐱𝐲𝐲𝐱 = -𝐱(𝐲𝐲)𝐱 = -𝐱𝐱 = -1.
```
The pseudoscalar ``𝐈 = 𝐱𝐲𝐳`` also squares to ``-1``, so its
inverse is ``𝐈^{-1} = -𝐈 = 𝐳𝐲𝐱``.

### Duality

``𝐈`` has a special property in three dimensions: since moving it past any
grade-1 vector costs ``(-1)^{n-1} = (-1)^2 = +1`` sign changes, it commutes
with every element of the algebra.  Left- and right-multiplication by ``𝐈``
are therefore identical, and we can unambiguously write ``𝐈\,𝐯 = 𝐯\,𝐈``.

The [Hodge dual](https://en.wikipedia.org/wiki/Hodge_star_operator) maps
grade-``k`` elements to grade-``(n-k)`` elements.[^hodge]  For grade-1 vectors
in ℝ³ the reverse is trivial (``\widetilde{𝐯} = 𝐯``), and explicit computation
gives
```math
\star 𝐱 = 𝐲𝐳, \qquad \star 𝐲 = 𝐳𝐱, \qquad \star 𝐳 = 𝐱𝐲.
```
For the grade-2 bivectors the reverse introduces a sign, and one finds
```math
\star(𝐲𝐳) = 𝐱, \qquad \star(𝐳𝐱) = 𝐲, \qquad \star(𝐱𝐲) = 𝐳.
```
Since ``𝐈^2 = -1`` in ℝ³, we have ``𝐈^{-1} = -𝐈``, so the map
``𝐯 \mapsto 𝐈^{-1}𝐯`` gives the *negatives* of the Hodge duals:
```math
𝐈^{-1}𝐱 = 𝐳𝐲 = -\star 𝐱, \qquad
𝐈^{-1}𝐲 = 𝐱𝐳 = -\star 𝐲, \qquad
𝐈^{-1}𝐳 = 𝐲𝐱 = -\star 𝐳.
```
The two natural maps ``𝐯 \mapsto \pm \star 𝐯`` correspond exactly to the two
sign conventions for quaternion multiplication.

[^hodge]: The Hodge dual is defined in general by the property
    ``a \wedge \star b = (a \mid b)\, 𝐈``, where ``(a \mid b)`` is the
    symmetric bilinear form naturally induced by the metric on grade-``r``
    elements.  For a geometric algebra this form is characterized by
    ```math
    (a \mid b) = \langle a\, \widetilde{b} \rangle_0,
    ```
    where ``\langle \cdot \rangle_0`` extracts the grade-0 (scalar) part.
    The reverse ``\widetilde{b}`` compensates for the reordering cost of
    extracting a scalar from a product of two same-grade blades: without it,
    one picks up an extra factor of ``(-1)^{r(r-1)/2}``.  Three alternative
    expressions are all equal:
    ```math
    (a \mid b)
    = \langle a\, \widetilde{b} \rangle_0
    = \langle \widetilde{a}\, b \rangle_0
    = \langle b\, \widetilde{a} \rangle_0
    = \langle \widetilde{b}\, a \rangle_0,
    ```
    because ``\langle MN \rangle_0 = \langle NM \rangle_0`` (cyclic property of
    the scalar part) and ``\langle \widetilde{M} \rangle_0 = \langle M
    \rangle_0``.  With this bilinear form, the formula
    ``\star A = \widetilde{A}\, 𝐈`` can be verified to satisfy the defining
    property for arbitrary grade and arbitrary signature ``(p, q)``.


## The even subalgebra: Quaternions

Products of an *even* number of vectors — grades 0 and 2 — form a closed
subalgebra, the *even subalgebra*, spanned by ``\{𝟏, 𝐱𝐲, 𝐱𝐳, 𝐲𝐳\}``.
This four-dimensional algebra is precisely the quaternions.[^1]

[^1]: Elements built from an *odd* number of vectors produce reflections
      (rather than rotations) when used in the sandwich product ``𝐐\, 𝐯\,
      𝐐^{-1}``.  This is why rotations are represented exclusively by
      even-grade elements.  See [DoranLasenby_2010](@citet) for details.

### The canonical bivector basis: ``𝐢``, ``𝐣``, ``𝐤``

We define the quaternion units as the duals of the Cartesian basis vectors
under ``𝐈^{-1}``:
```math
𝐢 = 𝐈^{-1}𝐱 = 𝐳𝐲, \qquad
𝐣 = 𝐈^{-1}𝐲 = 𝐱𝐳, \qquad
𝐤 = 𝐈^{-1}𝐳 = 𝐲𝐱.
```
Each of these squares to ``-1`` (as expected for a bivector), and one can
verify the standard cyclic rules:
```math
𝐢^2 = 𝐣^2 = 𝐤^2 = 𝐢𝐣𝐤 = -\boldsymbol{1},
\qquad
𝐢𝐣 = 𝐤, \quad 𝐣𝐤 = 𝐢, \quad 𝐤𝐢 = 𝐣.
```
We will see below that this choice ensures ``𝐢`` generates right-handed
rotations about ``𝐱``, ``𝐣`` about ``𝐲``, and ``𝐤`` about ``𝐳``.

Using ``𝐈`` instead of ``𝐈^{-1}`` would give the alternative assignment
``𝐢 = 𝐲𝐳``, ``𝐣 = 𝐳𝐱``, ``𝐤 = 𝐱𝐲``, which flips the sign of every
multiplication: ``𝐢𝐣𝐤 = +\boldsymbol{1}``.  This is the convention used in
much of the robotics and aerospace literature.[^2]

[^2]: The two quaternion conventions are sometimes called the *Hamilton*
      convention (``𝐢𝐣𝐤 = -1``, used here and in most mathematics and
      physics) and the *Shuster* or *JPL* convention (``𝐢𝐣𝐤 = +1``, common
      in spacecraft attitude control).  See [SommerEtAl_2018](@citet) for a
      thorough discussion of the differences and how they arise.

### Components and storage

A general quaternion is written
```math
𝐐 = w\,\boldsymbol{1} + x\,𝐢 + y\,𝐣 + z\,𝐤,
```
and stored as the tuple `(w, x, y, z)`, matching the component names used
throughout the package.


## The reverse and the quaternion conjugate

The *reverse* ``\widetilde{𝐐}`` of a geometric-algebra element reverses the
order of every vector factor in each basis element.  For a grade-0 element
the reverse is the identity; for a grade-2 bivector it introduces a sign flip
(since swapping two anticommuting vectors costs a minus sign):
```math
\widetilde{𝐢} = \widetilde{𝐳𝐲} = 𝐲𝐳 = -𝐢, \quad
\widetilde{𝐣} = -𝐣, \quad
\widetilde{𝐤} = -𝐤.
```
For a general quaternion this gives
```math
\widetilde{𝐐} = w - x\,𝐢 - y\,𝐣 - z\,𝐤,
```
which is exactly the *quaternion conjugate* `conj(Q)` in the code.

Two important consequences follow.  First,
```math
𝐐\,\widetilde{𝐐} = w^2 + x^2 + y^2 + z^2
```
is a *pure scalar* — a non-negative real number for real quaternions.  In
analogy with `abs2` for complex numbers, this product is `abs2(Q)` in the
code.  Second, every nonzero quaternion has an inverse
```math
𝐐^{-1} = \frac{\widetilde{𝐐}}{𝐐\,\widetilde{𝐐}} = \frac{\texttt{conj}(Q)}{\texttt{abs2}(Q)}.
```

!!! note "Real vs. complex quaternions"

    For real quaternions, ``𝐐\,\widetilde{𝐐}`` is a positive real scalar.
    For complex quaternions, ``𝐐\,\widetilde{𝐐}`` is still a *scalar*
    element of the algebra (no ``𝐢``, ``𝐣``, or ``𝐤`` component), but it
    is complex-valued — it is the *spinor norm* ``\sum_i z_i^2`` rather than
    the Euclidean norm ``\sum_i |z_i|^2``.  The distinction matters for
    Lorentz boosts; it is discussed on the spacetime algebra page.


## Rotors and the double cover of SO(3)

A *rotor* is a unit quaternion — a quaternion satisfying
``𝐑\,\widetilde{𝐑} = \boldsymbol{1}``.  Any unit quaternion can be written
```math
𝐑 = \exp\!\left(\frac{\rho}{2}\,\hat{𝔯}\right)
  = \cos\frac{\rho}{2} + \hat{𝔯}\,\sin\frac{\rho}{2},
```
where ``\rho`` is a real angle and ``\hat{𝔯} = \hat{r}_x 𝐢 + \hat{r}_y 𝐣 +
\hat{r}_z 𝐤`` is a unit pure-vector quaternion.

The group of unit quaternions is ``\mathrm{Spin}(3) \cong \mathrm{SU}(2)``,
topologically the 3-sphere ``\mathbb{S}^3``.  It is a *double cover* of the
rotation group ``\mathrm{SO}(3)``: both ``𝐑`` and ``-𝐑`` represent the same
rotation (as we will see below), but they are distinct quaternions
representing distinct *spinors*.


## The vector–quaternion isomorphism and rotations

### The isomorphism

The duality relations ``𝐢 = 𝐈^{-1}𝐱``, etc., provide a natural
identification between ℝ³ and the space of pure-vector quaternions:
```math
𝐱 \leftrightarrow 𝐢, \qquad
𝐲 \leftrightarrow 𝐣, \qquad
𝐳 \leftrightarrow 𝐤.
```
This is more than a vector-space isomorphism.  For any two pure quaternions
``𝐩 = p_x 𝐢 + p_y 𝐣 + p_z 𝐤`` and ``𝐪 = q_x 𝐢 + q_y 𝐣 + q_z 𝐤``, direct
expansion of the geometric product gives
```math
𝐩\,𝐪 = -(𝐩 \cdot 𝐪) + (𝐩 \times 𝐪),
```
where the dot product contributes a scalar and the cross product contributes
a pure quaternion via the same identification.  This is the exact analogue of
the geometric product for vectors, ``𝐯𝐰 = 𝐯 \cdot 𝐰 + 𝐯 \wedge 𝐰``, with
Hodge duality mapping ``𝐯 \wedge 𝐰 \leftrightarrow 𝐯 \times 𝐰``.  In
particular, the antisymmetric part gives
```math
\frac{𝐩\,𝐪 - 𝐪\,𝐩}{2} = 𝐩 \times 𝐪,
```
which is precisely the Lie bracket of the cross-product algebra on ℝ³.  The
map is therefore an isomorphism of Lie algebras.

Under this identification, a vector ``𝐯 = v_x 𝐱 + v_y 𝐲 + v_z 𝐳`` is
represented by the pure quaternion ``v_x 𝐢 + v_y 𝐣 + v_z 𝐤``, and we freely
pass between the two descriptions throughout the package.

### The rotation formula

A rotor ``𝐑`` acts on a vector ``𝐯`` (written as a pure quaternion) by the
*sandwich product*
```math
𝐯' = 𝐑\, 𝐯\, 𝐑^{-1} = 𝐑\, 𝐯\, \widetilde{𝐑}.
```
To see that this is a right-handed rotation through angle ``\rho`` about the
axis ``\hat{𝔯}``, decompose ``𝐯 = 𝐯_\parallel + 𝐯_\perp`` into parts
parallel and perpendicular to ``\hat{𝔯}``.

- ``𝐯_\parallel`` commutes with ``\hat{𝔯}`` (parallel vectors commute), so
  it commutes with ``𝐑`` and passes through unchanged.
- ``𝐯_\perp`` anticommutes with ``\hat{𝔯}`` (orthogonal vectors
  anticommute), so the two factors of ``𝐑`` act non-trivially.  Expanding:

```math
\begin{aligned}
𝐑\, 𝐯_\perp\, 𝐑^{-1}
&= \left(\cos\tfrac{\rho}{2} + \sin\tfrac{\rho}{2}\,\hat{𝔯}\right)
   𝐯_\perp
   \left(\cos\tfrac{\rho}{2} - \sin\tfrac{\rho}{2}\,\hat{𝔯}\right) \\
&= \cos^2\!\tfrac{\rho}{2}\; 𝐯_\perp
   + \sin\tfrac{\rho}{2}\cos\tfrac{\rho}{2}\,[\hat{𝔯}, 𝐯_\perp]
   - \sin^2\!\tfrac{\rho}{2}\;\hat{𝔯}\, 𝐯_\perp\, \hat{𝔯} \\
&= \cos\rho\; 𝐯_\perp + \sin\rho\; \hat{𝔯} \times 𝐯_\perp,
\end{aligned}
```

which is precisely the right-handed rotation of ``𝐯_\perp`` through angle
``\rho`` about ``\hat{𝔯}``.  Putting the two parts together:

```math
𝐑\, 𝐯\, 𝐑^{-1}
= 𝐯_\parallel + \cos\rho\; 𝐯_\perp + \sin\rho\; \hat{𝔯} \times 𝐯_\perp.
```

The presence of *two* factors of ``𝐑`` in the sandwich explains two things.
First, the rotation angle is twice the quaternion half-angle: each factor
of ``𝐑 = \cos(\rho/2) + \ldots`` contributes a half-angle, and together
they give the full angle ``\rho``.  Second, negating ``𝐑 \to -𝐑`` leaves the
sandwich unchanged, confirming the double cover.

### Rotations about coordinate axes

As a concrete check, the rotor for a 90° rotation about ``𝐳`` is
```math
𝐑 = \exp\!\left(\frac{\pi}{4} 𝐤\right)
  = \frac{1}{\sqrt{2}} + \frac{1}{\sqrt{2}}\, 𝐤.
```
Applying it to ``𝐢`` (corresponding to ``𝐱``):
```math
𝐑\, 𝐢\, 𝐑^{-1}
= \left(\tfrac{1}{\sqrt{2}} + \tfrac{1}{\sqrt{2}} 𝐤\right)
  𝐢
  \left(\tfrac{1}{\sqrt{2}} - \tfrac{1}{\sqrt{2}} 𝐤\right)
= 𝐣,
```
confirming a right-handed rotation ``𝐱 \to 𝐲`` about ``𝐳``.

### Efficient computation

In the package, the function-call syntax `R(v)` implements the sandwich
``𝐑\, 𝐯\, \widetilde{𝐑}`` efficiently with a dedicated formula — roughly
twice as fast as the three-multiplication sequence `R * v * conj(R)`:
```julia
Q * v * conj(Q) ≈ Q(v)   # both correct; Q(v) is ~2× faster
```
This applies to both `Quaternion` and `Rotor` objects acting on a `QuatVec`.

## Further reading

The geometric algebra perspective on quaternions and rotors is developed in
depth in [DoranLasenby_2010](@citet).

```@bibliography
Pages = ["geometric_algebra.md"]
Canonical = false
```
