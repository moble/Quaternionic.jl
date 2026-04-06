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
𝐯 𝐰 = 𝐯 \cdot 𝐰 + 𝐯 \wedge 𝐰,
```
where the dot product is the usual scalar inner product and the wedge product
is the antisymmetric outer product (the [exterior
product](https://en.wikipedia.org/wiki/Exterior_algebra)).  The geometric
product is linear, associative, and distributive, and satisfies
```math
𝐯 𝐯 = \| 𝐯 \|^2.
```
Two key consequences follow immediately.  Parallel vectors commute:
``𝐯 𝐰 = 𝐰 𝐯`` when ``𝐯 \parallel 𝐰``.  Orthogonal vectors anticommute:
``𝐯 𝐰 = -𝐰 𝐯`` when ``𝐯 \perp 𝐰``.

The full algebra over ℝ³ has dimension ``2^3 = 8``.  A basis is provided by
products of the Cartesian basis vectors, grouped by *grade* (number of
vector factors):
```math
\begin{array}{ll}
\text{grade 0 (scalar):}   & \boldsymbol{1}, \\[4pt]
\text{grade 1 (vectors):}  & 𝐱,\; 𝐲,\; 𝐳, \\[4pt]
\text{grade 2 (bivectors):}& 𝐱𝐲,\; 𝐱𝐳,\; 𝐲𝐳, \\[4pt]
\text{grade 3 (pseudoscalar):} & 𝐈 = 𝐱𝐲𝐳.
\end{array}
```
Each bivector squares to ``-1``.  For example,
```math
𝐱𝐲𝐱𝐲 = -𝐱𝐲𝐲𝐱 = -𝐱(𝐲𝐲)𝐱 = -𝐱𝐱 = -1.
```
The pseudoscalar ``𝐈 = 𝐱𝐲𝐳`` also squares to ``-1``, and its
inverse is ``𝐈^{-1} = 𝐳𝐲𝐱 = -𝐱𝐲𝐳``.


## The even subalgebra: Quaternions

Products of an *even* number of vectors — grades 0 and 2 — form a
closed subalgebra, the *even subalgebra*, spanned by ``\{𝟏, 𝐱𝐲,
𝐱𝐳, 𝐲𝐳\}``.  This four-dimensional algebra is precisely the
quaternions.[^1]

[^1]: Elements built from an *odd* number of vectors produce reflections
      (rather than rotations) when used in the sandwich product ``𝐐\, 𝐯\,
      𝐐^{-1}``.  This is why rotations are represented exclusively by
      even-grade elements.  See [DoranLasenby_2010](@citet) for details.

### Assigning ``𝐢``, ``𝐣``, ``𝐤`` to bivectors

The familiar quaternion units ``𝐢``, ``𝐣``, ``𝐤`` each square to ``-1``
and satisfy ``𝐢𝐣𝐤 = -\boldsymbol{1}``.  One natural bivector assignment is
``𝐢 = 𝐲𝐳``, ``𝐣 = 𝐳𝐱``, ``𝐤 = 𝐱𝐲``, but this gives ``𝐢𝐣 = (𝐲𝐳)(𝐳𝐱) =
𝐲𝐱 = -𝐤``, so ``𝐢𝐣𝐤 = +\boldsymbol{1}`` — the wrong sign.

The correct assignment uses the *reversed* bivectors:
```math
\begin{aligned}
𝐢 &= 𝐳𝐲 = -𝐲𝐳, \\
𝐣 &= 𝐱𝐳 = -𝐳𝐱, \\
𝐤 &= 𝐲𝐱 = -𝐱𝐲.
\end{aligned}
```
These are equivalently the [Hodge
duals](https://en.wikipedia.org/wiki/Hodge_star_operator) of the Cartesian
basis vectors under the pseudoscalar inverse:
```math
𝐢 = 𝐈^{-1}𝐱, \qquad
𝐣 = 𝐈^{-1}𝐲, \qquad
𝐤 = 𝐈^{-1}𝐳.
```
With this convention one can verify the standard rules:
```math
𝐢^2 = 𝐣^2 = 𝐤^2 = 𝐢𝐣𝐤 = -\boldsymbol{1},
```
and the cyclic multiplication table
```math
𝐢𝐣 = 𝐤, \qquad 𝐣𝐤 = 𝐢, \qquad 𝐤𝐢 = 𝐣.
```
We will see below that this choice ensures ``𝐢`` generates right-handed
rotations about ``𝐱``, ``𝐣`` about ``𝐲``, and ``𝐤`` about ``𝐳``.

### Components and storage

A general quaternion is written
```math
𝐐 = w\,\boldsymbol{1} + x\,𝐢 + y\,𝐣 + z\,𝐤,
```
and stored as the tuple `(w, x, y, z)`, matching the component names used
throughout the package.  The norm is the standard Euclidean norm on ℝ⁴:
```math
\| 𝐐 \| = \sqrt{w^2 + x^2 + y^2 + z^2}.
```


## The reverse and the quaternion conjugate

The *reverse* ``\widetilde{𝐐}`` of a geometric-algebra element reverses the
order of every vector factor in each basis blade.  For a grade-0 element
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

Two important consequences follow.  First, the norm squared is recovered by
```math
𝐐\,\widetilde{𝐐} = w^2 + x^2 + y^2 + z^2 = \| 𝐐 \|^2,
```
a pure scalar.  Second, every nonzero quaternion has an inverse
```math
𝐐^{-1} = \frac{\widetilde{𝐐}}{\| 𝐐 \|^2} = \frac{\texttt{conj}(𝐐)}{\texttt{abs2}(𝐐)}.
```

!!! note "Real vs. complex quaternions"
    For real quaternions (``w, x, y, z \in \mathbb{R}``), the norm above
    equals both the *Euclidean* norm ``\sum |z_i|^2`` and the *spinor* norm
    ``\sum z_i^2``.  For complex quaternions these coincide only when all
    components are real or purely imaginary; the distinction matters for the
    spacetime algebra and Lorentz boosts.


## Rotors and the double cover of SO(3)

A *rotor* is a unit quaternion — a quaternion with ``\| 𝐑 \| = 1``, so that
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

There is a natural vector-space isomorphism between ℝ³ and the space of
pure-vector quaternions (those with ``w = 0``):
```math
𝐱 \leftrightarrow 𝐢, \qquad
𝐲 \leftrightarrow 𝐣, \qquad
𝐳 \leftrightarrow 𝐤.
```
Via the Hodge-duality relations ``𝐢 = 𝐈^{-1}𝐱``, etc., this isomorphism is
in fact an *algebra* isomorphism (not merely a vector-space one).  It allows
us to treat vectors as if they were pure quaternions, and vice versa — a
convention used extensively throughout this package.

Under this identification, a vector ``𝐯 = v_x 𝐱 + v_y 𝐲 + v_z 𝐳`` is
represented by the pure quaternion ``v_x 𝐢 + v_y 𝐣 + v_z 𝐤``.

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

As a concrete check, ``𝐤 = \exp(\rho/2\; 𝐤)`` at ``\rho = \pi/2`` gives the
rotor for a 90° rotation about ``𝐳``:
```math
𝐑 = \exp\!\left(\frac{\pi}{4} 𝐤\right)
  = \frac{1}{\sqrt{2}} + \frac{1}{\sqrt{2}}\, 𝐤.
```
Applying it to ``𝐱`` (i.e., `imx`):
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
