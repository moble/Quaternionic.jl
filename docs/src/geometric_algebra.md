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


## The geometric algebra over ``ℝ³``

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
``𝐯 𝐰 = -𝐰 𝐯`` when ``𝐯 ⟂ 𝐰``.  In particular, computations with
an orthonormal basis become simple exercises in counting the parity of
permutations.

The full algebra over ``ℝ³`` has dimension ``2^3 = 8``.  A basis is
provided by products of the Cartesian basis vectors, grouped by
*grade* (number of vector factors):
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

``𝐈`` has a special property in three dimensions: since moving it
past any grade-1 vector costs ``(-1)^{n-1} = (-1)^2 = +1`` sign
changes, it commutes with every element of the algebra.  Left- and
right-multiplication by ``𝐈`` are therefore identical, and we can
unambiguously write ``𝐈\,𝐯 = 𝐯\,𝐈``.

### The reverse

The *reverse* of a multivector, denoted ``\widetilde{A}``, is obtained
by reversing the order of the vector factors in each basis blade:
```math
\widetilde{e_{i_1} e_{i_2} \cdots e_{i_r}}
  = e_{i_r} \cdots e_{i_2}\, e_{i_1}.
```
Reordering the ``r`` vectors back to standard order requires
``r(r-1)/2`` transpositions, each costing a sign change, so a
grade-``r`` blade picks up a factor of ``(-1)^{r(r-1)/2}``.  In ``ℝ³``:
- Grades 0 and 1 are unchanged: ``\widetilde{𝐯} = 𝐯``.
- Grades 2 and 3 negate: for example
```math
\widetilde{𝐱𝐲} = 𝐲𝐱 = -𝐱𝐲, \qquad \widetilde{𝐈} = 𝐳𝐲𝐱 = -𝐱𝐲𝐳 = -𝐈.
```

### Duality

The [Hodge dual](https://en.wikipedia.org/wiki/Hodge_star_operator)
maps grade-``k`` elements to grade-``(n-k)`` elements.  The
formula[^1]
```math
\star b = \widetilde{b}\, 𝐈
```
holds for blades of any grade and any metric signature.  For grade-1
vectors ``\widetilde{𝐯} = 𝐯``, so right-multiplication by ``𝐈``
suffices:
```math
\star 𝐱 = 𝐱\, 𝐈 = 𝐲𝐳, \qquad
\star 𝐲 = 𝐲\, 𝐈 = 𝐳𝐱, \qquad
\star 𝐳 = 𝐳\, 𝐈 = 𝐱𝐲.
```
For grade-2 bivectors the reverse introduces a sign (e.g.,
``\widetilde{𝐲𝐳} = 𝐳𝐲 = -𝐲𝐳``), giving
```math
\star(𝐲𝐳) = (𝐳𝐲)(𝐱𝐲𝐳) = 𝐱, \qquad
\star(𝐳𝐱) = 𝐲, \qquad
\star(𝐱𝐲) = 𝐳.
```

[^1]: The Hodge dual is defined by linearity and the property that for
    grade-``r`` elements ``a`` and ``b``, we have ``a \wedge \star b =
    (a \mid b)\, 𝐈``, where ``(a \mid b)`` is the symmetric bilinear
    form given by
    ```math
    (a \mid b) = \langle a\, \widetilde{b} \rangle_0,
    ```
    with ``\langle \cdot \rangle_0`` extracting the scalar part.  The
    reverse ensures symmetry; without it one picks up a sign
    ``(-1)^{r(r-1)/2}`` from reordering.  The formula ``\star b =
    \widetilde{b}\, 𝐈`` satisfies the defining property for any grade
    and signature:
    ```math
    a \wedge \star b
    = \langle a\, (\star b) \rangle_n
    = \langle a\, \widetilde{b}\, 𝐈 \rangle_n
    = \langle a\, \widetilde{b} \rangle_0\, 𝐈
    = (a \mid b)\, 𝐈,
    ```
    where the first step uses ``a \wedge c = \langle ac \rangle_n``
    when the grades of ``a`` and ``c`` sum to ``n``, and the third
    uses the fact that right-multiplication by ``𝐈`` shifts grade by
    ``n``, so only the scalar part of ``a\widetilde{b}`` contributes.


## Reflections and rotations

One of the most important applications of the geometric product is to
represent reflections and rotations via simple products with other
elements of the algebra.  For a unit vector ``𝐧`` (``𝐧^2 = 1``), the
sandwich map ``𝐯 \mapsto -𝐧\,𝐯\,𝐧`` reflects ``𝐯`` through the
hyperplane orthogonal to ``𝐧`` — or, we might say, "along" ``𝐧``.
Decompose ``𝐯 = 𝐯_\parallel + 𝐯_\perp``: the parallel component
commutes with ``𝐧`` (parallel vectors commute), giving
``-𝐧\,𝐯_\parallel\,𝐧 = -𝐯_\parallel``; the perpendicular component
anticommutes, giving ``-𝐧\,𝐯_\perp\,𝐧 = +𝐯_\perp``.  So the
component along ``𝐧`` is flipped and the hyperplane component is
preserved — the correct reflection.

The [Cartan–Dieudonné
theorem](https://en.wikipedia.org/wiki/Cartan%E2%80%93Dieudonn%C3%A9_theorem)
states that every rotation of ``\mathbb{R}^n`` is a composition of an
even number of reflections.  Composing two reflections along unit
vectors ``𝐦`` and ``𝐧``:
```math
-𝐦\bigl(-𝐧\,𝐯\,𝐧\bigr)𝐦 = (𝐦𝐧)\,𝐯\,(𝐦𝐧)^{-1}.
```
The product ``R = 𝐦𝐧`` is an even-grade element satisfying ``R
\widetilde{R} = 1``.  For unit vectors separated by angle ``\alpha``
this product equals ``\cos\alpha + \sin\alpha\,\hat{B}``, where
``\hat{B}`` is the unit bivector of the plane spanned by ``𝐦`` and
``𝐧``.  The effect of ``R\,𝐯\,R^{-1}`` is to rotate ``𝐯`` by an
angle ``2\alpha`` in the plane of ``\hat{B}``.  Setting ``\alpha =
\vartheta/2``:
```math
R = \exp\!\left(\frac{\vartheta}{2}\,\hat{B}\right)
  = \cos\frac{\vartheta}{2} + \sin\frac{\vartheta}{2}\,\hat{B},
```
and ``R\,𝐯\,R^{-1} = R\,𝐯\,\widetilde{R}`` rotates the component of
``𝐯`` in the plane of ``\hat{B}`` by angle ``\vartheta``.

The set of all such "rotors" ``R`` forms a group under multiplication,
called the *spin group*: ``\mathrm{Spin}(3) ≅ \mathrm{SU}(2)``.
The factor of ``1/2`` is a hallmark of the well known double-cover
structure over the rotation group ``\mathrm{SO}(3)``.

Each plane is represented by either of two unit bivectors that differ
by a sign.  For example, the ``𝐱𝐲``-plane is represented by either
``𝐱𝐲`` or ``𝐲𝐱 = -𝐱𝐲``.  We must determine which sign gives a
*right-handed* rotation.  Consider ``R = \exp(\frac{\vartheta}{2}\,𝐱𝐲) =
\cos\frac{\vartheta}{2} + \sin\frac{\vartheta}{2}\,(𝐱𝐲)`` acting on ``𝐱``.  We can readily calculate
```math
R\,𝐱\,\widetilde{R}
= \left(\cos^2\!\frac{\vartheta}{2} - \sin^2\!\frac{\vartheta}{2}\right)\,𝐱
  - 2\sin\frac{\vartheta}{2}\cos\frac{\vartheta}{2}\,𝐲
= \cos(\vartheta)\,𝐱 - \sin(\vartheta)\,𝐲.
```
A right-handed rotation by ``\vartheta`` about ``𝐳`` (curling the
right hand from ``𝐱`` toward ``𝐲``, thumb along ``+𝐳``) maps ``𝐱
\to \cos\vartheta\,𝐱 + \sin\vartheta\,𝐲``.  The formula above has a
negative ``\sin`` term: ``\exp(\theta\,𝐱𝐲)`` rotates ``𝐱`` toward
``-𝐲``, i.e., *clockwise* in the ``𝐱𝐲``-plane (left-handed about
``𝐳``).  Right-handed rotation by ``\vartheta`` about ``+𝐳``
requires the opposite sign: ``\exp(\frac{\vartheta}{2}\,𝐲𝐱)``.  The
same test on the remaining coordinate planes — checking whether
``\exp(\vartheta\,𝐲𝐳)`` maps ``𝐲`` toward ``+𝐳`` or ``-𝐳``, and
whether ``\exp(\vartheta\,𝐳𝐱)`` maps ``𝐳`` toward ``+𝐱`` or
``-𝐱`` — gives in each case the same conclusion: the bivector with
*reversed* index order is the right-handed one.  The results are:
```math
\begin{array}{ll}
\text{right-handed rotation about } 𝐱: & \exp\left(\tfrac{\vartheta}{2}\,𝐳𝐲\right), \\[4pt]
\text{right-handed rotation about } 𝐲: & \exp\left(\tfrac{\vartheta}{2}\,𝐱𝐳\right), \\[4pt]
\text{right-handed rotation about } 𝐳: & \exp\left(\tfrac{\vartheta}{2}\,𝐲𝐱\right).
\end{array}
```

## The even subalgebra: Quaternions

Products of an *even* number of vectors — grades 0 and 2 — form a
closed subalgebra, the *even subalgebra*, spanned by ``\{𝟏, 𝐱𝐲,
𝐲𝐳, 𝐳𝐱\}``.[^2]  This four-dimensional algebra is precisely
isomorphic to the algebra of quaternions.

[^2]: Elements built from an *odd* number of vectors do not form a
    closed subalgebra; two such elements will multiply to give a
    multivector with only *even* grades.  Odd-grade elements are
    useful for representing reflections, but are not used frequently.

In quaternion terminology, the bivector basis elements are canonically
presented as the generators of rotations about the specific axes.
Comparing with the results at the end of the previous section, the
right-handed generators are
```math
\begin{aligned}
𝐢 &= -\star 𝐱 = 𝐳𝐲, \\
𝐣 &= -\star 𝐲 = 𝐱𝐳, \\
𝐤 &= -\star 𝐳 = 𝐲𝐱.
\end{aligned}
```
Each of these squares to ``-1`` (as expected for a unit bivector), and
one can verify the standard cyclic rules:
```math
𝐢^2 = 𝐣^2 = 𝐤^2 = 𝐢𝐣𝐤 = -\boldsymbol{1},
\qquad
𝐢𝐣 = 𝐤, \quad 𝐣𝐤 = 𝐢, \quad 𝐤𝐢 = 𝐣.
```
The opposite sign choice, ``𝐢' = \star 𝐱 = 𝐲𝐳``, produces
left-handed rotations and satisfies ``𝐢'𝐣'𝐤' = +\boldsymbol{1}``.
This is the convention used in much of the engineering literature.[^3]

[^3]: The two quaternion conventions are sometimes called the
      *Hamilton* convention (``𝐢𝐣𝐤 = -1``, used here and in most of
      mathematics and physics) and the *Shuster* or *JPL* convention
      (``𝐢'𝐣'𝐤' = +1``, common in engineering).  However, note that
      this difference may also be compensated by a different choice of
      the "sandwich" product: ``𝐯' = \widetilde{𝐑}\, 𝐯\, 𝐑``
      instead of ``𝐯' = 𝐑 𝐯\, \widetilde{𝐑}`` as we use.  This may
      also be degenerate with different choices of active vs. passive
      rotations.  See [SommerEtAl_2018](@citet) for a thorough
      discussion of the differences and how they arise.

A general quaternion can be written
```math
𝐐 = w\,\boldsymbol{1} + x\,𝐢 + y\,𝐣 + z\,𝐤,
```
and is represented in this package by the `Quaternion` type.
Specifically, ``(w, x, y, z)`` are the names of the components used
throughout this package, and the symbols `𝐢`, `𝐣`, and `𝐤` are
exported with precisely the meaning given here:
```jldoctest
julia> using Quaternionic

julia> 𝐢^2 == 𝐣^2 == 𝐤^2 == 𝐢*𝐣*𝐤 == -1
true

julia> 𝐢*𝐣 == 𝐤 && 𝐣*𝐤 == 𝐢 && 𝐤*𝐢 == 𝐣
true
```

For this general quaternion, the reverse is
```math
\widetilde{𝐐} = w - x\,𝐢 - y\,𝐣 - z\,𝐤,
```
which is exactly the *quaternion conjugate* `conj(Q)` in the code.
Two important consequences follow.  First,
```math
𝐐\,\widetilde{𝐐} = w^2 + x^2 + y^2 + z^2
```
is a *pure scalar* — a non-negative real number for real quaternions.
In analogy with `abs2` for complex numbers, this product is `abs2(Q)`
in the code, and `abs(Q)` is its square root.  Second, every nonzero
quaternion has an inverse
```math
𝐐^{-1} = \frac{\widetilde{𝐐}}{𝐐\,\widetilde{𝐐}}
= \frac{\texttt{conj(Q)}}{\texttt{abs2(Q)}}
= \texttt{inv(Q)}.
```

!!! note "Real vs. complex quaternions"

    For real quaternions, ``𝐐\,\widetilde{𝐐}`` is a positive real scalar.
    For complex quaternions, ``𝐐\,\widetilde{𝐐}`` is still a *scalar*
    element of the algebra (no ``𝐢``, ``𝐣``, or ``𝐤`` component), but it
    is complex-valued.  In particular, this is the *spinor norm*
    ``\sum_i z_i^2`` rather than the Euclidean norm ``\sum_i |z_i|^2``.
    The distinction matters for Lorentz boosts, as discussed on the
    [spacetime algebra page](spacetime_algebra.md).

The general quaternion is represented in this package by the
`Quaternion` type.  Unit quaternions (rotors, which have
``𝐐\,\widetilde{𝐐} = 1``) are represented by the `Rotor` type.  And
"pure-vector" quaternions (having ``w=0``) are represented by the
`QuatVec` type.

Note that the function-call syntax `R(v)` implements the sandwich
``𝐑\, 𝐯\, \widetilde{𝐑}`` efficiently with a dedicated formula that
is roughly twice as fast as the multiplication sequence `R * v *
conj(R)`:
```julia
R * v * conj(R) ≈ R(v)   # both correct; R(v) is ~2× faster
```
This applies to both `Quaternion` and `Rotor` objects acting on a
`QuatVec`.


## The vector–quaternion isomorphism and rotations

We have already demonstrated an isomorphism between the vector space
of ``ℝ³`` and the space of "pure-vector" quaternions via the
(negative) Hodge dual:
```math
𝐯 \mapsto 𝐕 = -\star 𝐯.
```
Here, ``𝐕`` represents the corresponding quaternion.  The inverse
map is given by the same formula:[^4]
```math
𝐕 \mapsto 𝐯 = -\star 𝐕.
```
This isomorphism is a vector-space isomorphism.  But, importantly, it
is also compatible with rotation in the sense that rotating ``𝐯`` is
the same as mapping to ``𝐕``, rotating ``𝐕``, and then mapping back
to the original vector space.  That is,
```math
R 𝐯 \widetilde{R} = -\star \left( R (-\star 𝐯) \widetilde{R} \right).
```

[^4]: This is true for vectors in ``ℝ³``.  There may be an additional
    sign for different grades or spaces: for a grade-``r`` element in
    an ``n``-dimensional space with signature ``(p,q)``, we have
    ``\star \star b = (-1)^{r(n-r)} \chi(q) b``.  Here, ``\chi`` is
    the parity function which is ``+1`` for even ``q`` and ``-1`` for
    odd ``q``.

Thus, we can freely move back and forth between the two
representations of vectors.


## Further reading

The geometric algebra perspective on quaternions and rotors is developed in
depth in [DoranLasenby_2010](@citet).

```@bibliography
Pages = ["geometric_algebra.md"]
Canonical = false
```
