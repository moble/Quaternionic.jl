# Quaternions and the Spacetime Algebra

The [geometric algebra page](geometric_algebra.md) showed that *real*
quaternions arise as the even subalgebra of the geometric algebra over ℝ³.
This page extends that picture to four-dimensional Minkowski spacetime,
where the even subalgebra is isomorphic to the *complexified* quaternions
ℍ(ℂ), and Lorentz boosts appear alongside rotations as unit elements of
that algebra.

The central new feature is the *spinor norm*: for complex quaternions,
the GA reverse gives `𝐐𝐐̃ = w² + x² + y² + z²` (a complex scalar), not
the Euclidean `|w|² + |x|² + |y|² + |z|²`.  Both rotation and boost
rotors are unit elements under the spinor norm — but not under the Euclidean
norm.

The basic rotation machinery (reflections, the Cartan–Dieudonné theorem,
the sandwich product) carries over unchanged from the GA page; we focus here
on what is genuinely new.


## The spacetime algebra

The *spacetime algebra* (STA) is the geometric algebra over four-dimensional
Minkowski space with metric signature ``(-,+,+,+)``.  The basis vectors
``𝐭, 𝐱, 𝐲, 𝐳`` satisfy
```math
𝐭^2 = -1, \qquad 𝐱^2 = 𝐲^2 = 𝐳^2 = +1,
```
and all off-diagonal products anticommute.  The full algebra has dimension
``2^4 = 16``, with basis elements grouped by grade:

```math
\begin{array}{lll}
\text{grade 0:} & 1 \text{ element} & \boldsymbol{1} \\[4pt]
\text{grade 1:} & 4 \text{ elements} & 𝐭,\; 𝐱,\; 𝐲,\; 𝐳 \\[4pt]
\text{grade 2:} & 6 \text{ elements} & 𝐭𝐱,\; 𝐭𝐲,\; 𝐭𝐳,\quad 𝐱𝐲,\; 𝐱𝐳,\; 𝐲𝐳 \\[4pt]
\text{grade 3:} & 4 \text{ elements} & 𝐭𝐱𝐲,\; 𝐭𝐱𝐳,\; 𝐭𝐲𝐳,\; 𝐱𝐲𝐳 \\[4pt]
\text{grade 4:} & 1 \text{ element} & 𝐈 = 𝐭𝐱𝐲𝐳
\end{array}
```

The grade-2 elements split into two qualitatively different groups.  The
three *spatial bivectors* ``𝐱𝐲``, ``𝐱𝐳``, ``𝐲𝐳`` involve only spatial
directions and square to ``-1``, exactly as in ℝ³.  The three *boost
bivectors* ``𝐭𝐱``, ``𝐭𝐲``, ``𝐭𝐳`` each involve the timelike direction,
and they square to ``+1`` because of the sign in ``𝐭^2 = -1``:
```math
(𝐭𝐱)^2 = 𝐭𝐱𝐭𝐱 = -𝐭𝐭𝐱𝐱 = -(𝐭^2)(𝐱^2) = -(-1)(+1) = +1.
```
This sign difference is what distinguishes boosts from rotations.


## The pseudoscalar in four dimensions

The pseudoscalar ``𝐈 = 𝐭𝐱𝐲𝐳`` has several properties that differ
markedly from its three-dimensional counterpart:

**It still squares to ``-1``.**  In an ``n``-dimensional space with metric
``g``, the pseudoscalar satisfies ``𝐈^2 = (-1)^{n(n-1)/2} \det(g)``.  For
the STA: ``n = 4`` and ``\det(g) = (-1)(+1)^3 = -1``, giving
``(-1)^6 \times (-1) = -1``. ✓

**Its reverse is itself.**  The reverse of a grade-``r`` blade picks up a
sign ``(-1)^{r(r-1)/2}``.  For ``r = 4``: ``(-1)^{4 \cdot 3/2} = (-1)^6 = +1``,
so
```math
\widetilde{𝐈} = +𝐈.
```
Contrast the three-dimensional pseudoscalar (grade 3), for which
``(-1)^{3 \cdot 2/2} = -1``, giving ``\widetilde{𝐈}_{3\mathrm{D}} = -𝐈_{3\mathrm{D}}``.

**It does not commute with everything.**  Moving a grade-``n`` pseudoscalar
past a grade-``s`` element costs ``(-1)^{s(n-1)}`` sign changes.  In three
dimensions ``(-1)^{2s} = +1`` always, so ``𝐈_{3\mathrm{D}}`` commutes with
everything.  In four dimensions ``(-1)^{3s}``:
```math
\begin{array}{ll}
\text{grade 0 (scalars):} & (-1)^0 = +1 \quad \text{(commutes)} \\[4pt]
\text{grade 1 (vectors):} & (-1)^3 = -1 \quad \text{(anti-commutes)} \\[4pt]
\text{grade 2 (bivectors):} & (-1)^6 = +1 \quad \text{(commutes)} \\[4pt]
\text{grade 3 (trivectors):} & (-1)^9 = -1 \quad \text{(anti-commutes)}
\end{array}
```
In particular, ``𝐈`` commutes with every element of the even subalgebra
(grades 0 and 2), which is the property we will need.


## The even subalgebra ≅ ℍ(ℂ)

Products of an even number of STA basis vectors span a closed subalgebra of
dimension ``1 + 6 + 1 = 8`` (grades 0, 2, and 4).  Eight real dimensions is
exactly the dimension of the complexified quaternions ``\mathbb{H}(\mathbb{C})``,
and the identification is explicit.

**Spatial bivectors → quaternion units.**  The three spatial bivectors that
generated right-handed rotations in ℝ³ are unchanged:
```math
𝐢 = 𝐳𝐲, \qquad 𝐣 = 𝐱𝐳, \qquad 𝐤 = 𝐲𝐱.
```
They still square to ``-1`` and satisfy ``𝐢𝐣𝐤 = -\boldsymbol{1}``.

**Pseudoscalar → imaginary unit.**  Since ``𝐈^2 = -1`` and ``𝐈`` commutes
with all bivectors and scalars, it plays the role of the imaginary unit
``\mathbf{i} = \sqrt{-1}`` on the even subalgebra.  We write ``𝒾`` for this
imaginary unit when we want to emphasize its algebraic role:
```math
𝒾 \;\leftrightarrow\; 𝐈 = 𝐭𝐱𝐲𝐳.
```

**Boost bivectors → ``𝒾 \times`` quaternion units.**  The three boost
bivectors are products of ``𝐈`` with the spatial bivectors.  For example:
```math
𝐈\,𝐢 = (𝐭𝐱𝐲𝐳)(𝐳𝐲)
= 𝐭𝐱𝐲(𝐳𝐳)𝐲
= 𝐭𝐱𝐲𝐲
= 𝐭𝐱(𝐲^2)
= 𝐭𝐱,
```
and similarly ``𝐈\,𝐣 = 𝐭𝐲`` and ``𝐈\,𝐤 = 𝐭𝐳``.  In the quaternion
picture these are ``𝒾𝐢``, ``𝒾𝐣``, ``𝒾𝐤``.  Consistent with
``(𝐭𝐱)^2 = +1``, we have ``(𝒾𝐢)^2 = 𝒾^2 𝐢^2 = (-1)(-1) = +1``. ✓

A general element of the even subalgebra is therefore a quaternion
``𝐐 = w\,\boldsymbol{1} + x\,𝐢 + y\,𝐣 + z\,𝐤`` with *complex*
coefficients ``(w, x, y, z) \in \mathbb{C}``, represented by the
`Quaternion{Complex{T}}` type in this package.


## The reverse and the spinor norm

For a general complex quaternion ``𝐐 = w + x\,𝐢 + y\,𝐣 + z\,𝐤``, the
reverse is computed grade by grade.  Grades 0 and 2 behave the same as for
real quaternions (grades 0 unchanged, grade 2 negated), and the
pseudoscalar grade-4 part (if present) is unchanged by the reverse
(``\widetilde{𝐈} = +𝐈``).  For an element of the even subalgebra this gives
```math
\widetilde{𝐐} = w - x\,𝐢 - y\,𝐣 - z\,𝐤 = \texttt{conj}(Q),
```
exactly as for real quaternions.

The *spinor norm* is
```math
𝐐\,\widetilde{𝐐} = w^2 + x^2 + y^2 + z^2,
```
a complex scalar.  For real quaternions ``z_i \in \mathbb{R}``, this equals
the Euclidean norm ``\sum |z_i|^2``.  For complex quaternions the two differ:
the spinor norm is complex-valued in general, while the Euclidean norm is
always real and non-negative.

The code functions `abs2(Q)` and `abs(Q)` compute the spinor norm and its
(complex) square root.  This is why `abs2` for a boost rotor is not `1` in
the ordinary sense — it equals `1` as a complex number, but the *modulus*
``|``abs2``(Q)|`` can be greater than 1.

!!! note "Real vs. complex quaternions"

    For real quaternions, ``𝐐\,\widetilde{𝐐} = \sum r_i^2 = \sum |r_i|^2``,
    so the spinor and Euclidean norms coincide.  This is why the distinction
    never arises on the GA page.  For complex quaternions they diverge, and
    only the spinor norm has the correct geometric meaning for Lorentz
    transformations.


## Rotors in the STA

### Rotation rotors

A rotation rotor is an element of the even subalgebra with real components
and spinor norm 1.  It has the same form as on the GA page:
```math
R = \exp\!\left(\frac{\vartheta}{2}\,\hat{B}\right)
  = \cos\frac{\vartheta}{2} + \sin\frac{\vartheta}{2}\,\hat{B},
```
where ``\hat{B} \in \{𝐳𝐲, 𝐱𝐳, 𝐲𝐱\}`` is a unit spatial bivector.
The components are real, and the spinor norm equals the Euclidean norm:
``R\widetilde{R} = \cos^2(\vartheta/2) + \sin^2(\vartheta/2) = 1``. ✓

### Boost rotors

A boost in the ``𝐱``-direction by rapidity ``\varphi`` uses the boost
bivector ``𝐭𝐱 = 𝒾\,𝐢``:
```math
B = \exp\!\left(-\frac{\varphi}{2}\,𝐭𝐱\right)
  = \exp\!\left(-\frac{𝒾\varphi}{2}\,𝐢\right)
  = \cosh\frac{\varphi}{2} - 𝒾\sinh\frac{\varphi}{2}\,𝐢.
```
The sign convention ``-\varphi/2`` makes a positive rapidity correspond to
a boost in the positive ``x``-direction.[^boost]

The components of ``B`` are ``(w, x, y, z) = (\cosh(\varphi/2),\; -i\sinh(\varphi/2),\; 0,\; 0)``
(where ``i = \sqrt{-1}`` is the ordinary complex unit).  The spinor norm is
```math
B\,\widetilde{B}
= \cosh^2\!\frac{\varphi}{2} + \left(-i\sinh\frac{\varphi}{2}\right)^{\!2}
= \cosh^2\!\frac{\varphi}{2} - \sinh^2\!\frac{\varphi}{2}
= 1, \qquad \checkmark
```
while the Euclidean norm is ``\cosh^2(\varphi/2) + \sinh^2(\varphi/2) = \cosh\varphi \neq 1``
for ``\varphi \neq 0``.  This concretely shows why the spinor norm is the
physically correct notion of "unit" for Lorentz transformations.

[^boost]: The sign can be verified by expanding the Lorentz transformation
    ``\Lambda\,𝐯\,\widetilde{\Lambda}`` on the 4-vector ``𝐯 = t\,𝐭 + x\,𝐱``
    and confirming that a positive rapidity boosts ``t \to \cosh\varphi\,t +
    \sinh\varphi\,x`` and ``x \to \sinh\varphi\,t + \cosh\varphi\,x``.

### Phase factors are not Lorentz rotors

A pure phase ``\exp(𝒾\varphi) = \cos\varphi + 𝒾\sin\varphi`` is a complex
scalar (grade 0 of the even subalgebra).  Its reverse equals itself (grade 0
is unchanged), so its spinor norm is
```math
\exp(𝒾\varphi)\,\widetilde{\exp(𝒾\varphi)} = \exp(𝒾\varphi)^2 = \exp(2𝒾\varphi),
```
which is not equal to ``1`` in general.  Phase factors therefore do *not*
belong to the Lorentz rotor group, even though they are unit elements under
the Euclidean norm (``|\exp(𝒾\varphi)| = 1``).


## Further reading

The spacetime algebra and its application to relativistic physics are
developed in detail in [DoranLasenby_2010](@citet).

```@bibliography
Pages = ["spacetime_algebra.md"]
Canonical = false
```
