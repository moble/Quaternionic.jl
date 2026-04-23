# Quaternions and the Spacetime Algebra

The [geometric algebra page](geometric_algebra.md) showed that *real*
quaternions ``ℍ`` arise as the even subalgebra of the geometric
algebra over ``ℝ³``.  This page extends that picture to
four-dimensional Minkowski spacetime ``ℝ^{3,1}``, where we simply
include a fourth basis vector ``𝐭`` which squares to ``-1``.  The
geometric algebra over this space is a simple extension of the
three-dimensional case, producing the "Spacetime Algebra" (STA).

*By sheer coincidence* the even subalgebra happens to be isomorphic to
the *complexified* quaternions ``ℍ(ℂ)`` — sometimes also called
"biquaternions".  The isomorphism also happens to be very simple: the
unit imaginary scalar in the complexified quaternions corresponds to
the pseudoscalar ``𝐈 = 𝐭𝐱𝐲𝐳`` in the spacetime algebra, and the
complex components of the bivectors represent timelike bivectors,
which correspond to ``𝐈`` times spatial bivectors.  That is, we just
think of the unit imaginary ``i`` as being equivalent to the
pseudoscalar ``𝐈``.

The basic rotation machinery (reflections, the Cartan–Dieudonné
theorem, the sandwich product) carries over unchanged from the
three-dimensional case.  And just as rotations are represented by unit
quaternions in the three-dimensional case, general Lorentz
transformations are represented by unit complexified quaternions in
this more general algebra.

One crucial new feature is the *spinor norm*: for complex quaternions,
the GA reverse gives ``𝐐𝐐̃ = w² + x² + y² + z²``, which is a
*"complex"* scalar, not the Euclidean norm ``|w|² + |x|² + |y|² +
|z|²`` which will always be real.  The surprising point is that the
normalization condition ``𝐐𝐐̃ = 1`` is now *two real* conditions,
because this says that the "imaginary" part is zero.  This reduces the
eight real degrees of freedom in a complex quaternion down to six,
which is the correct number for a Lorentz transformation in four
dimensions.


## The spacetime algebra

The *spacetime algebra* (STA) is the geometric algebra over
four-dimensional Minkowski space, and we use the metric signature
``{-}{+}{+}{+}``.  The basis vectors ``𝐭, 𝐱, 𝐲, 𝐳`` satisfy
```math
𝐭^2 = -1, \qquad 𝐱^2 = 𝐲^2 = 𝐳^2 = +1,
```
and all other products anticommute.  The full algebra has dimension
``2^4 = 16``, with basis elements grouped by grade:

```math
\begin{array}{lll}
\text{grade 0:} & 1 \text{ element} & \boldsymbol{1} \\[4pt]
\text{grade 1:} & 4 \text{ elements} & 𝐭,\; 𝐱,\; 𝐲,\; 𝐳 \\[4pt]
\text{grade 2:} & 6 \text{ elements} & 𝐭𝐱,\; 𝐭𝐲,\; 𝐭𝐳,\; 𝐱𝐲,\; 𝐱𝐳,\; 𝐲𝐳 \\[4pt]
\text{grade 3:} & 4 \text{ elements} & 𝐭𝐱𝐲,\; 𝐭𝐱𝐳,\; 𝐭𝐲𝐳,\; 𝐱𝐲𝐳 \\[4pt]
\text{grade 4:} & 1 \text{ element} & 𝐈 = 𝐭𝐱𝐲𝐳
\end{array}
```

The grade-2 elements split into two qualitatively different groups.
The three *spatial bivectors* ``𝐱𝐲``, ``𝐱𝐳``, ``𝐲𝐳`` involve
only spatial directions and square to ``-1``, exactly as in ``ℝ³``.
The three *boost bivectors* ``𝐭𝐱``, ``𝐭𝐲``, ``𝐭𝐳`` each involve
the timelike direction, and they square to ``+1`` because of the sign
in ``𝐭^2 = -1``:
```math
(𝐭𝐱)^2 = 𝐭𝐱𝐭𝐱 = -𝐭𝐭𝐱𝐱 = -(𝐭^2)(𝐱^2) = -(-1)(+1) = +1.
```
This sign difference is what distinguishes boosts from rotations.


## The pseudoscalar in four dimensions

The pseudoscalar ``𝐈 = 𝐭𝐱𝐲𝐳`` has several properties that differ
markedly from its three-dimensional counterpart:

**It still squares to ``-1``.**  Though the inclusion of an extra
basis factor gives the permutation a different sign, the fact that
``𝐭² = -1`` means that the overall result is unchanged.

**Its reverse is itself.**  The reverse of a grade-``r`` blade picks
up a sign ``(-1)^{r(r-1)/2}``.  For ``r = 4``, this is ``(-1)^6 =
+1``, so
```math
\widetilde{𝐈} = +𝐈.
```
Contrast the (grade 3) three-dimensional pseudoscalar ``𝐈_{3}``, for
which ``\widetilde{𝐈}_{3} = -𝐈_{3}``.

**It does not commute with everything.**  Whereas ``𝐈_{3}`` commutes
with everything in its three-dimensional algebra, ``𝐈`` does not
commute with everything in the four-dimensional algebra.  In
particular, ``𝐈`` *anti*-commutes with vectors and trivectors.
However, it does commute with every element of the even subalgebra
(grades 0, 2, and 4), which is an important property.


## The even subalgebra: Complexified quaternions

Products of an *even* number of vectors — grades 0, 2, and 4 — form a
closed subalgebra, the *even subalgebra*, spanned by ``\{1, 𝐱𝐲,
𝐲𝐳, 𝐳𝐱, 𝐭𝐱, 𝐭𝐲, 𝐭𝐳, 𝐈\}``.  This algebra happens to be
precisely isomorphic to the one spanned by ``\{1, 𝐱𝐲, 𝐲𝐳, 𝐳𝐱,
i𝐱𝐲, i𝐲𝐳, i𝐳𝐱, i\}``, where ``i`` is the unit imaginary — which
is of course precisely the algebra of *complexified* quaternions
``ℍ(ℂ)``.

**Pseudoscalar → imaginary unit.**  Since ``𝐈^2 = -1`` and ``𝐈``
commutes with all scalars, bivectors, and the pseudoscalar, it plays
the role of the imaginary unit ``i = \sqrt{-1}`` on the even
subalgebra.  Thus, we have this basic identification:
```math
i \;\leftrightarrow\; 𝐈 = 𝐭𝐱𝐲𝐳.
```
We will abuse notation and terminology to write complex numbers
interchangeably as both ``a + i b`` and ``a + 𝐈 b``, and speak of
linear combinations of ``\{1,i\}`` and ``\{1,𝐈\}`` as complex
numbers.

**Spatial bivectors → quaternion units.**  The three spatial bivectors
that generated right-handed rotations in ``ℝ³`` are unchanged:
```math
𝐢 = 𝐳𝐲, \qquad 𝐣 = 𝐱𝐳, \qquad 𝐤 = 𝐲𝐱.
```
They still square to ``-1`` and satisfy ``𝐢𝐣𝐤 = -\boldsymbol{1}``.

**Boost bivectors → ``i \times`` quaternion units.**  The three boost
bivectors are products of ``𝐈`` with the spatial bivectors.  For
example:
```math
𝐈\,𝐢 = (𝐭𝐱𝐲𝐳)(𝐳𝐲)
= 𝐭𝐱𝐲(𝐳𝐳)𝐲
= 𝐭𝐱𝐲𝐲
= 𝐭𝐱(𝐲^2)
= 𝐭𝐱,
```
and similarly ``𝐈\,𝐣 = 𝐭𝐲`` and ``𝐈\,𝐤 = 𝐭𝐳``.  In the
complexified quaternion picture these are ``i𝐢``, ``i𝐣``, ``i𝐤``.
Consistent with ``(𝐭𝐱)^2 = +1``, we have ``(i𝐢)^2 = i^2 𝐢^2 =
(-1)(-1) = +1``, and so on.

A general element of the even subalgebra is therefore a quaternion
```math
\begin{aligned}
𝐐 &= w\,\boldsymbol{1} + x\,𝐢 + y\,𝐣 + z\,𝐤 \\
&= \Re\{w\}\,\boldsymbol{1} + \Re\{x\}\,𝐳𝐲 + \Re\{y\}\,𝐱𝐳 + \Re\{z\}\,𝐲𝐱
+ \Im\{x\}\,𝐭𝐱 + \Im\{y\}\,𝐭𝐲 + \Im\{z\}\,𝐭𝐳 + \Im\{w\}\,𝐈
\end{aligned}
```
with *complex* coefficients ``(w, x, y, z) \in \mathbb{C}``,
represented by the `Quaternion{Complex{T}}` type in this package.


## The reverse and the spinor norm

As usual in Geometric Algebra, the *reverse* of a multivector is
computed by reversing the order of all products.  Obviously, the
reverse of a scalar is itself, and — as mentioned above — the reverse
of the pseudoscalar is itself.  The reverse applied to a bivector is
just negation.  Therefore, the reverse applied to a general complex
quaternion is exactly the same as for a real quaternion: the scalar
part is unchanged, and the vector part is negated.  The complex
components themselves are not conjugated:
```math
\widetilde{𝐐} = w - x\,𝐢 - y\,𝐣 - z\,𝐤 = \texttt{conj(Q)},
```
exactly as for real quaternions.

Now, we have to be careful to distinguish between two different
notions of "norm" for complex quaternions: the *spinor norm* and the
*Euclidean norm*.  The spinor norm,
```math
𝐐\,\widetilde{𝐐} = w^2 + x^2 + y^2 + z^2,
```
is the one that arises from the GA reverse, and it is the physically
meaningful notion of "unit" for Lorentz transformations.  Notably, the
result is a *complex* number.  On the other hand, the Euclidean norm,
```math
𝐐\,𝐐^\dagger = |w|^2 + |x|^2 + |y|^2 + |z|^2,
```
is always a positive real number, and does not have a direct geometric
meaning in this context.  The relevant condition that defines a spinor
(in this case, a Lorentz rotor) is that the *spinor* norm equals 1:
```math
𝐐\,\widetilde{𝐐} = 1 + 0i.
```
The code functions `abs2(Q)` and `abs(Q)` compute the *spinor* norm
and its (complex) square root.

One possibly surprising consequence of this is that a pure phase
``\exp(𝐈\alpha)`` is not generally a Lorentz rotor.  Its spinor norm
is
```math
\exp(𝐈\alpha)\,\widetilde{\exp(𝐈\alpha)} = \exp(𝐈\alpha)^2 = \exp(2𝐈\alpha),
```
which is not equal to ``1`` in general.  Phase factors therefore do *not*
belong to the Lorentz rotor group, even though they are unit elements under
the Euclidean norm (``|\exp(i\varphi)| = 1``).



!!! note "Real vs. complex quaternions"

    For real quaternions, ``𝐐\,\widetilde{𝐐} = \sum r_i^2 = \sum |r_i|^2``,
    so the spinor and Euclidean norms coincide.  This is why the distinction
    never arises on the GA page.  For complex quaternions they diverge, and
    only the spinor norm has the correct geometric meaning for Lorentz
    transformations.


## Rotors in the STA

A rotation rotor is an element of the even subalgebra with purely real
components and spinor norm 1.  It has the same form as in the
three-dimensional case:
```math
R = \exp\!\left(\frac{\vartheta}{2}\,\hat{B}\right)
  = \cos\frac{\vartheta}{2} + \sin\frac{\vartheta}{2}\,\hat{B},
```
where ``\hat{B} \in \{𝐳𝐲, 𝐱𝐳, 𝐲𝐱\}`` is a unit spatial bivector.
The components are real, and the spinor norm happens to equal the
Euclidean norm: ``R\widetilde{R} = \cos^2(\vartheta/2) +
\sin^2(\vartheta/2) = 1``.

A boost in the ``𝐯``-direction by rapidity ``\varphi`` uses the boost
bivector ``𝐭𝐯``:
```math
B = \exp\!\left(-\frac{\varphi}{2}\,𝐭𝐯\right)
  = \cosh\frac{\varphi}{2} - \sinh\frac{\varphi}{2}\,𝐭𝐯.
```
The sign convention ``-\varphi/2`` makes a positive rapidity
correspond to a boost in the positive ``𝐯``-direction.[^boost]  The
spinor norm is
```math
B\,\widetilde{B}
= \cosh^2\!\frac{\varphi}{2} + \left(-𝐭𝐯\sinh\frac{\varphi}{2}\right)^{\!2}
= \cosh^2\!\frac{\varphi}{2} - \sinh^2\!\frac{\varphi}{2}
= 1,
```
while the Euclidean norm is ``\cosh^2(\varphi/2) + \sinh^2(\varphi/2)
= \cosh\varphi \neq 1`` for ``\varphi \neq 0``.  This concretely shows
why the spinor norm is the physically correct notion of "unit" for
Lorentz transformations.

[^boost]: The sign can be verified by expanding the Lorentz transformation
    ``\Lambda\,𝐯\,\widetilde{\Lambda}`` on the 4-vector ``𝐯 = t\,𝐭 + x\,𝐱``
    and confirming that a positive rapidity boosts ``t \to \cosh\varphi\,t +
    \sinh\varphi\,x`` and ``x \to \sinh\varphi\,t + \cosh\varphi\,x``.

## Further reading

The spacetime algebra and its application to relativistic physics are
developed in detail in [DoranLasenby_2010](@citet).

```@bibliography
Pages = ["spacetime_algebra.md"]
Canonical = false
```
