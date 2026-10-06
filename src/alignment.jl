@doc raw"""
    align(a⃗::AbstractArray{<:QuatVec}, b⃗::AbstractArray{<:QuatVec}, [w])

Solve [Wahba's problem](https://en.wikipedia.org/wiki/Wahba%27s_problem), finding the rotor
`R` that best rotates the set of points `b⃗` onto a corresponding set of points `a⃗`, so
that `R * b⃗[i] * R⁻¹ ≈ a⃗[i]`.

Here, `a⃗` and `b⃗` must be equally sized arrays of `QuatVec`s with real components.  If
present, `w` must be an equally sized array of real numbers; if not, it is taken to be 1.
We define the loss function
```math
L(ℛ) ≔ Σᵢ wᵢ ‖a⃗ᵢ - ℛ b⃗ᵢ‖²
```
where ``ℛ`` is a rotation operator, and return the `Rotor` corresponding to the optimal
``ℛ`` that minimizes this function.

Since `R` and `-R` represent the same rotation, the sign of the result is fixed by
convention: its scalar part is positive, or, if the scalar part is zero, its first nonzero
component is positive.

Note that it is possible that the points do not uniquely determine a rotation — as when one
or both sets of points is rotationally symmetric.  In that case, the loss function ``L(ℛ)``
will still be minimized and the points will still be optimally aligned by the output rotor,
but that rotor will not be unique.


# Notes

[MarkleyCrassidis_2014](@citet) say that "Davenport’s method remains the best method for
solving Wahba’s problem".  This method provides the optimal quaternion as the dominant
eigenvector (the one with the largest eigenvalue) of a certain matrix.  We start by defining
the supplementary matrix
```math
S ≔ Σᵢ wᵢ a⃗ᵢ b⃗ᵢᵀ
```
and vector
```math
s⃗ ≔ \begin{bmatrix}
S₂₃-S₃₂ \\
S₃₁-S₁₃ \\
S₁₂-S₂₁
\end{bmatrix}.
```
Then the key matrix is
```math
M ≔ \begin{bmatrix}
\mathrm{tr}S & -s⃗ᵀ \\
-s⃗ & S + Sᵀ - (\mathrm{tr}S)\, I₃
\end{bmatrix}.
```
If ``(w, x, y, z)`` is the dominant eigenvector of this matrix, the result is the rotor ``w
+ x𝐢 + y𝐣 + z𝐤``.  (Markley and Crassidis order the components with the scalar last and
use the opposite sign for ``s⃗``, because their convention for quaternion multiplication
differs from the one used here.)  It is possible for this matrix to have degenerate
eigenvalues, corresponding to cases where the points do not uniquely determine the rotation,
as described above.

This method works with any real float type, and the result has the same element type as the
input.  It computes the eigen-decomposition of `M`, which LAPACK provides only for
`Float16`, `Float32`, and `Float64`; for anything else we rely on `GenericLinearAlgebra`, as
[`from_rotation_matrix`](@ref) does, and which `Quaternionic` depends on directly so that
you do not have to load it yourself.

# References

* [MarkleyCrassidis_2014](@cite) F. L. Markley and J. L. Crassidis, _Fundamentals of
  Spacecraft Attitude Determination and Control_ (Springer, New York, 2014)

"""
function align(a⃗::AbstractArray{<:QuatVec}, b⃗::AbstractArray{<:QuatVec}, w::AbstractArray{<:Real})
    # This is Eq. (5.11) from Markley and Crassidis.  Each term of the sum is built from
    # static three-vectors, so that it stays on the stack.  The parentheses prevent the
    # product from being parsed as a three-argument `*`, which would allocate.
    S = sum(
        (w[i] * SVector(a⃗[i][2], a⃗[i][3], a⃗[i][4])) *
            transpose(SVector(b⃗[i][2], b⃗[i][3], b⃗[i][4]))
        for i in eachindex(a⃗, b⃗, w)
    )
    return _align_Wahba(S)
end

function align(a⃗::AbstractArray{<:QuatVec}, b⃗::AbstractArray{<:QuatVec})
    # This is Eq. (5.11) from Markley and Crassidis
    S = sum(
        SVector(a⃗[i][2], a⃗[i][3], a⃗[i][4]) *
            transpose(SVector(b⃗[i][2], b⃗[i][3], b⃗[i][4]))
        for i in eachindex(a⃗, b⃗)
    )
    return _align_Wahba(S)
end

function _align_Wahba(S)
    # This is Eq. (5.17) from Markley and Crassidis, modified to suit our conventions by
    # flipping the sign of their vector `z` (the `s⃗` of the docstring), and moving the
    # final dimension to the first dimension.  The elements are given in column-major
    # order.
    M = Symmetric(SMatrix{4,4}(
        S[1,1]+S[2,2]+S[3,3], S[3,2]-S[2,3], S[1,3]-S[3,1], S[2,1]-S[1,2],
        S[3,2]-S[2,3], S[1,1]-S[2,2]-S[3,3], S[1,2]+S[2,1], S[1,3]+S[3,1],
        S[1,3]-S[3,1], S[1,2]+S[2,1], -S[1,1]+S[2,2]-S[3,3], S[2,3]+S[3,2],
        S[2,1]-S[1,2], S[1,3]+S[3,1], S[2,3]+S[3,2], -S[1,1]-S[2,2]+S[3,3]
    ))
    # This extracts the dominant eigenvector, and interprets it as a Rotor.  Note that
    # `dominant_eigenvector` selects by eigen*value* rather than by position, because
    # `eigen` guarantees no particular ordering across the various backends, and it also
    # handles element types that LAPACK cannot.
    return positive_hemisphere(rotor(dominant_eigenvector(M)...))
end


@doc raw"""
    align(A::AbstractArray{<:Rotor}, B::AbstractArray{<:Rotor}, [w])

Find the `Rotor` `R` that best rotates the set of rotors `B` onto a
corresponding set `A`, so that `R * B[i] ≈ A[i]`, by minimizing the distance
between the first set and the rotated second set.

Here, `A` and `B` must be equally sized arrays of `Rotor`s.  If present, `w`
must be an equally sized array of real numbers; if not, it is taken to be 1.  We
define the loss function
```math
L(R) ≔ Σᵢ wᵢ |Aᵢ - R Bᵢ|²
```
where ``R`` is a `Rotor`, and return the optimal ``R`` that minimizes this
function.

Note that it is possible that the input data do not uniquely determine a rotor,
which will happen when the sum below is zero.  When this happens, the result
will contain `NaN`s, but no error will be raised.  When the sum is very close to
— but not exactly — zero, the accuracy of the result will be limited.  However,
the loss function will not depend strongly on the result in that case.

Be aware that this function _is_ sensitive to the signs of the input
quaternions.  See the [`unflip`](@ref) function for one way to avoid problems
related to signs.


## Notes

We can ensure that the loss function is minimized by multiplying ``R`` by an
exponential, differentiating with respect to the argument of the exponential,
and setting that argument to 0.  This derivative should be 0 at the minimum.  We
have
```math
∂ⱼ Σᵢ wᵢ |Aᵢ - \exp[vⱼ] R Bᵢ|²  →  -2 ⟨ eⱼ R Σᵢ wᵢ Bᵢ Āᵢ ⟩₀
```
where → denotes taking ``vⱼ→0``, the symbol ``⟨⟩₀`` denotes taking the scalar
part, and ``eⱼ`` is the unit quaternionic vector in the ``j`` direction.  The
only way for this quantity to be zero for each choice of ``j`` is if
```math
R Σᵢ wᵢ Bᵢ Āᵢ
```
is itself a pure scalar.  This, in turn, can only happen if either (1) the sum
is 0 or (2) if ``R`` is proportional to the _conjugate_ of the sum:
```math
R ∝ Σᵢ wᵢ Aᵢ B̄ᵢ
```
Now, since we want ``R`` to be a rotor, we simply define it to be the normalized
sum.  The positive normalization makes the scalar in question positive, which
gives the minimum of the loss function rather than the maximum.

"""
function align(A::AbstractArray{<:Rotor}, B::AbstractArray{<:Rotor}, w::AbstractArray{<:Real})
    rotor(sum(w[i] * A[i] * conj(B[i]) for i in eachindex(A, B, w)))
end

function align(A::AbstractArray{<:Rotor}, B::AbstractArray{<:Rotor})
    rotor(sum(A[i] * conj(B[i]) for i in eachindex(A, B)))
end
