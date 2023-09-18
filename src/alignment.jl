@doc raw"""
    align(a⃗, b⃗, [w])

Solve [Wahba's problem](https://en.wikipedia.org/wiki/Wahba%27s_problem),
finding a rotation that aligns the set of points `a⃗` to a corresponding set of
points `b⃗` by minimizing the distance between the first set and the rotated
second set.

Here, `a⃗` and `b⃗` must be equally sized arrays of `QuatVec`s.  If present, `w`
must be an equally sized array of real numbers; if not, it is taken to be 1.
We define the loss function
```math
L(ℛ) ≔ Σᵢ wᵢ ‖a⃗ᵢ - ℛ b⃗ᵢ‖²
```
where ``ℛ`` is a rotation operator, and return the quaternion corresponding to
the optimal ``ℛ`` that minimizes this function.

Note that it is possible that the points do not uniquely determine a rotation —
as when one or both sets of points is rotationally symmetric.  In that case, the
loss function ``L(ℛ)`` will still be minimized and the points will still be
optimally aligned by the output quaternion, but that quaternion will not be
unique.


# Notes

In their book [_Fundamentals of Spacecraft Attitude Determination and Control_
(2014)](https://doi.org/10.1007/978-1-4939-0802-8), Markley and Crassidis say
that "Davenport’s method remains the best method for solving Wahba’s problem".
This method provides the optimal quaternion as the dominant eigenvector (the one
with the largest eigenvalue) of a certain matrix.  We start by defining the
supplementary matrix
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
S + Sᵀ - (\mathrm{tr}S)\, I₃ & s⃗ \\
s⃗ᵀ & \mathrm{tr}S
\end{bmatrix}
```
It is possible for this matrix to have degenerate eigenvalues, corresponding to
cases where the points do not uniquely determine the rotation, as described
above.

"""
function align(a⃗::AbstractArray{<:QuatVec}, b⃗::AbstractArray{<:QuatVec}, w::AbstractArray{<:Real})
    # This is Eq. (5.11) from Markley and Crassidis
    S = sum(w[i] * vec(a⃗[i]) * vec(b⃗[i])' for i in eachindex(a⃗, b⃗, w))
    return _align_Wahba(S)
end

function align(a⃗::AbstractArray{<:QuatVec}, b⃗::AbstractArray{<:QuatVec})
    # This is Eq. (5.11) from Markley and Crassidis
    S = sum(vec(a⃗[i]) * vec(b⃗[i])' for i in eachindex(a⃗, b⃗))
    return _align_Wahba(S)
end

function _align_Wahba(S)
    # This is Eq. (5.17) from Markley and Crassidis, modified to suit our
    # conventions by flipping the sign of ``z``, and moving the final dimension
    # to the first dimension.
    M = Symmetric([
            S[1,1]+S[2,2]+S[3,3]      S[3,2]-S[2,3]         S[1,3]-S[3,1]           S[2,1]-S[1,2]
                S[3,2]-S[2,3]      S[1,1]-S[2,2]-S[3,3]     S[1,2]+S[2,1]           S[1,3]+S[3,1]
                S[1,3]-S[3,1]         S[1,2]+S[2,1]      -S[1,1]+S[2,2]-S[3,3]      S[2,3]+S[3,2]
                S[2,1]-S[1,2]         S[1,3]+S[3,1]         S[2,3]+S[3,2]       -S[1,1]-S[2,2]+S[3,3]
    ])
    # This extracts the dominant eigenvector, and interprets it as a Rotor.  In
    # particular, note that the _last_ eigenvector output by `eigen` (the 4th)
    # has the largest eigenvalue.
    return rotor(eigen(M, 4:4).vectors[:, 1]...)
end


@doc raw"""
    align(A, B, [w])

Find a `Rotor` that aligns the set of rotors `A` to a corresponding set `B` by
minimizing the distance between the first set and the rotated second set.

Here, `A` and `B` must be equally sized arrays of `AbstractQuaternion`s.  If
present, `w` must be an equally sized array of real numbers; if not, it is
taken to be 1.  We define the loss function
```math
L(R) ≔ Σᵢ wᵢ |Aᵢ - R Bᵢ|²
```
where ``R`` is a `Rotor`, and return the quaternion corresponding to the optimal
``R`` that minimizes this function.

Note that it is possible that the input data do not uniquely determine a rotor,
which will happen when sum below is zero.  When this happens, the result will
contain `NaN`s, but no error will be raised.  When the sum is very close to —
but not exactly — zero, the accuracy of the result will be limited.  However,
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
∂ⱼ Σᵢ wᵢ |Aᵢ - \exp[vⱼ] R Bᵢ|²  →  2 ⟨ eⱼ R Σᵢ wᵢ Bᵢ Āᵢ ⟩₀
```
where → denotes taking ``vⱼ→0``, the symbol ``⟨⟩₀`` denotes taking the scalar
part, and ``eⱼ`` is the unit quaternionic vector in the ``j`` direction.  The
only way for this quantity to be zero for each choice of ``j`` is if
```math
R Σᵢ wᵢ Bᵢ Āᵢ
```
is itself a pure scalar.  This, in turn, can only happen if either (1) the sum
is 0 or (2) if ``R`` is proportional to the _conjugate_ of the sum:
```math
R ∝ Σᵢ wᵢ Aᵢ B̄ᵢ
```
Now, since we want ``R`` to be a rotor, we simply define it to be the normalized
sum.

"""
function align(A, B, w)
    rotor(sum(w[i] * A[i] * conj(B[i]) for i in eachindex(A, B, w)))
end

function align(A, B)
    rotor(sum(A[i] * conj(B[i]) for i in eachindex(A, B)))
end
