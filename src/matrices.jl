# LinearAlgebra's generic algorithms work for matrices of quaternions when they keep the
# factors of every product in order and rely on `dot(p, q) == conj(p) * q`.  These include
# products, `inv`, `\` and `/`, `lu`, `cholesky`, `qr`, least squares, `svd` (supplied by
# GenericLinearAlgebra), and `eigen` of `Hermitian` matrices.  A few functions instead rely
# on identities that hold only when the entries of a matrix commute, and would silently
# return wrong results for quaternions, so the methods in this file throw an `ArgumentError`
# instead:
#
#   * `det`, `logdet`, and `logabsdet` multiply the pivots of an LU factorization (or the
#     diagonal entries of a triangular matrix), and that product depends on the order of the
#     factors and is not multiplicative.  Issue #119 proposes an implementation.
#   * For a vector `x`, `transpose(x) * A` and `transpose(x) / A` are computed as
#     `transpose(transpose(A) * x)` and `transpose(transpose(A) \ x)`, but (AB)ᵀ ≠ BᵀAᵀ for
#     quaternion matrices.  The same goes for `transpose(x) * A * y` and
#     `muladd(transpose(x), A, z)`.  The conjugate transposes, as in `x' * A`, are correct,
#     because conjugation does reverse products.  `transpose` itself cannot throw, because
#     GenericLinearAlgebra uses transposed rows of quaternion matrices correctly.
#   * `inv(transpose(A))` is computed as `transpose(inv(A))`, and the inverse of a
#     `Symmetric` matrix is assumed to be symmetric, for the same reason.
#   * Solving with the factorization `lu!(A)` of a `Tridiagonal` matrix, which uses
#     LinearAlgebra's specialized solver, and `ldlt` of a `SymTridiagonal` matrix.  Other
#     solves with these matrices use a dense factorization instead (see below).
#
# LinearAlgebra and StaticArrays define methods of some of these functions for particular
# types of matrices, so each of those types needs a method here as well, to avoid
# ambiguities.  Where such a method is needed only to resolve an ambiguity and the operation
# is correct, it computes the result instead of throwing.

function noncommutative_error(operation, alternative)
    throw(ArgumentError(
        "$operation assumes that the entries of the matrix commute, but quaternions do "
        * "not, so the result would be wrong.  $alternative"
    ))
end

## Determinants

const determinant_alternative = (
    "A determinant of matrices of quaternions, such as the Study determinant, is not "
    * "implemented yet; see https://github.com/moble/Quaternionic.jl/issues/119."
)

for f ∈ (:det, :logdet, :logabsdet)
    for M ∈ (
        :AbstractMatrix, :Diagonal, :Bidiagonal, :Tridiagonal, :SymTridiagonal,
        :UpperTriangular, :LowerTriangular, :UnitUpperTriangular, :UnitLowerTriangular,
        :Symmetric, :Hermitian, :UpperHessenberg, :LU,
    )
        @eval LinearAlgebra.$f(::$M{<:AbstractQuaternion}) =
            noncommutative_error($("`$f(A)`"), determinant_alternative)
    end
    @eval LinearAlgebra.$f(::StaticMatrix{<:Any,<:Any,<:AbstractQuaternion}) =
        noncommutative_error($("`$f(A)`"), determinant_alternative)
end
for M ∈ (:UpperTriangular, :LowerTriangular)
    @eval LinearAlgebra.logabsdet(::$M{<:AbstractQuaternion,<:StaticMatrix}) =
        noncommutative_error("`logabsdet(A)`", determinant_alternative)
end

## Transposed vectors

const TransposedQuaternionVector = Transpose{<:AbstractQuaternion,<:AbstractVector}
const TransposedStaticQuaternionVector = Transpose{<:AbstractQuaternion,<:StaticVector}
const StaticTriangular = Union{
    UpperTriangular{<:Any,<:StaticMatrix}, LowerTriangular{<:Any,<:StaticMatrix},
    UnitUpperTriangular{<:Any,<:StaticMatrix}, UnitLowerTriangular{<:Any,<:StaticMatrix},
}
const transposed_product_alternative =
    "Use `permutedims(x) * A` or, for the conjugate transpose, `x' * A`."
const transposed_quotient_alternative =
    "Use `permutedims(x) / A` or, for the conjugate transpose, `x' / A`."

for M ∈ (
    :AbstractMatrix,
    :(Transpose{<:Any,<:Adjoint{<:Any,<:AbstractVector}}),
    :(Adjoint{<:Any,<:Transpose{<:Any,<:AbstractVector}}),
    # Julia 1.10 has specialized methods for these types.
    :Diagonal, :(LinearAlgebra.AbstractTriangular), :(Transpose{<:Any,<:AbstractVector}),
)
    @eval Base.:*(::TransposedQuaternionVector, ::$M) =
        noncommutative_error("`transpose(x) * A`", transposed_product_alternative)
end
Base.:*(::TransposedStaticQuaternionVector, ::StaticTriangular) =
    noncommutative_error("`transpose(x) * A`", transposed_product_alternative)

const Triangular = (
    :UpperTriangular, :LowerTriangular, :UnitUpperTriangular, :UnitLowerTriangular
)
for M ∈ (
    :AbstractMatrix, :Diagonal, :Bidiagonal, :UpperHessenberg,
    :(Union{UpperTriangular,LowerTriangular}),
    :(Union{UnitUpperTriangular,UnitLowerTriangular}),
    (:($T{<:Any,<:$W}) for T ∈ Triangular, W ∈ (:Adjoint, :Transpose))...,
    :(Adjoint{<:Any,<:AbstractMatrix}),
    :(Adjoint{T,<:UpperHessenberg{T}} where {T}),
    :(Transpose{T,<:UpperHessenberg{T}} where {T}),
    # Julia 1.10 has specialized methods for these types.
    :(Adjoint{<:Any,<:Bidiagonal}), :(Transpose{<:Any,<:Bidiagonal}),
)
    @eval Base.:/(::TransposedQuaternionVector, ::$M) =
        noncommutative_error("`transpose(x) / A`", transposed_quotient_alternative)
end
Base.:/(::TransposedStaticQuaternionVector, ::StaticTriangular) =
    noncommutative_error("`transpose(x) / A`", transposed_quotient_alternative)

Base.:*(::TransposedQuaternionVector, ::AbstractMatrix, ::AbstractVector) =
    noncommutative_error("`transpose(x) * A * y`", "Use `transpose(x) * (A * y)`.")
# LinearAlgebra computes these two products correctly; the methods only resolve ambiguities
# with the method above.
Base.:*(x::TransposedQuaternionVector, A::Diagonal, y::AbstractVector) = x * (A * y)
Base.:*(
    x::TransposedQuaternionVector,
    A::Union{Adjoint{<:Any,<:AbstractMatrix},Transpose{<:Any,<:AbstractMatrix}},
    y::AbstractVector,
) = x * (A * y)

Base.muladd(
    ::TransposedQuaternionVector, ::AbstractMatrix, ::Union{Number,AbstractVecOrMat}
) = noncommutative_error(
    "`muladd(transpose(x), A, z)`", "Compute `permutedims(x) * A` first."
)

## Inverses

Base.inv(::Transpose{<:AbstractQuaternion,<:AbstractMatrix}) =
    noncommutative_error("`inv(transpose(A))`", "Use `inv(copy(transpose(A)))`.")
Base.inv(::Symmetric{<:AbstractQuaternion,<:StridedMatrix}) = noncommutative_error(
    "`inv(A)` for a `Symmetric` matrix",
    "The inverse of a symmetric matrix of quaternions need not be symmetric; use "
    * "`inv(Matrix(A))`."
)

## Tridiagonal matrices

# LinearAlgebra's specialized solvers for `Tridiagonal` and `SymTridiagonal` matrices
# multiply some factors in the wrong order for quaternions.  So `lu` and `factorize` of
# these matrices factorize a dense copy instead, `\` with a `SymTridiagonal` matrix uses
# that factorization rather than `ldlt`, and `ldiv!` with either type solves with a dense
# copy.  Every solve then goes through a correct factorization, at a cost of O(n³) rather
# than O(n).  Each method is defined for the two types separately, because a method for
# their union would be ambiguous with LinearAlgebra's methods for each of them.
for M ∈ (:Tridiagonal, :SymTridiagonal)
    @eval begin
        LinearAlgebra.lu(
            A::$M{<:AbstractQuaternion},
            pivot::Union{RowMaximum,NoPivot}=RowMaximum();
            kwargs...
        ) = LinearAlgebra.lu(Matrix(A), pivot; kwargs...)
        LinearAlgebra.factorize(A::$M{<:AbstractQuaternion}) = LinearAlgebra.lu(A)
        LinearAlgebra.ldiv!(A::$M{<:AbstractQuaternion}, B::AbstractVecOrMat) =
            copyto!(B, LinearAlgebra.lu(A) \ B)
    end
end
Base.:\(A::SymTridiagonal{<:AbstractQuaternion}, B::AbstractVecOrMat) =
    LinearAlgebra.lu(A) \ B

# `lu!` still computes the factors in the storage of the `Tridiagonal` matrix, and those
# factors are correct, but LinearAlgebra's solvers that use them are not.  `ldlt` assumes
# that transposition reverses products, which it does not for quaternions.
const tridiagonal_alternative = "Use `lu(A)`, which factorizes a dense copy of `A`."

const tridiagonal_lu_operation =
    "Solving with the factorization `lu!(A)` of a `Tridiagonal` matrix"

LinearAlgebra.ldiv!(
    ::LU{T,Tridiagonal{T,V}}, ::AbstractVecOrMat
) where {T<:AbstractQuaternion,V} =
    noncommutative_error(tridiagonal_lu_operation, tridiagonal_alternative)
for W ∈ (:AdjointFactorization, :TransposeFactorization)
    @eval LinearAlgebra.ldiv!(
        ::LinearAlgebra.$W{<:Any,<:LU{T,Tridiagonal{T,V}}}, ::AbstractVecOrMat
    ) where {T<:AbstractQuaternion,V} =
        noncommutative_error(tridiagonal_lu_operation, tridiagonal_alternative)
end
LinearAlgebra.ldlt!(::SymTridiagonal{<:AbstractQuaternion}) =
    noncommutative_error("`ldlt` of a `SymTridiagonal` matrix", tridiagonal_alternative)

## Cholesky factorization on Julia 1.10

@static if VERSION < v"1.11"
    # COV_EXCL_START
    # Julia 1.10's `cholesky!` handles only `Hermitian` matrices of real or complex numbers;
    # for other element types, its generic method calls itself until the stack overflows.
    # This is the method of later versions, which works for quaternions.  (Coverage is
    # measured with the latest Julia, where this method is not defined.)
    function LinearAlgebra.cholesky!(
        A::Hermitian{<:AbstractQuaternion}, ::LinearAlgebra.NoPivot=LinearAlgebra.NoPivot();
        check::Bool=true
    )
        Triangle = A.uplo == 'U' ? UpperTriangular : LowerTriangular
        C, info = LinearAlgebra._chol!(A.data, Triangle)
        check && LinearAlgebra.checkpositivedefinite(info)
        return LinearAlgebra.Cholesky(C.data, A.uplo, info)
    end
    # COV_EXCL_STOP
end
