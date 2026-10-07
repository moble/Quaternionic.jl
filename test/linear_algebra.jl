# The dot product, and linear algebra with vectors and matrices of quaternions
#
# The references are written out entry by entry, or use the complex representation χ, which
# replaces each quaternion `q = a + b𝐣` (with `a = w + x𝑖` and `b = y + z𝑖`) by the complex
# matrix `[a b; -b̄ ā]`.  χ is an algebra homomorphism, so complex LAPACK gives independent
# reference values for singular values and eigenvalues.

@testsnippet ComplexRepresentation begin
    using LinearAlgebra
    χ(q::AbstractQuaternion) =
        (a = complex(q[1], q[2]); b = complex(q[3], q[4]); [a b; -conj(b) conj(a)])
    χ(A::AbstractMatrix) =
        reduce(vcat, [reduce(hcat, [χ(A[i, j]) for j ∈ axes(A, 2)]) for i ∈ axes(A, 1)])
    # The products of the entries in the order in which they are written
    mulref(X::AbstractMatrix, Y::AbstractMatrix) =
        [sum(X[i, k] * Y[k, j] for k ∈ axes(X, 2)) for i ∈ axes(X, 1), j ∈ axes(Y, 2)]
    relerr(a, b) = norm(a - b) / max(norm(b), 1)
end


@testitem "dot: values for each type of quaternion" tags=[:unit, :fast] begin
    using LinearAlgebra: dot
    import Quaternionic: componentdot
    p, q = quaternion(1.2, -0.7, 0.5, 0.3), quaternion(0.3, 0.8, -0.4, 1.1)
    v, w = quatvec(0.3, -0.6, 0.2), quatvec(-1.1, 0.4, 0.9)
    R, S = rotor(0.5, 0.3, -0.2, 0.4), rotor(-0.8, 0.6, 0.5, -0.7)
    pc = Quaternion{ComplexF64}(1 + 2im, 0.5, -im, 0.3 - 0.2im)
    qc = Quaternion{ComplexF64}(0.2, 1 - im, 0.4im, -0.6)

    # General quaternions: `conj(p) * q`, as LinearAlgebra assumes for numbers
    @test p ⋅ q == conj(p) * q
    @test dot(p, q) == conj(p) * q
    @test R ⋅ S ≈ conj(R) * S
    @test p ⋅ v == conj(p) * v
    @test v ⋅ p == conj(v) * p
    @test pc ⋅ qc == conj(pc) * qc
    @test dot(2.5, q) == 2.5 * q
    @test dot(q, 2.5) == conj(q) * 2.5

    # Two `QuatVec`s: the dot product of the vector parts, which is the scalar part of the
    # general product
    @test v ⋅ w === 0.3 * -1.1 + -0.6 * 0.4 + 0.2 * 0.9
    @test v ⋅ w ≈ real(conj(v) * w)
    @test 𝐢 ⋅ 𝐢 == 1 && 𝐢 ⋅ 𝐣 == 0
    @test (QuatVec{ComplexF64}(0, 1, im, 2) ⋅ QuatVec{ComplexF64}(0, im, 1, 0)) == 2im

    # The scalar part of `p ⋅ q` is the sum of the products of the components, which
    # `unflip` and `slerp` use to choose a hemisphere.
    @test real(p ⋅ q) ≈ componentdot(p, q)
    @test real(R ⋅ S) ≈ componentdot(R, S)
    @test componentdot(R, S) < 0
    @test slerp(R, S, 0.3; unflip=true) ≈ slerp(R, -S, 0.3)
    @test unflip([R, S]) == [R, -S]
end


@testitem "dot: vectors of quaternions" tags=[:unit, :fast] begin
    using LinearAlgebra: dot
    using Random: Xoshiro
    rng = Xoshiro(1)
    x, y = randn(rng, QuaternionF64, 5), randn(rng, QuaternionF64, 5)
    ref = sum(conj(x[i]) * y[i] for i ∈ 1:5)
    @test dot(x, y) ≈ ref
    @test x' * y ≈ ref
    @test x' * y ≉ y' * x
    @test y' * x ≈ conj(ref)
    # `transpose` does not conjugate, and `transpose(x) * y` keeps the order of the factors.
    @test transpose(x) * y ≈ sum(x .* y)
    @test x .* transpose(y) == [x[i] * y[j] for i ∈ 1:5, j ∈ 1:5]
    # Vectors of `QuatVec`s use the dot product of vectors element by element.
    v, w = randn(rng, QuatVecF64, 5), randn(rng, QuatVecF64, 5)
    @test dot(v, w) ≈ sum(vec(v[i]) ⋅ vec(w[i]) for i ∈ 1:5)
end


@testitem "Linear algebra: generic algorithms with matrices of quaternions" tags=[:unit, :validation] setup=[ComplexRepresentation] begin
    using LinearAlgebra
    using Random: Xoshiro
    rng = Xoshiro(2)
    tol = 1e-10
    for n ∈ (2, 3, 6)
        A, B = randn(rng, QuaternionF64, n, n), randn(rng, QuaternionF64, n, n)
        x, y = randn(rng, QuaternionF64, n), randn(rng, QuaternionF64, n)
        H = A' * A + I
        M, b = randn(rng, QuaternionF64, n + 3, n), randn(rng, QuaternionF64, n + 3)

        @test χ(A * B) ≈ χ(A) * χ(B)
        @test χ(A') ≈ χ(A)'
        @test relerr(x' * A * y, sum(conj(x[i]) * A[i, j] * y[j] for i ∈ 1:n, j ∈ 1:n)) < tol
        @test relerr(dot(x, A, y), x' * A * y) < tol
        @test relerr(vec(collect(x' * A)), vec(mulref(permutedims(conj.(x)), A))) < tol
        @test relerr(permutedims(x) * A, mulref(permutedims(x), A)) < tol
        @test relerr(transpose(A) * B, mulref(permutedims(A), B)) < tol
        @test relerr(inv(A) * A, Matrix(I, n, n)) < tol
        @test relerr(A * (A \ x), x) < tol
        @test relerr((x' / A) * A, x') < tol
        F = lu(A)
        @test relerr(F.L * F.U, A[F.p, :]) < tol
        C = cholesky(Hermitian(H))
        @test relerr(C.U' * C.U, H) < tol
        Fq = qr(M)
        @test relerr(Matrix(Fq.Q) * Fq.R, M) < tol
        z = M \ b
        @test relerr(M' * (M * z), M' * b) < tol
        Fs = svd(A)
        @test relerr(Fs.U * Diagonal(Fs.S) * Fs.Vt, A) < tol
        @test sort(repeat(svdvals(A), 2)) ≈ sort(svdvals(χ(A)))
        @test opnorm(A) ≈ opnorm(χ(A))
        @test cond(A) ≈ cond(χ(A))
        @test rank(x * y') == 1
        @test relerr(pinv(M) * M, Matrix(I, n, n)) < tol
        @test sort(repeat(eigvals(Hermitian(H)), 2)) ≈ eigvals(Hermitian(χ(H)))
        Fe = eigen(Hermitian(H))
        @test relerr(H * Fe.vectors, Fe.vectors * Diagonal(Fe.values)) < tol
    end
end


@testitem "Linear algebra: operations that assume commuting entries throw" tags=[:unit, :fast] setup=[ComplexRepresentation] begin
    using LinearAlgebra
    using StaticArrays: @SMatrix, @SVector
    using Random: Xoshiro
    rng = Xoshiro(3)
    n = 4
    A = randn(rng, QuaternionF64, n, n) + 4n * I
    d, e = randn(rng, QuaternionF64, n) .+ 3, randn(rng, QuaternionF64, n - 1)
    x, y = randn(rng, QuaternionF64, n), randn(rng, QuaternionF64, n)
    S = @SMatrix [quaternion(1.0, 0.2, 0.3, 0.4) quaternion(0.1, 0.5, -0.2, 0.3);
                  quaternion(-0.3, 0.1, 0.6, 0.2) quaternion(2.0, -0.4, 0.1, 0.7)]
    s = @SVector [quaternion(0.3, 0.1, -0.2, 0.5), quaternion(-0.4, 0.2, 0.1, 0.3)]
    matrices = (
        A, Diagonal(d), Bidiagonal(d, e, :U), Tridiagonal(e, d, e), SymTridiagonal(d, e),
        UpperTriangular(A), LowerTriangular(A), UnitUpperTriangular(A), UnitLowerTriangular(A),
        Symmetric(A), Hermitian(A' * A), UpperHessenberg(A), transpose(A), A', view(A, :, :),
        UpperTriangular(A'), LowerTriangular(transpose(A)), S, UpperTriangular(S),
    )

    # Determinants
    for M ∈ matrices, f ∈ (det, logdet, logabsdet)
        @test_throws ArgumentError f(M)
    end
    @test_throws ArgumentError det(lu(A))
    @test_throws "issues/119" det(A)

    # Transposed vectors
    for M ∈ matrices
        size(M, 1) == n || continue
        @test_throws ArgumentError transpose(x) * M
        @test_throws ArgumentError transpose(x) / M
    end
    @test_throws ArgumentError transpose(s) * S
    @test_throws ArgumentError transpose(s) * UpperTriangular(S)
    @test_throws ArgumentError transpose(s) / UpperTriangular(S)
    @test_throws ArgumentError transpose(x) * A * y
    @test_throws ArgumentError muladd(transpose(x), A, transpose(y))
    # The alternatives that the messages suggest are correct.
    @test relerr(permutedims(x) * A, mulref(permutedims(x), A)) < 1e-10
    @test relerr((permutedims(x) / A) * A, permutedims(x)) < 1e-10
    @test transpose(x) * (A * y) ≈ sum(x[i] * A[i, j] * y[j] for i ∈ 1:n, j ∈ 1:n)
    # These products are correct, and are computed rather than rejected.
    D = Diagonal(d)
    @test transpose(x) * D * y ≈ sum(x[i] * d[i] * y[i] for i ∈ 1:n)
    @test transpose(x) * A' * y ≈ sum(x[i] * conj(A[j, i]) * y[j] for i ∈ 1:n, j ∈ 1:n)
    @test transpose(x) * transpose(A) * y ≈ sum(x[i] * A[j, i] * y[j] for i ∈ 1:n, j ∈ 1:n)

    # Inverses
    @test_throws ArgumentError inv(transpose(A))
    @test_throws ArgumentError inv(Symmetric(A))
    @test relerr(inv(copy(transpose(A))) * permutedims(A), Matrix(I, n, n)) < 1e-10
    @test relerr(inv(Matrix(Symmetric(A))) * Matrix(Symmetric(A)), Matrix(I, n, n)) < 1e-10

    # Tridiagonal solvers
    T, ST = Tridiagonal(e, d, e), SymTridiagonal(d, e)
    for M ∈ (T, ST)
        @test_throws ArgumentError M \ x
        @test_throws ArgumentError x' / M
        @test_throws ArgumentError inv(M)
        @test_throws ArgumentError ldiv!(copy(M), copy(x))
        @test relerr(Matrix(M) * (Matrix(M) \ x), x) < 1e-10
    end
    @test_throws ArgumentError ldlt(ST)
end
