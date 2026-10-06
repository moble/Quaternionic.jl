"""
    from_float_array(A)

Reinterpret a float array as an array of quaternions

The input array must have an initial dimension whose size is 4, because
successive indices in that dimension will be considered successive components
of the output quaternion.

Note that this returns a view of the original data [via
`reinterpret(reshape,...)`] only if the base type of the input array
`isbitstype`; otherwise, a new array of `Quaternion`s must be created, and the
memory copied.

See also [`to_float_array`](@ref).

"""
function from_float_array(A::AbstractArray{T}) where {T<:Number}
    isbitstype(T) ? from_float_array(Val(true), A) : from_float_array(Val(false), A)
end
function from_float_array(::Val{true}, A::AbstractArray{T}) where {T<:Number}
    reinterpret(reshape, Quaternion{T}, A)
end
function from_float_array(::Val{false}, A::AbstractArray{T}) where {T<:Number}
    @assert size(A, 1)==4 "First dimension of `A` must be 4, not $(size(A, 1))"
    Q = Array{Quaternion{T}}(undef, size(A)[2:end])
    @inbounds for (i, j) in zip(eachindex(Q), Base.Iterators.partition(eachindex(A), 4))
        @views Q[i] = Quaternion{T}(A[j])
    end
    Q
end


"""
    to_float_array(A)
    to_float_array(q)

View a quaternion array as an array of numbers

The output array will have an extra initial dimension whose size is 4, because successive
indices in that dimension correspond to successive components of the quaternion.  The
element type of the output is the base type of the input quaternions; it is not converted to
a float type.

Note that this returns a view of the original data only if the base type of the input array
`isbitstype`; otherwise, a new array of that type must be created, and the memory copied.

When the argument is a single quaternion `q`, the result is a new `Vector` of
its four components, converted to a float type as by `float(q)`.

See also [`from_float_array`](@ref).
"""
function to_float_array(A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Number}
    isbitstype(T) ? to_float_array(Val(true), A) : to_float_array(Val(false), A)
end
function to_float_array(::Val{true}, A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Number}
    reinterpret(reshape, T, A)
end
function to_float_array(::Val{false}, A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Number}
    F = Array{T}(undef, (4, size(A)...))
    @inbounds for (i, j) in zip(eachindex(A), Base.Iterators.partition(eachindex(F), 4))
        @views F[j] .= components(A[i])
    end
    F
end
# WORKAROUND for an Enzyme bug (EnzymeAD/Enzyme.jl#ISSUE_LTS).  The components are copied
# with `Vector` rather than `collect`, because Enzyme (as of 0.13.211, on Julia 1.10)
# crashes in reverse mode on `collect` of an `SVector` inside a closure called by another
# closure.  The result is the same.  Once that bug is fixed, this can be `collect` again.
to_float_array(q::AbstractQuaternion) = Vector(components(float(q)))


@doc raw"""
    to_euler_angles(R)

Open Pandora's Box.

If somebody is trying to make you use Euler angles, tell them no, and walk away, and go and
tell your mum.

You don't want to use Euler angles.  They are awful.  Stay away.  It's one thing to convert
from Euler angles to quaternions; at least you're moving in the right direction.  But to go
the other way?!  It's just not right.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α𝐤/2) * exp(β𝐣/2) * exp(γ𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so this
represents an initial rotation (in the positive sense) through an angle ``γ`` about the axis
``z``, followed by a rotation through ``β`` about the axis ``y``, and a final rotation
through ``α`` about the axis ``z``.  This is equivalent to performing an initial rotation
through ``α`` about the axis ``z``, followed by a rotation through ``β`` about the *rotated*
axis ``y'``, followed by a rotation through ``γ`` about the *twice-rotated* axis ``z''``.
The angles are naturally assumed to be in radians.

The outputs from this function are in these ranges:

  - α ∈ (-2π, 2π]
  - β ∈ [0, π]
  - γ ∈ (-2π, 2π]

This is redundant, and inconsistent with more standard conventions for Euler angles on
``\mathrm{Spin}(3)``:

  - α ∈ [0, 2π)
  - β ∈ [0, π]
  - γ ∈ [0, 4π)

But since these angles will usually be fed into periodic functions, it usually won't matter.
If you really need angles in those ranges, you can always post-process the output of this
function.

NOTE: Before opening an issue reporting something "wrong" with this function, be sure to
read all of [this page](https://github.com/moble/quaternion/wiki/Euler-angles-are-horrible),
*especially* the very last section about opening issues or pull requests.

!!! warning "Derivatives at the Euler singularities"
    The Euler angles are not differentiable functions of the quaternion components at the
    Euler singularities (β=0 or β=π), where only α+γ or α-γ, respectively, is defined.
    Automatic differentiation of this function at exactly those points returns `NaN` for
    some of the derivatives, and those `NaN`s propagate into the derivatives of anything
    computed from the angles, even when that result is itself smooth.

Quaternions with complex components (`Lorentz` rotors) are not supported, and raise an
`ArgumentError`.

# Returns
- `αβγ::SVector{3,T}`

# Raises
- `AllHell` if you try to actually use Euler angles, when you could have been using
  quaternions like a sensible person.

# See Also
- [`from_euler_angles`](@ref): Create quaternion from Euler angles
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
"""
function to_euler_angles(q::AbstractQuaternion)
    q = float(q)
    a0 = 2atan(hypot(q[2], q[3]), hypot(q[1], q[4]))
    a1 = atan(q[4], q[1])
    a2 = atan(-q[2], q[3])
    SVector(a1+a2, a0, a1-a2)
end
function to_euler_angles(q::AbstractQuaternion{<:Complex})
    throw(ArgumentError(complex_conversion_message(to_euler_angles)))
end

# The message for conversions that are defined only for real rotations.
function complex_conversion_message(f)
    "`$f` is not defined for quaternions with complex components (such as `Lorentz` " *
    "rotors), because they do not represent rotations of real three-dimensional space.  " *
    "Note that `to_rotation_matrix` of a `Lorentz` rotor returns its complex-orthogonal " *
    "action on bivectors."
end


"""
    from_euler_angles(α, β, γ)
    from_euler_angles(αβγ)

Come over from the dark side.

Return the `Rotor` corresponding to the Euler angles `α`, `β`, and `γ`, which may also be
given as a single collection `αβγ` of three angles.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α𝐤/2) * exp(β𝐣/2) * exp(γ𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so this
represents an initial rotation (in the positive sense) through an angle ``γ`` about the axis
``z``, followed by a rotation through ``β`` about the axis ``y``, and a final rotation
through ``α`` about the axis ``z``.  This is equivalent to performing an initial rotation
through ``α`` about the axis ``z``, followed by a rotation through ``β`` about the *rotated*
axis ``y'``, followed by a rotation through ``γ`` about the *twice-rotated* axis ``z''``.
The angles are naturally assumed to be in radians.

NOTE: Before opening an issue reporting something "wrong" with this function, be sure to
read all of [this page](https://github.com/moble/quaternion/wiki/Euler-angles-are-horrible),
*especially* the very last section about opening issues or pull requests.

# See Also
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
"""
function from_euler_angles(α, β, γ)
    rotor(
        cos(β/2)*cos((α+γ)/2),
        -sin(β/2)*sin((α-γ)/2),
        sin(β/2)*cos((α-γ)/2),
        cos(β/2)*sin((α+γ)/2)
    )
end
from_euler_angles(αβγ) = from_euler_angles(αβγ...)


"""
    to_euler_phases(q)

Convert input quaternion to complex phases of Euler angles

Interpreting the input quaternion as a rotation (though its normalization scales out), we
can define the complex Euler phases from the Euler angles (α, β, γ) as

    zₐ ≔ exp(i*α)
    zᵦ ≔ exp(i*β)
    zᵧ ≔ exp(i*γ)

These are more useful geometric quantities than the angles themselves — being involved in
computing spherical harmonics and Wigner's 𝔇 matrices — and can be computed from the
components of the corresponding quaternion algebraically (without the use of transcendental
functions).

Note that `to_euler_phases!(z, q)` is supported for backwards compatibility, but because
this function returns an `SVector`, there is probably no advantage to the in-place approach.

Integer quaternions are converted to a float type first.  The normalization of the input
scales out, as long as the squares of its components neither overflow nor underflow.
Quaternions with complex components (`Lorentz` rotors) are not supported, and raise an
`ArgumentError`.

!!! warning "Incorrect derivatives near singularities"
    Note that the derivatives of these phases with respect to the quaternion components are
    ill-defined at the Euler singularities (specifically, when β=0 or β=π).  The problematic
    cases will just return 1 for the phase factor, which will autodiff to zero.
    Fortunately, this should not pollute downstream results if they are grounded in
    geometric quantities.  Specifically, Wigner's 𝔇 matrices can be computed from these
    phases without any issues, because the ``d`` factors should be zero anyway, and the
    product rule will save us.

# Returns
- `z::SVector{Complex{T}}`: complex phases (zₐ, zᵦ, zᵧ) in that order.

# See Also
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function to_euler_phases(R::AbstractQuaternion{T}) where {T}
    # These are not computed with `hypot`, whose derivative is NaN at zero; the derivative
    # of `√` at zero gives the finite results described in the docstring.
    a = R[1]^2 + R[4]^2
    b = R[2]^2 + R[3]^2
    sqrta = √a
    sqrtb = √b
    if iszerovalue(sqrta)
        zp = one(Complex{T})
    else
        zp = Complex{T}(R[1], R[4]) / sqrta  # exp[i(α+γ)/2]
    end
    if iszerovalue(sqrtb)
        zm = one(Complex{T})
    else
        zm = Complex{T}(R[3], -R[2]) / sqrtb  # exp[i(α-γ)/2]
    end
    SVector(
        zp * zm,  # exp[iα]
        Complex{T}((a - b), 2 * sqrta * sqrtb) / (a + b),  # exp[iβ]
        zp * conj(zm),  # exp[iγ]
    )
end
to_euler_phases(R::AbstractQuaternion{<:Integer}) = to_euler_phases(float(R))
function to_euler_phases(R::AbstractQuaternion{<:Complex})
    throw(ArgumentError(complex_conversion_message(to_euler_phases)))
end


"""
    to_euler_phases!(z, R)

Store the Euler phases of the quaternion `R` in the vector `z`, and return `z`.

This is equivalent to `z[:] = to_euler_phases(R)`, and is provided for backwards
compatibility.  See [`to_euler_phases`](@ref) for details.
"""
function to_euler_phases!(z::Array{Complex{T}}, R::AbstractQuaternion{T}) where {T}
    z[:] = to_euler_phases(R)
    z
end


"""
    from_euler_phases(zₐ, zᵦ, zᵧ)
    from_euler_phases(z)

Return the `Rotor` corresponding to these Euler phases.

The complex Euler phases are defined in terms of the Euler angles (α, β, γ) as

    zₐ ≔ exp(i*α)
    zᵦ ≔ exp(i*β)
    zᵧ ≔ exp(i*γ)

These are more useful geometric quantities than the angles themselves — being involved in
computing spherical harmonics and Wigner's 𝔇 matrices — and can be computed from the
components of the corresponding quaternion algebraically (without the use of transcendental
functions).

# Parameters

- `z`: a collection (such as a `Vector` or `SVector`) of three complex numbers, representing
  the complex phases (zₐ, zᵦ, zᵧ) in that order.

# Returns
- `R::Rotor{T}`

# See Also
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function from_euler_phases(zₐ, zᵦ, zᵧ)
    zb = csqrt(zᵦ)  # exp[iβ/2]
    zp = csqrt(zₐ * zᵧ)  # exp[i(α+γ)/2]
    zm = csqrt(zₐ * conj(zᵧ))  # exp[i(α-γ)/2]
    # This comparison is equivalent to one of `abs`, but `abs2` has no singular derivative
    # at zero, which reverse-mode automatic differentiation would otherwise record.
    if abs2(zₐ - zp * zm) > abs2(zₐ + zp * zm)
        zp = -zp
    end
    rotor(real(zb) * real(zp), -imag(zb) * imag(zm), imag(zb) * real(zm), real(zb) * imag(zp))
end
from_euler_phases(z) = from_euler_phases(z...)


"""
    to_spherical_coordinates(q)

Return the spherical coordinates corresponding to this quaternion.

We can treat the quaternion as a transformation taking the ``z`` axis to some direction
``n̂``.  This direction can be described in terms of spherical coordinates (θ, ϕ).  Here, we
use the convention commonly used in physics: θ represents the "polar angle" between the
``z`` axis and the direction ``n̂``, while ϕ represents the "azimuthal angle" between the
``x`` axis and the projection of ``n̂`` into the ``x``-``y`` plane.  Both angles are given
in radians, and the result is an `SVector` of (θ, ϕ).

The outputs are in these ranges:

  - θ ∈ [0, π]
  - ϕ ∈ (-2π, 2π]

As for [`to_euler_angles`](@ref), ϕ is not reduced to a range of length 2π, so `q` and `-q`
may give values of ϕ that differ by 2π.  Quaternions with complex components (`Lorentz`
rotors) are not supported, and raise an `ArgumentError`.

!!! warning "Derivatives at the poles"
    The spherical coordinates are not differentiable functions of the quaternion components
    at the poles (θ=0 or θ=π).  Automatic differentiation of this function at exactly those
    points returns zero for the derivatives of θ, and `NaN` for some of the derivatives of
    ϕ.

"""
function to_spherical_coordinates(q::Q) where {Q<:AbstractQuaternion}
    q = float(q)
    # This is the same expression as in `to_euler_angles`, which is accurate near the poles
    # and independent of the normalization of `q`.  At exactly the poles, the sum of squares
    # replaces `hypot`; it has the same value (zero) there, but a zero derivative rather
    # than a NaN.
    v = iszerovalue(q[2]) && iszerovalue(q[3]) ? q[2]^2 + q[3]^2 : hypot(q[2], q[3])
    s = iszerovalue(q[1]) && iszerovalue(q[4]) ? q[1]^2 + q[4]^2 : hypot(q[1], q[4])
    a0 = 2atan(v, s)
    a1 = atan(q[4], q[1])
    a2 = atan(-q[2], q[3])
    SVector(a0, a1+a2)
end
function to_spherical_coordinates(q::AbstractQuaternion{<:Complex})
    throw(ArgumentError(complex_conversion_message(to_spherical_coordinates)))
end


"""
    from_spherical_coordinates(θ, ϕ)
    from_spherical_coordinates(θϕ)

Return a rotor corresponding to these spherical coordinates, which may also be given as a
single collection `θϕ` of two angles.

Considering (θ, ϕ) as a point ``n̂`` on the sphere, we can also construct a quaternion that
rotates the ``z`` axis onto that point.  Here, we use the convention commonly used in
physics: θ represents the "polar angle" between the ``z`` axis and the direction ``n̂``,
while ϕ represents the "azimuthal angle" between the ``x`` axis and the projection of ``n̂``
into the ``x``-``y`` plane.  Both angles must be given in radians.

"""
function from_spherical_coordinates(θ, ϕ)
    # `sin` and `cos` are called separately, rather than through `sincos`, because Enzyme
    # cannot take second derivatives of `sincos`.
    sϕ, cϕ = sin(ϕ/2), cos(ϕ/2)
    sθ, cθ = sin(θ/2), cos(θ/2)
    rotor(cθ*cϕ, -sθ*sϕ, sθ*cϕ, cθ*sϕ)
end
from_spherical_coordinates(θϕ) = from_spherical_coordinates(θϕ...)


# Implementation notes for `dominant_eigenvector`, which `from_rotation_matrix` and `align`
# call on 4×4 matrices with static storage.  The methods are split on whether LAPACK can
# handle the element type:
#
#   * For `Float16`, `Float32`, and `Float64`, `eigen(M, n:n)` asks LAPACK for the largest
#     eigenpair alone, which is both cheaper and unambiguous.  That range form has no method
#     for other element types.  LAPACK computes `Float16` in `Float32`, so the result is
#     converted back to the input's element type.  For `Float32` and `Float64`, the LAPACK
#     call is isolated in `dominant_eigenvector_lapack`, which takes a plain `Matrix` and
#     the `uplo` character of the `Symmetric` wrapper.  Automatic differentiation packages
#     cannot trace LAPACK, so the extensions give that function derivative rules, which use
#     the formulas in `dominant_eigenvector_pushforward` and `dominant_eigenvector_pullback`
#     below.
#
#   * Otherwise we take the full decomposition and select by `argmax` of the eigenvalues
#     rather than by position.  Selecting `vectors[:, end]` would be wrong: the *only*
#     guarantee `eigen` makes about ordering is what the underlying implementation chooses
#     to provide, and the generic implementations disagree.  `GenericLinearAlgebra` sorts
#     ascending, but `GenericSchur` — which takes precedence for `Symmetric` once it is
#     loaded, and which arrives indirectly with packages such as `DoubleFloats` — does not.
#     Since the eigenvalues of the Bar-Itzhack matrix are `(-1, -1, -1, 3)`, picking the
#     wrong column does not merely lose accuracy; it returns a rotor a full `π/2` away.
#     Selecting by value makes the result independent of which packages happen to be loaded,
#     and scanning four eigenvalues costs nothing beside the decomposition itself.
#
# The generic method converts static storage to a dense `Matrix` first with `_dense`,
# because the generic backends have no `eigen` method for static storage, and on Julia 1.6
# through 1.8 neither does LAPACK's range form.
#
# Element types that wrap a float, such as the dual numbers of ForwardDiff, are handled by
# `refined_dominant_eigenvector` instead.  The generic `eigen` would propagate derivatives
# through every eigenpair, and for an exact rotation matrix the Bar-Itzhack matrix has the
# triply degenerate eigenvalue -1, so first derivatives of the dominant eigenvector were
# correct but second derivatives (nested duals) were not, and were NaN at the identity.
"""
    dominant_eigenvector(M::Symmetric)

Return the unit eigenvector of the real symmetric matrix `M` that belongs to its largest
eigenvalue.

This is the computation at the heart of [`from_rotation_matrix`](@ref) and [`align`](@ref),
each of which builds a symmetric 4×4 matrix whose dominant eigenvector is the optimal rotor.
The sign of the result is arbitrary, as for any eigenvector; those functions choose it
afterward.

The result is a vector with the element type of `M`.  `Float16`, `Float32`, and `Float64`
matrices are handled by LAPACK, and other real element types (such as `BigFloat`) by a
generic eigendecomposition, which needs a package such as `GenericLinearAlgebra` (loaded by
this package).  The function can be differentiated by the AD packages that this package
supports.  For 4×4 matrices whose elements are the number types of AD packages (such as
ForwardDiff's dual numbers), the eigenvector of the values is refined by Newton steps, so
that higher derivatives are also correct; they have been checked through fourth order.
Derivatives assume that the largest eigenvalue is simple.  Complex matrices are not
supported.
"""
function dominant_eigenvector(M::Symmetric{T}) where {T}
    if size(M) == (4, 4) && typeof(value(first(M))) !== T
        return refined_dominant_eigenvector(M)
    end
    λ, V = try
        eigen(_dense(M))
    catch err
        _eigen_failed(T, err)
    end
    V[:, argmax(λ)]
end
function dominant_eigenvector(M::Symmetric{T}) where {T<:Union{Float32,Float64}}
    dominant_eigenvector_lapack(Matrix(parent(M)), M.uplo)
end
function dominant_eigenvector(M::Symmetric{Float16})
    n = size(M, 1)
    v = eigen(_dense(M), n:n).vectors[:, 1]
    eltype(v) === Float16 ? v : Float16.(v)
end
function dominant_eigenvector(M::Symmetric{<:Complex})
    throw(ArgumentError(
        "`from_rotation_matrix` and `align` need matrices or vectors with real "
        * "elements; complex (and `Lorentz`) input is not supported."
    ))
end

# The unit eigenvector of `Symmetric(A, Symbol(uplo))` belonging to its largest eigenvalue,
# computed by LAPACK.  This is the only place where `from_rotation_matrix` and `align` call
# LAPACK for `Float32` and `Float64`, and its arguments are plain arrays, so that a single
# derivative rule for each automatic-differentiation package covers both functions.
function dominant_eigenvector_lapack(
    A::Matrix{T}, uplo::Char
) where {T<:Union{Float32,Float64}}
    n = size(A, 1)
    eigen(Symmetric(A, Symbol(uplo)), n:n).vectors[:, 1]
end

# The first `n = length(v)` entries of the solution `x` of the bordered system
#
#     [S - λI  v] [x]   [b]
#     [ vᵀ     0] [μ] = [0],
#
# which is nonsingular when `λ` is a simple eigenvalue of the symmetric matrix `S` with unit
# eigenvector `v`.  Then `x` is the solution of `(S - λI) x = b` orthogonal to `v`, for any
# `b` orthogonal to `v`.
function bordered_solve(S, λ, v, b)
    n = length(v)
    K = [S - λ * LinearAlgebra.I  v; transpose(v)  zero(λ)]
    (K \ [b; zero(λ)])[1:n]
end

# The derivative of `v = dominant_eigenvector_lapack(A, uplo)` along the tangent `Ȧ` of `A`,
# given `v` itself.  Differentiating `S v = λ v` and `vᵀv = 1`, where `S` is the symmetric
# matrix and `λ = vᵀ S v`, gives `(S - λI) v̇ = λ̇ v - Ṡ v` with `vᵀ v̇ = 0`, so `v̇` is
# minus the bordered solve with `b = Ṡ v`.  The component of `b` along `v` changes only
# `λ̇`, and is removed first.  Only the triangle of `Ȧ` named by `uplo` is used, as for `A`.
function dominant_eigenvector_pushforward(A, uplo::Char, v, Ȧ)
    S = Matrix(Symmetric(A, Symbol(uplo)))
    λ = v ⋅ (S * v)
    b = Matrix(Symmetric(Ȧ, Symbol(uplo))) * v
    -bordered_solve(S, λ, v, b - (v ⋅ b) * v)
end

# The cotangent of `A` in `v = dominant_eigenvector_lapack(A, uplo)`, given `v` and its
# cotangent `v̄`.  The bordered matrix is symmetric, so the pullback applies the same solve
# as the pushforward, sign included, to `v̄`, giving `w`.  The cotangent of the full
# symmetric matrix is `G = w vᵀ`.  Each off-diagonal entry of the triangle named by `uplo`
# stands for two entries of the symmetric matrix, and each diagonal entry for one, so the
# cotangent of `A` is `triu(G + Gᵀ) - Diagonal(G)` for `'U'`, or the `tril` analogue for
# `'L'`, and zero in the other triangle.
function dominant_eigenvector_pullback(A, uplo::Char, v, v̄)
    S = Matrix(Symmetric(A, Symbol(uplo)))
    λ = v ⋅ (S * v)
    w = -bordered_solve(S, λ, v, v̄ - (v ⋅ v̄) * v)
    G = w * transpose(v)
    H = G + transpose(G)
    triangle = uplo == 'U' ? LinearAlgebra.triu(H) : LinearAlgebra.tril(H)
    triangle - LinearAlgebra.Diagonal(G)
end

# The 5×5 matrix `[B v; vᵀ 0]` for a 4×4 matrix `B` and a 4-vector `v`, built element by
# element in column-major order.  Concatenation would be simpler to read, but it hits an
# ambiguous `vcat` method for some element types, such as those of ReverseDiff.
function bordered_matrix(B::SMatrix{4,4,T}, v::SVector{4,T}) where {T}
    SMatrix{5,5,T}(
        B[1,1], B[2,1], B[3,1], B[4,1], v[1],
        B[1,2], B[2,2], B[3,2], B[4,2], v[2],
        B[1,3], B[2,3], B[3,3], B[4,3], v[3],
        B[1,4], B[2,4], B[3,4], B[4,4], v[4],
        v[1], v[2], v[3], v[4], zero(T)
    )
end

# The dominant eigenvector of a 4×4 symmetric matrix `M` whose element type `T` wraps a
# float, such as a (possibly nested) dual number.  The eigenvector `v₀` of the matrix of
# values, `map(value, M)`, is computed by the float methods above; it is exact for the
# values but carries no derivatives.  Then two Newton steps on `F(v, λ) = (M v - λ v, (vᵀv -
# 1)/2)` in the element type `T` fill in the derivatives.  As a power series in the
# perturbation that the derivatives describe, the error of `v₀` begins at first order, and
# each Newton step doubles the order at which the error begins.  So one step gives exact
# first derivatives, and two steps give exact derivatives through third order.  The Jacobian
# of `F` is `[M-λI -v; vᵀ 0]`, which is nonsingular when the dominant eigenvalue is simple,
# as it is for any (nearly) orthogonal input to `from_rotation_matrix`.  This avoids
# differentiating the degenerate eigenvectors that the generic `eigen` would otherwise
# differentiate.
function refined_dominant_eigenvector(M::Symmetric{T}) where {T}
    A = SMatrix{4,4,T}(M)
    v = SVector{4,T}(dominant_eigenvector(Symmetric(map(value, A))))
    λ = v ⋅ (A * v)
    for _ ∈ 1:2
        r = A * v - λ * v
        # The Newton step `(δv, δλ)` solves `[A-λI  -v; vᵀ  0] (δv, δλ) = -F`, which is
        # the symmetric system below for `(δv, -δλ)`.
        δ = bordered_matrix(A - λ * LinearAlgebra.I, v) \
            SVector{5,T}(-r[1], -r[2], -r[3], -r[4], (1 - v ⋅ v) / 2)
        v += SVector{4,T}(δ[1], δ[2], δ[3], δ[4])
        λ -= δ[5]
    end
    v
end

# Give a usable error for element types that cannot support an
# eigen-decomposition at all.  Symbolic types are the realistic case:
# `Symmetric{Symbolics.Num}` gets several steps into the QR iteration
# before dying with `TypeError: non-boolean (Num) used in boolean
# context`, which tells the caller nothing about what they did wrong.
# `AbstractFloat` and `Integer` are passed through untouched — LAPACK
# and `GenericLinearAlgebra` between them cover every float, and
# integers are promoted to floats by `eigen` — so a failure for those
# is a genuine bug, and more useful seen unfiltered.
@noinline function _eigen_failed(::Type{T}, err) where {T}
    (T <: AbstractFloat || T <: Integer) && rethrow(err)
    # Only the first line of the original error: some of them
    # (Symbolics') carry several paragraphs of their own advice, which
    # would bury ours.
    summary = first(split(sprint(showerror, err), '\n'))
    throw(ArgumentError("""
        Cannot compute the eigen-decomposition needed here for element
        type `$T`.

        Finding a dominant eigenvector is inherently iterative, so it
        needs a floating-point element type.  Integers are promoted to
        floats automatically, but symbolic types — and any other type
        that cannot answer the comparisons the iteration performs —
        cannot be.

        Convert the inputs to a float type first.  `float.(x)` is
        usually enough; use `BigFloat` or `Double64` if you need more
        precision than `Float64` provides.

        The underlying failure was `$summary`."""))
end

# Materialize static storage for `eigen`.  This is a separate helper rather than another
# `dominant_eigenvector` method so that `Symmetric{Float64,<:SMatrix}` — which would match
# both an element-type method and a storage-type method — cannot become an ambiguity.
#
# The generic backends implement `eigen` only against dense storage:
# `eigen(::Symmetric{Double64,<:SMatrix})` raises a `MethodError`, even though the same
# matrix decomposes fine once it is dense.  On Julia 1.6 through 1.8, the LAPACK range form
# `eigen(::Symmetric{Float64,<:SMatrix}, n:n)` fails in the same way.  A 4×4 copy is
# negligible beside an eigendecomposition, and it is what lets `from_rotation_matrix` and
# `align` — which build an `SMatrix` — accept the full range of float types.
_dense(M::Symmetric) = M
_dense(M::Symmetric{<:Any,<:SMatrix}) = Symmetric(Matrix(M))


"""
    from_rotation_matrix(ℛ)

Convert 3x3 rotation matrix to quaternion.

Assuming the 3x3 matrix `ℛ` rotates a vector `v` according to

    v' = ℛ * v,

we can also express this rotation in terms of a quaternion `R` such that

    v' = R * v * R⁻¹.

This function returns that quaternion, using Bar-Itzhack's algorithm (version 3) to allow
for non-orthogonal matrices.  [J. Guidance, Vol. 23, No. 6, p.
1085](http://dx.doi.org/10.2514/2.4654)

Since `R` and `-R` represent the same rotation, the sign of the result is fixed by
convention: its scalar part is positive, or, if the scalar part is zero, its first nonzero
component is positive.  The result is therefore discontinuous at rotations through an angle
of π, where the scalar part changes sign.  See [`unflip`](@ref) for making a sequence of
rotors continuous.

This works for any real element type, and for float input the result has the same element
type as the input.  It computes the eigen-decomposition of a 4×4 symmetric matrix built from
`ℛ`, which LAPACK provides only for `Float16`, `Float32`, and `Float64`; for anything else
we rely on `GenericLinearAlgebra`, which `Quaternionic` depends on directly.  The
eigen-decomposition makes this function far more expensive than
[`to_rotation_matrix`](@ref).  Matrices with complex elements are not supported.

"""
function from_rotation_matrix(ℛ::AbstractMatrix)
    @assert size(ℛ) == (3, 3)
    @inbounds begin
        # Compute 3K₃ according to Eq. (2) of Bar-Itzhack.  We will just be looking for the
        # eigenvector with the largest eigenvalue, so scaling by a strictly positive number
        # (3, in this case) won't change that.  Only the upper triangle is used, and the
        # elements are given in column-major order.
        z = zero(ℛ[1,1])
        K₃3 = Symmetric(SMatrix{4,4}(
            ℛ[1,1]-ℛ[2,2]-ℛ[3,3], z, z, z,
            ℛ[2,1]+ℛ[1,2], ℛ[2,2]-ℛ[1,1]-ℛ[3,3], z, z,
            ℛ[3,1]+ℛ[1,3], ℛ[3,2]+ℛ[2,3], ℛ[3,3]-ℛ[1,1]-ℛ[2,2], z,
            ℛ[2,3]-ℛ[3,2], ℛ[3,1]-ℛ[1,3], ℛ[1,2]-ℛ[2,1], ℛ[1,1]+ℛ[2,2]+ℛ[3,3]
        ))

        # Compute the *dominant* eigenvector (the one with the largest eigenvalue)
        de = dominant_eigenvector(K₃3)

        # Convert it into a quaternion
        R = rotor(de[4], -de[1], -de[2], -de[3])
    end
    positive_hemisphere(R)
end

# Choose the sign of the rotor `R`, which represents the same rotation as `-R`, so that its
# first nonzero component is positive.  The comparisons use only the values of the
# components, so that dual numbers make the same choice as the floats they wrap.
function positive_hemisphere(R::Rotor)
    for c ∈ components(R)
        v = value(c)
        if !iszero(v)
            return v < 0 ? -R : R
        end
    end
    R
end


"""
    to_rotation_matrix(q)

Convert quaternion to 3x3 rotation matrix.

Assuming the quaternion `R` rotates a vector `v` according to

    v' = R * v * R⁻¹,

we can also express this rotation in terms of a 3x3 matrix `ℛ` such that

    v' = ℛ * v.

This function returns that matrix, as an `SMatrix{3,3}`.  The normalization of `q` scales
out, as long as `abs2(q)` neither overflows nor underflows.

For a `Lorentz` rotor (a quaternion with complex components), the result is the
complex-orthogonal matrix in SO(3,ℂ) that represents the action of the rotor on bivectors.

"""
function to_rotation_matrix(q::Q) where {Q<:AbstractQuaternion}
    n = inv(abs2(q))
    # The elements are given in column-major order.
    SMatrix{3,3}(
        1 - 2*(q[3]^2 + q[4]^2) * n,
        2*(q[2]*q[3] + q[4]*q[1]) * n,
        2*(q[2]*q[4] - q[3]*q[1]) * n,
        2*(q[2]*q[3] - q[4]*q[1]) * n,
        1 - 2*(q[2]^2 + q[4]^2) * n,
        2*(q[3]*q[4] + q[2]*q[1]) * n,
        2*(q[2]*q[4] + q[3]*q[1]) * n,
        2*(q[3]*q[4] - q[2]*q[1]) * n,
        1 - 2*(q[2]^2 + q[3]^2) * n
    )
end
