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

View a quaternion array as an array of numbers

The output array will have an extra initial dimension whose size is 4, because
successive indices in that dimension correspond to successive components of the
quaternion.

Note that this returns a view of the original data only if the base type of the
input array `isbitstype`; otherwise, a new array of that type must be created,
and the memory copied.

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
to_float_array(q::AbstractQuaternion) = collect(float(q).components)


"""
    to_euler_angles(R)

Open Pandora's Box.

If somebody is trying to make you use Euler angles, tell them no, and walk
away, and go and tell your mum.

You don't want to use Euler angles.  They are awful.  Stay away.  It's one
thing to convert from Euler angles to quaternions; at least you're moving in
the right direction.  But to go the other way?!  It's just not right.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α𝐤/2) * exp(β𝐣/2) * exp(γ𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation (in the positive sense) through an angle
``γ`` about the axis ``z``, followed by a rotation through ``β`` about the axis
``y``, and a final rotation through ``α`` about the axis ``z``.  This is
equivalent to performing an initial rotation through ``α`` about the axis
``z``, followed by a rotation through ``β`` about the *rotated* axis ``y'``,
followed by a rotation through ``γ`` about the *twice-rotated* axis ``z''``.
The angles are naturally assumed to be in radians.

NOTE: Before opening an issue reporting something "wrong" with this function,
be sure to read all of [this
page](https://github.com/moble/quaternion/wiki/Euler-angles-are-horrible),
*especially* the very last section about opening issues or pull requests.

# Returns
- `αβγ::Vector{T}`

# Raises
- `AllHell` if you try to actually use Euler angles, when you could have been
  using quaternions like a sensible person.

# See Also
- [`from_euler_angles`](@ref): Create quaternion from Euler angles
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
"""
function to_euler_angles(q::AbstractQuaternion)
    q = float(q)
    a0 = 2acos(√((q[1]^2+q[4]^2)/abs2(q)))
    a1 = atan(q[4], q[1])
    a2 = atan(-q[2], q[3])
    @SVector [a1+a2, a0, a1-a2]
end


"""
    from_euler_angles(α, β, γ)

Come over from the dark side.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α𝐤/2) * exp(β𝐣/2) * exp(γ𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation (in the positive sense) through an angle
``γ`` about the axis ``z``, followed by a rotation through ``β`` about the axis
``y``, and a final rotation through ``α`` about the axis ``z``.  This is
equivalent to performing an initial rotation through ``α`` about the axis
``z``, followed by a rotation through ``β`` about the *rotated* axis ``y'``,
followed by a rotation through ``γ`` about the *twice-rotated* axis ``z''``.
The angles are naturally assumed to be in radians.

NOTE: Before opening an issue reporting something "wrong" with this function,
be sure to read all of [this
page](https://github.com/moble/quaternion/wiki/Euler-angles-are-horrible),
*especially* the very last section about opening issues or pull requests.

# See Also
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
"""
function from_euler_angles(α, β, γ)
    Rotor(
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

Interpreting the input quaternion as a rotation (though its normalization
scales out), we can define the complex Euler phases from the Euler angles (α,
β, γ) as

    zₐ ≔ exp(i*α)
    zᵦ ≔ exp(i*β)
    zᵧ ≔ exp(i*γ)

These are more useful geometric quantites than the angles themselves — being
involved in computing spherical harmonics and Wigner's 𝔇 matrices — and can be
computed from the components of the corresponding quaternion algebraically
(without the use of transcendental functions).

Note that `to_euler_phases!(z, q)` is supported for backwards compatibility,
but because this function returns an `SVector`, there is probably no advantage
to the in-place approach.

# Returns
- `z::SVector{Complex{T}}`: complex phases (zₐ, zᵦ, zᵧ) in that order.

# See Also
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function to_euler_phases(R::AbstractQuaternion{T}) where {T}
    a = R[1]^2 + R[4]^2
    b = R[2]^2 + R[3]^2
    sqrta = √a
    sqrtb = √b
    if iszero(sqrta)
        zp = one(Complex{T})
    else
        zp = Complex{T}(R[1], R[4]) / sqrta  # exp[i(α+γ)/2]
    end
    if iszero(sqrtb)
        zm = one(Complex{T})
    else
        zm = Complex{T}(R[3], -R[2]) / sqrtb  # exp[i(α-γ)/2]
    end
    @SVector [
        zp * zm,  # exp[iα]
        Complex{T}((a - b), 2 * sqrta * sqrtb) / (a + b),  # exp[iβ]
        zp * conj(zm),  # exp[iγ]
    ]
end


function to_euler_phases!(z::Array{Complex{T}}, R::AbstractQuaternion{T}) where {T}
    z[:] = to_euler_phases(R)
    z
end


"""
    from_euler_phases(zₐ, zᵦ, zᵧ)
    from_euler_phases(z)

Return the quaternion corresponding to these Euler phases.

Interpreting the input quaternion as a rotation (though its normalization
scales out), we can define the complex Euler phases from the Euler angles (α,
β, γ) as

    zₐ ≔ exp(i*α)
    zᵦ ≔ exp(i*β)
    zᵧ ≔ exp(i*γ)

These are more useful geometric quantites than the angles themselves — being
involved in computing spherical harmonics and Wigner's 𝔇 matrices — and can be
computed from the components of the corresponding quaternion algebraically
(without the use of transcendental functions).

# Parameters

- `z::Vector{Complex{T}}`: complex vector of length 3, representing the complex
  phases (zₐ, zᵦ, zᵧ) in that order.

# Returns
- `R::Quaternion{T}`

# See Also
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function from_euler_phases(zₐ, zᵦ, zᵧ)
    zb = √(zᵦ)  # exp[iβ/2]
    zp = √(zₐ * zᵧ)  # exp[i(α+γ)/2]
    zm = √(zₐ * conj(zᵧ))  # exp[i(α-γ)/2]
    if abs(zₐ - zp * zm) > abs(zₐ + zp * zm)
        zp *= -1
    end
    Rotor(real(zb) * real(zp), -imag(zb) * imag(zm), imag(zb) * real(zm), real(zb) * imag(zp))
end
from_euler_phases(z) = from_euler_phases(z...)


"""
    to_spherical_coordinates(q)

Return the spherical coordinates corresponding to this quaternion.

We can treat the quaternion as a transformation taking the ``z`` axis to some
direction ``n̂``.  This direction can be described in terms of spherical
coordinates (θ, ϕ).  Here, we use the standard commonly used in physics: θ
represents the "polar angle" between the ``z`` axis and the direction ``n̂``,
while ϕ represents the "azimuthal angle" between the ``x`` axis and the
projection of ``n̂`` into the ``x``-``y`` plane.  Both angles are given in
radians.

"""
function to_spherical_coordinates(q::Q) where {Q<:AbstractQuaternion}
    q = float(q)
    a0 = 2acos(√((q[1]^2+q[4]^2)/abs2(q)))
    a1 = atan(q[4], q[1])
    a2 = atan(-q[2], q[3])
    @SVector [a0, a1+a2]
end


"""
    from_spherical_coordinates(θ, ϕ)

Return a rotor corresponding to these spherical coordinates.

Considering (θ, ϕ) as a point ``n̂`` on the sphere, we can also construct a
quaternion that rotates the ``z`` axis onto that point.  Here, we use the
standard commonly used in physics: θ represents the "polar angle" between the
``z`` axis and the direction ``n̂``, while ϕ represents the "azimuthal angle"
between the ``x`` axis and the projection of ``n̂`` into the ``x``-``y`` plane.
Both angles must be given in radians.

"""
function from_spherical_coordinates(θ, ϕ)
    sϕ, cϕ = sincos(ϕ/2)
    sθ, cθ = sincos(θ/2)
    Rotor(cθ*cϕ, -sθ*sϕ, sθ*cϕ, cθ*sϕ)
end
from_spherical_coordinates(θϕ) = from_spherical_coordinates(θϕ...)


dominant_eigenvector(M::Symmetric{T,SMatrix{4,4,T,16}}) where T = eigen(M).vectors[:, 4]
dominant_eigenvector(M::Symmetric{Float64,SMatrix{4,4,Float64,16}}) = eigen(M, 4:4).vectors[:, 1]
dominant_eigenvector(M::Symmetric{Float32,SMatrix{4,4,Float32,16}}) = eigen(M, 4:4).vectors[:, 1]


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

!!! note
    If you want to use this function for matrices with elements of types other than
    `Float64` or `Float32`, you will need to (install and) import `GenericLinearAlgebra`
    first.  The reason is that this function computes the eigen-decomposition of `ℛ`, which
    is only available for more generic float types via that package.

"""
function from_rotation_matrix(ℛ::AbstractMatrix)
    @assert size(ℛ) == (3, 3)
    @inbounds begin
        # Compute 3K₃ according to Eq. (2) of Bar-Itzhack.  We will just be looking for the
        # eigenvector with the largest eigenvalue, so scaling by a strictly positive number
        # (3, in this case) won't change that.
        K₃3 = Symmetric(@SMatrix[
            ℛ[1,1]-ℛ[2,2]-ℛ[3,3]  ℛ[2,1]+ℛ[1,2]         ℛ[3,1]+ℛ[1,3]         ℛ[2,3]-ℛ[3,2];
            0                     ℛ[2,2]-ℛ[1,1]-ℛ[3,3]  ℛ[3,2]+ℛ[2,3]         ℛ[3,1]-ℛ[1,3];
            0                     0                     ℛ[3,3]-ℛ[1,1]-ℛ[2,2]  ℛ[1,2]-ℛ[2,1];
            0                     0                     0               ℛ[1,1]+ℛ[2,2]+ℛ[3,3]
        ])

        # Compute the *dominant* eigenvector (the one with the largest eigenvalue)
        de = dominant_eigenvector(K₃3)

        # Convert it into a quaternion
        R = Rotor(de[4], -de[1], -de[2], -de[3])
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

This function returns that matrix.

"""
function to_rotation_matrix(q::Q) where {Q<:AbstractQuaternion}
    n = inv(abs2(q))
    @SMatrix [
        1 - 2*(q[3]^2 + q[4]^2) * n  2*(q[2]*q[3] - q[4]*q[1]) * n  2*(q[2]*q[4] + q[3]*q[1]) * n ;
        2*(q[2]*q[3] + q[4]*q[1]) * n  1 - 2*(q[2]^2 + q[4]^2) * n  2*(q[3]*q[4] - q[2]*q[1]) * n ;
        2*(q[2]*q[4] - q[3]*q[1]) * n  2*(q[3]*q[4] + q[2]*q[1]) * n  1 - 2*(q[2]^2 + q[3]^2) * n
    ]
end
