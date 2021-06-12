"""
    as_quat_array(A)

View a real array as an array of quaternions

The input array must have an initial dimension whose size is divisible by four
(or better yet *is* 4), because successive indices in that last dimension will
be considered successive components of the output quaternion.
"""
as_quat_array(A::AbstractArray{T}) where {T<:Real} = reinterpret(reshape, Quaternion{T}, A)


"""
    as_float_array(A)

View a quaternion array as an array of real numbers

This function is fast because no data is copied; the returned quantity is just
a "view" of the original.  The output view will have an extra initial dimension
(of size 4), but is otherwise the same shape as the input array.
"""
as_float_array(A::AbstractArray{Quaternion{T}}) where {T<:AbstractFloat} = reinterpret(reshape, T, A)
as_float_array(q::Quaternion) = collect(float(q).components)


"""
    to_euler_angles(R)

Open Pandora's Box.

If somebody is trying to make you use Euler angles, tell them no, and walk
away, and go and tell your mum.

You don't want to use Euler angles.  They are awful.  Stay away.  It's one
thing to convert from Euler angles to quaternions; at least you're moving in
the right direction.  But to go the other way?!  It's just not right.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α 𝐤/2) * exp(β 𝐣/2) * exp(γ 𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation about the ``z`` axis (in the positive sense)
through an angle γ, followed by a rotation about the ``y`` axis by β, and a
final rotation about the ``z`` axis by α.  This is equivalent to performing an
initial rotation about ``z`` by α, followed by a rotation about *the rotated*
``y'`` axis by β, followed by a rotation about *the twice-rotated* ``z''`` axis
by γ.  The angles are naturally in radians.

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
function to_euler_angles(q::Quaternion)
    q = float(q)
    a0 = 2acos(√((q.w^2+q.z^2)/abs2(q)))
    a1 = atan(q.z, q.w)
    a2 = atan(-q.x, q.y)
    [a1+a2, a0, a1-a2]
end


"""
    from_euler_angles(α, β, γ)

Improve your life drastically.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(α 𝐤/2) * exp(β 𝐣/2) * exp(γ 𝐤/2)

where 𝐣 and 𝐤 rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation about the ``z`` axis (in the positive sense)
through an angle γ, followed by a rotation about the ``y`` axis by β, and a
final rotation about the ``z`` axis by α.  This is equivalent to performing an
initial rotation about ``z`` by α, followed by a rotation about *the rotated*
``y'`` axis by β, followed by a rotation about *the twice-rotated* ``z''`` axis
by γ.  The angles naturally must be in radians for this to make any sense.

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
    Quaternion(
        cos(β/2)*cos((α+γ)/2),
        -sin(β/2)*sin((α-γ)/2),
        sin(β/2)*cos((α-γ)/2),
        cos(β/2)*sin((α+γ)/2)
    )
end
from_euler_angles(αβγ) = from_euler_angles(αβγ...)


function to_euler_phases!(z::Array{Complex{T}}, R::Quaternion{T}) where {T}
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
    z[1] = zp * zm  # exp[iα]
    z[2] = Complex{T}((a - b), 2 * sqrta * sqrtb) / (a + b)  # exp[iβ]
    z[3] = zp * conj(zm)  # exp[iγ]
    z
end


"""
    to_euler_phases(q)
    to_euler_phases!(z, q)

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

# Returns
- `z::Vector{Complex{T}}`: complex phases (zₐ, zᵦ, zᵧ) in that order.

# See Also
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function to_euler_phases(R::Quaternion{T}) where {T}
    z = Array{Complex{T}}(undef, 3)
    to_euler_phases!(z, R)
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
    Quaternion(zb.re * zp.re, -zb.im * zm.im, zb.im * zm.re, zb.re * zp.im)
end
from_euler_phases(z) = from_euler_phases(z...)
                       

"""
    to_spherical_coordinates(q)

Return the spherical coordinates corresponding to this quaternion.

We can treat the quaternion as a transformation taking the ``z`` axis to some direction ``n̂``.  This
direction can be described in terms of spherical coordinates (θ, ϕ).  Here, we use the standard
commonly used in physics: θ represents the "polar angle" between the ``z`` axis and the direction
``n̂``, while ϕ represents the "azimuthal angle" between the ``x`` axis and the projection of ``n̂``
into the ``x``-``y`` plane.  Both angles are given in radians.

"""
function to_spherical_coordinates(q::Quaternion)
    q = float(q)
    a0 = 2acos(√((q.w^2+q.z^2)/abs2(q)))
    a1 = atan(q.z, q.w)
    a2 = atan(-q.x, q.y)
    [a0, a1+a2]
end


"""
    from_spherical_coordinates(θ, ϕ)

Return a quaternion corresponding to these spherical coordinates.

Considering (θ, ϕ) as a point ``n̂`` on the sphere, we can also construct a quaternion that rotates
the ``z`` axis onto that point.  Here, we use the standard commonly used in physics: θ represents
the "polar angle" between the ``z`` axis and the direction ``n̂``, while ϕ represents the "azimuthal
angle" between the ``x`` axis and the projection of ``n̂`` into the ``x``-``y`` plane.  Both angles
must be given in radians.

"""
function from_spherical_coordinates(θ, ϕ)
    sϕ, cϕ = sincos(ϕ/2)
    sθ, cθ = sincos(θ/2)
    Quaternion(cθ*cϕ, -sθ*sϕ, sθ*cϕ, cθ*sϕ)
end
from_spherical_coordinates(θϕ) = from_spherical_coordinates(θϕ...)


"""
    from_rotation_matrix(ℛ)

Convert 3x3 rotation matrix to quaternion.

Assuming the 3x3 matrix `ℛ` rotates a vector `v` according to

    v' = ℛ * v,

we can also express this rotation in terms of a quaternion `R` such that

    v' = R * v * R⁻¹.

This function returns that quaternion, using Bar-Itzhack's algorithm to allow for non-orthogonal
matrices.  [J. Guidance, Vol. 23, No. 6, p. 1085](http://dx.doi.org/10.2514/2.4654)

"""
function from_rotation_matrix(ℛ)
    K = Array{eltype(ℛ), 2}(undef, (4, 4))
    K[1, 1] = (ℛ[1, 1] - ℛ[2, 2] - ℛ[3, 3])/3
    K[1, 2] = (ℛ[2, 1] + ℛ[1, 2])/3
    K[1, 3] = (ℛ[3, 1] + ℛ[1, 3])/3
    K[1, 4] = (ℛ[2, 3] - ℛ[3, 2])/3
    K[2, 2] = (ℛ[2, 2] - ℛ[1, 1] - ℛ[3, 3])/3
    K[2, 3] = (ℛ[3, 2] + ℛ[2, 3])/3
    K[2, 4] = (ℛ[3, 1] - ℛ[1, 3])/3
    K[3, 3] = (ℛ[3, 3] - ℛ[1, 1] - ℛ[2, 2])/3
    K[3, 4] = (ℛ[1, 2] - ℛ[2, 1])/3
    K[4, 4] = (ℛ[1, 1] + ℛ[2, 2] + ℛ[3, 3])/3
    H = Hermitian(K)

    # compute the *dominant* (largest eigenvalue) eigenvector
    eigenvec = eigen(transpose(H), 4:4).vectors[:, 1]

    # convert it into a quaternion
    Quaternion(eigenvec[4], -eigenvec[1], -eigenvec[2], -eigenvec[3])
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
function to_rotation_matrix(q)
    m = Array{eltype(q)}(undef, 3, 3)
    n = inv(abs2(q))
    m[1, 1] = 1 - 2*(q.y^2 + q.z^2) * n
    m[1, 2] = 2*(q.x*q.y - q.z*q.w) * n
    m[1, 3] = 2*(q.x*q.z + q.y*q.w) * n
    m[2, 1] = 2*(q.x*q.y + q.z*q.w) * n
    m[2, 2] = 1 - 2*(q.x^2 + q.z^2) * n
    m[2, 3] = 2*(q.y*q.z - q.x*q.w) * n
    m[3, 1] = 2*(q.x*q.z - q.y*q.w) * n
    m[3, 2] = 2*(q.y*q.z + q.x*q.w) * n
    m[3, 3] = 1 - 2*(q.x^2 + q.y^2) * n
    m
end
