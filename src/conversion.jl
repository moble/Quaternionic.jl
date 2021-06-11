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
as_float_array(A::AbstractArray{Quaternion{T}}) where {T} = reinterpret(reshape, T, A)


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
    a1 = atan(q.z, q.w)
    a2 = atan(-q.x, q.y)
    [a1+a2, 2*acos(√((q.w^2+q.z^2)/abs2vec(q))), a1-a2]
end
to_euler_angles(q) = to_euler_angles(q...)


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
function from_euler_phases(z)
    zb = √(z[2])  # exp[iβ/2]
    zp = √(z[1] * z[3])  # exp[i(α+γ)/2]
    zm = √(z[1] * conj(z[3]))  # exp[i(α-γ)/2]
    if abs(z[1] - zp * zm) > abs(z[1] + zp * zm)
        zp *= -1
    end
    Quaternion(zb.re * zp.re, -zb.im * zm.im, zb.im * zm.re, zb.re * zp.im)
end
                       
