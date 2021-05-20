"""
    as_quat_array(A)

View a real array as an array of quaternions

The input array must have an initial dimension whose size is
divisible by four (or better yet *is* 4), because successive
indices in that last dimension will be considered successive
components of the output quaternion.
"""
as_quat_array(A::AbstractArray{T}) where {T<:Real} = reinterpret(reshape, Quaternion{T}, A)


"""
    as_float_array(A)

View a quaternion array as an array of real numbers

This function is fast because no data is copied; the returned quantity is just a "view"
of the original.

The output view will have an extra initial dimension (of size 4), but is otherwise the
same shape as the input array.

"""
as_float_array(A::AbstractArray{Quaternion{T}}) where {T} = reinterpret(reshape, T, A)


"""
    to_euler_phases(q)

Convert input quaternion to complex phases of Euler angles

Returns
-------
z : complex array
    This array contains the complex phases (zₐ, zᵦ, zᵧ) in that order.

See Also
--------
from_euler_phases : Create quaternion from Euler phases
to_euler_angles : Convert quaternion to Euler angles
from_euler_angles : Create quaternion from Euler angles

Notes
-----
We define the Euler phases from the Euler angles (α, β, γ) as

    zₐ ≔ exp(i*α)
    zᵦ ≔ exp(i*β)
    zᵧ ≔ exp(i*γ)

These are more useful geometric quantites than the angles themselves — being
involved in computing spherical harmonics and Wigner's 𝔇 matrices — and can be
computed from the components of the corresponding quaternion algebraically
(without the use of transcendental functions).

"""
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


function to_euler_phases(R::Quaternion{T}) where {T}
    z = Array{Complex{T}}(undef, 3)
    to_euler_phases!(z, R)
    z
end
