"""
    from_float_array(A)

Reinterpret a real array as an array of quaternions

The input array must have an initial dimension whose size is 4, because
successive indices in that dimension will be considered successive components
of the output quaternion.

Note that this returns a view of the original data [via
`reinterpret(reshape,...)`] only if the base type of the input array
`isbitstype`; otherwise, a new array of `Quaternion`s must be created, and the
memory copied.

See also [`to_float_array`](@ref).

"""
function from_float_array(A::AbstractArray{T}) where {T<:Real}
    isbitstype(T) ? from_float_array(Val(true), A) : from_float_array(Val(false), A)
end
function from_float_array(::Val{true}, A::AbstractArray{T}) where {T<:Real}
    reinterpret(reshape, Quaternion{T}, A)
end
function from_float_array(::Val{false}, A::AbstractArray{T}) where {T<:Real}
    @assert size(A, 1)==4 "First dimension of `A` must be 4, not $(size(A, 1))"
    Q = Array{Quaternion{T}}(undef, size(A)[2:end])
    @inbounds for (i, j) in zip(eachindex(Q), Base.Iterators.partition(eachindex(A), 4))
        @views Q[i] = Quaternion{T}(A[j])
    end
    Q
end


"""
    to_float_array(A)

View a quaternion array as an array of real numbers

The output array will have an extra initial dimension whose size is 4, because
successive indices in that dimension correspond to successive components of the
quaternion.

Note that this returns a view of the original data only if the base type of the
input array `isbitstype`; otherwise, a new array of that type must be created,
and the memory copied.

See also [`from_float_array`](@ref).
"""
function to_float_array(A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Real}
    isbitstype(T) ? to_float_array(Val(true), A) : to_float_array(Val(false), A)
end
function to_float_array(::Val{true}, A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Real}
    reinterpret(reshape, T, A)
end
function to_float_array(::Val{false}, A::AbstractArray{<:AbstractQuaternion{T}}) where {T<:Real}
    F = Array{T}(undef, (4, size(A)...))
    @inbounds for (i, j) in zip(eachindex(A), Base.Iterators.partition(eachindex(F), 4))
        @views F[j] .= A[i].components[:]
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

    R = exp(??????/2) * exp(??????/2) * exp(??????/2)

where ???? and ???? rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation (in the positive sense) through an angle
``??`` about the axis ``z``, followed by a rotation through ``??`` about the axis
``y``, and a final rotation through ``??`` about the axis ``z``.  This is
equivalent to performing an initial rotation through ``??`` about the axis
``z``, followed by a rotation through ``??`` about the *rotated* axis ``y'``,
followed by a rotation through ``??`` about the *twice-rotated* axis ``z''``.
The angles are naturally assumed to be in radians.

NOTE: Before opening an issue reporting something "wrong" with this function,
be sure to read all of [this
page](https://github.com/moble/quaternion/wiki/Euler-angles-are-horrible),
*especially* the very last section about opening issues or pull requests.

# Returns
- `??????::Vector{T}`

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
    a0 = 2acos(???((q.w^2+q.z^2)/abs2(q)))
    a1 = atan(q.z, q.w)
    a2 = atan(-q.x, q.y)
    @SVector [a1+a2, a0, a1-a2]
end


"""
    from_euler_angles(??, ??, ??)

Come over from the dark side.

Assumes the Euler angles correspond to the quaternion `R` via

    R = exp(??????/2) * exp(??????/2) * exp(??????/2)

where ???? and ???? rotate about the fixed ``y`` and ``z`` axes, respectively, so
this reprents an initial rotation (in the positive sense) through an angle
``??`` about the axis ``z``, followed by a rotation through ``??`` about the axis
``y``, and a final rotation through ``??`` about the axis ``z``.  This is
equivalent to performing an initial rotation through ``??`` about the axis
``z``, followed by a rotation through ``??`` about the *rotated* axis ``y'``,
followed by a rotation through ``??`` about the *twice-rotated* axis ``z''``.
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
function from_euler_angles(??, ??, ??)
    Rotor(
        cos(??/2)*cos((??+??)/2),
        -sin(??/2)*sin((??-??)/2),
        sin(??/2)*cos((??-??)/2),
        cos(??/2)*sin((??+??)/2)
    )
end
from_euler_angles(??????) = from_euler_angles(??????...)


"""
    to_euler_phases(q)

Convert input quaternion to complex phases of Euler angles

Interpreting the input quaternion as a rotation (though its normalization
scales out), we can define the complex Euler phases from the Euler angles (??,
??, ??) as

    z??? ??? exp(i*??)
    z??? ??? exp(i*??)
    z??? ??? exp(i*??)

These are more useful geometric quantites than the angles themselves ??? being
involved in computing spherical harmonics and Wigner's ???? matrices ??? and can be
computed from the components of the corresponding quaternion algebraically
(without the use of transcendental functions).

Note that `to_euler_phases!(z, q)` is supported for backwards compatibility,
but because this function returns an `SVector`, there is probably no advantage
to the in-place approach.

# Returns
- `z::SVector{Complex{T}}`: complex phases (z???, z???, z???) in that order.

# See Also
- [`from_euler_phases`](@ref): Create quaternion from Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function to_euler_phases(R::AbstractQuaternion{T}) where {T}
    a = R[1]^2 + R[4]^2
    b = R[2]^2 + R[3]^2
    sqrta = ???a
    sqrtb = ???b
    if iszero(sqrta)
        zp = one(Complex{T})
    else
        zp = Complex{T}(R[1], R[4]) / sqrta  # exp[i(??+??)/2]
    end
    if iszero(sqrtb)
        zm = one(Complex{T})
    else
        zm = Complex{T}(R[3], -R[2]) / sqrtb  # exp[i(??-??)/2]
    end
    @SVector [
        zp * zm,  # exp[i??]
        Complex{T}((a - b), 2 * sqrta * sqrtb) / (a + b),  # exp[i??]
        zp * conj(zm),  # exp[i??]
    ]
end


function to_euler_phases!(z::Array{Complex{T}}, R::AbstractQuaternion{T}) where {T}
    z[:] = to_euler_phases(R)
    z
end


"""
    from_euler_phases(z???, z???, z???)
    from_euler_phases(z)

Return the quaternion corresponding to these Euler phases.

Interpreting the input quaternion as a rotation (though its normalization
scales out), we can define the complex Euler phases from the Euler angles (??,
??, ??) as

    z??? ??? exp(i*??)
    z??? ??? exp(i*??)
    z??? ??? exp(i*??)

These are more useful geometric quantites than the angles themselves ??? being
involved in computing spherical harmonics and Wigner's ???? matrices ??? and can be
computed from the components of the corresponding quaternion algebraically
(without the use of transcendental functions).

# Parameters

- `z::Vector{Complex{T}}`: complex vector of length 3, representing the complex
  phases (z???, z???, z???) in that order.

# Returns
- `R::Quaternion{T}`

# See Also
- [`to_euler_phases`](@ref): Convert quaternion to Euler phases
- [`to_euler_angles`](@ref): Convert quaternion to Euler angles
- [`from_euler_angles`](@ref): Create quaternion from Euler angles

"""
function from_euler_phases(z???, z???, z???)
    zb = ???(z???)  # exp[i??/2]
    zp = ???(z??? * z???)  # exp[i(??+??)/2]
    zm = ???(z??? * conj(z???))  # exp[i(??-??)/2]
    if abs(z??? - zp * zm) > abs(z??? + zp * zm)
        zp *= -1
    end
    Rotor(zb.re * zp.re, -zb.im * zm.im, zb.im * zm.re, zb.re * zp.im)
end
from_euler_phases(z) = from_euler_phases(z...)


"""
    to_spherical_coordinates(q)

Return the spherical coordinates corresponding to this quaternion.

We can treat the quaternion as a transformation taking the ``z`` axis to some
direction ``n??``.  This direction can be described in terms of spherical
coordinates (??, ??).  Here, we use the standard commonly used in physics: ??
represents the "polar angle" between the ``z`` axis and the direction ``n??``,
while ?? represents the "azimuthal angle" between the ``x`` axis and the
projection of ``n??`` into the ``x``-``y`` plane.  Both angles are given in
radians.

"""
function to_spherical_coordinates(q::Q) where {Q<:AbstractQuaternion}
    q = float(q)
    a0 = 2acos(???((q.w^2+q.z^2)/abs2(q)))
    a1 = atan(q.z, q.w)
    a2 = atan(-q.x, q.y)
    @SVector [a0, a1+a2]
end


"""
    from_spherical_coordinates(??, ??)

Return a rotor corresponding to these spherical coordinates.

Considering (??, ??) as a point ``n??`` on the sphere, we can also construct a
quaternion that rotates the ``z`` axis onto that point.  Here, we use the
standard commonly used in physics: ?? represents the "polar angle" between the
``z`` axis and the direction ``n??``, while ?? represents the "azimuthal angle"
between the ``x`` axis and the projection of ``n??`` into the ``x``-``y`` plane.
Both angles must be given in radians.

"""
function from_spherical_coordinates(??, ??)
    s??, c?? = sincos(??/2)
    s??, c?? = sincos(??/2)
    Rotor(c??*c??, -s??*s??, s??*c??, c??*s??)
end
from_spherical_coordinates(????) = from_spherical_coordinates(????...)


"""
    from_rotation_matrix(???)

Convert 3x3 rotation matrix to quaternion.

Assuming the 3x3 matrix `???` rotates a vector `v` according to

    v' = ??? * v,

we can also express this rotation in terms of a quaternion `R` such that

    v' = R * v * R?????.

This function returns that quaternion, using Bar-Itzhack's algorithm to allow
for non-orthogonal matrices.  [J. Guidance, Vol. 23, No. 6,
p. 1085](http://dx.doi.org/10.2514/2.4654)

"""
function from_rotation_matrix(???::AbstractMatrix)
    @assert size(???) == (3, 3)
    @inbounds begin
        K = Array{eltype(???), 2}(undef, (4, 4))
        K[1, 1] = (???[1, 1] - ???[2, 2] - ???[3, 3])/3
        K[1, 2] = (???[2, 1] + ???[1, 2])/3
        K[1, 3] = (???[3, 1] + ???[1, 3])/3
        K[1, 4] = (???[2, 3] - ???[3, 2])/3
        K[2, 2] = (???[2, 2] - ???[1, 1] - ???[3, 3])/3
        K[2, 3] = (???[3, 2] + ???[2, 3])/3
        K[2, 4] = (???[3, 1] - ???[1, 3])/3
        K[3, 3] = (???[3, 3] - ???[1, 1] - ???[2, 2])/3
        K[3, 4] = (???[1, 2] - ???[2, 1])/3
        K[4, 4] = (???[1, 1] + ???[2, 2] + ???[3, 3])/3
        H = Symmetric(K)

        # compute the *dominant* (largest eigenvalue) eigenvector
        eigenvec = eigen(H, 4:4).vectors[:, 1]

        # convert it into a quaternion
        R = Rotor(eigenvec[4], -eigenvec[1], -eigenvec[2], -eigenvec[3])
    end
    R
end


"""
    to_rotation_matrix(q)

Convert quaternion to 3x3 rotation matrix.

Assuming the quaternion `R` rotates a vector `v` according to

    v' = R * v * R?????,

we can also express this rotation in terms of a 3x3 matrix `???` such that

    v' = ??? * v.

This function returns that matrix.

"""
function to_rotation_matrix(q::Q) where {Q<:AbstractQuaternion}
    n = inv(abs2(q))
    @SMatrix [
        1 - 2*(q.y^2 + q.z^2) * n  2*(q.x*q.y - q.z*q.w) * n  2*(q.x*q.z + q.y*q.w) * n ;
        2*(q.x*q.y + q.z*q.w) * n  1 - 2*(q.x^2 + q.z^2) * n  2*(q.y*q.z - q.x*q.w) * n ;
        2*(q.x*q.z - q.y*q.w) * n  2*(q.y*q.z + q.x*q.w) * n  1 - 2*(q.x^2 + q.y^2) * n
    ]
end
