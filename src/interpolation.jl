# This requires a little more understanding of acting along certain axes of an
# array.  The only useful references I can find for that are
# https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes
# and
# https://discourse.julialang.org/t/crazy-allocations-using-cartesianindices/42262/2


@inline @fastmath function _inner_product(q1::Quaternion, q2::Quaternion)
    q1.w*q2.w + q1.x*q2.x + q1.y*q2.y + q1.z*q2.z
end


function _unflip!(q, Rpre, R, Rpost)
    @inbounds for Ipost in Rpost
        for i in R
            for Ipre in Rpre
                if _inner_product(q[Ipre, i-1, Ipost], q[Ipre, i, Ipost]) < 0
                    q[Ipre, i, Ipost] *= -1
                end
            end
        end
    end
end


function _unflip!(q, p, Rpre, R, Rpost)
    @inbounds for Ipost in Rpost
        let i = 1
            for Ipre in Rpre
                p[Ipre, i, Ipost] = q[Ipre, i, Ipost]
            end
        end
        for i in R
            for Ipre in Rpre
                if _inner_product(p[Ipre, i-1, Ipost], q[Ipre, i, Ipost]) < 0
                    p[Ipre, i, Ipost] = -q[Ipre, i, Ipost]
                else
                    p[Ipre, i, Ipost] = q[Ipre, i, Ipost]
                end
            end
        end
    end
end


function unflip!(q::AbstractArray{<:AbstractQuaternion}; dim::Integer=1)
    Rpre = CartesianIndices(size(q)[1:dim-1])
    R = 2:size(q, dim)
    Rpost = CartesianIndices(size(q)[dim+1:end])
    _unflip!(q, Rpre, R, Rpost)
end


"""
    unflip(q, [dim=1])
    unflip!(q, [dim=1])

Flip the signs of successive quaternions along dimension `dim` so that they are
as continuous as possible.

If `q` represents a series of rotations, the sign of each element is arbitrary.
However, for certain purposes — such as interpolation and differentiation — the
continuity of the quaternions matters, and so we want the *quaternions* to be
as continuous as possible without changing the *rotations* that they represent.

# Examples
```jldoctest
julia> q = [imx, -imx, imx, -imx];

julia> unflip(q)
4-element Vector{Quaternion{Int64}}:
 0 + 1𝐢 + 0𝐣 + 0𝐤
 0 + 1𝐢 + 0𝐣 + 0𝐤
 0 + 1𝐢 + 0𝐣 + 0𝐤
 0 + 1𝐢 + 0𝐣 + 0𝐤
```
"""
function unflip(q::AbstractArray{<:AbstractQuaternion}; dim::Integer=1)
    p = similar(q)
    Rpre = CartesianIndices(size(q)[1:dim-1])
    R = 2:size(q, dim)
    Rpost = CartesianIndices(size(q)[dim+1:end])
    _unflip!(q, p, Rpre, R, Rpost)
    p
end


# """
#     slerp(q₁, q₂, τ)

# Spherical linear interpolation of quaternions

# The result of a "slerp" is given by

#         (q₂ / q₁)^τ * q₁

# When `τ` is 0, this evaluates to `q₁`; when `τ` is 1, this evaluates to `q₂`;
# for any other values the result varies between the two.

# See also [`squad`](@ref)
# """
# function slerp(q₁::Rotor, q₂::Rotor, τ)
#     exp_s_log(τ, q₂ / q₁) * q₁
# end


# function slerp(qᵢ::Vector{Rotor}, tᵢ::Vector{<:AbstractFloat}, tₒ; unflip=false)

# end


# function squad(qᵢ::Vector{Rotor}, tᵢ::Vector{<:AbstractFloat}, tₒ; unflip=false)

# end

# function exp_s_log(s, q::Rotor)
#     q = float(q)
#     absolutevec = absvec(q)
#     if absolutevec ≤ eps(typeof(absolutevec))
#         if q.w < 0
#             # log(q) ≈ π𝐤
#             sin_πs, cos_πs = sincospi(oftype(absolutevec, s))
#             return Quaternion(cos_πs, 0, 0, sin_πs)
#         end
#         # log(q) ≈ 0
#         return one(q)
#     end
#     f1 = oftype(absolutevec, s) * atan(absolutevec, q.w)
#     sin_f1, cos_f1 = sincos(f1)
#     f2 = sin_f1 / absolutevec
#     Quaternion(cos_f1, f2*q.x, f2*q.y, f2*q.z)
# end
