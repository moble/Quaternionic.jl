"""
    distance(q₁, q₂; multiplicative=false, sqrt=false, rotation=false)

Measure the "distance" between two quaternions.

By default, this function just returns the natural measure in the *additive*
group of quaternions:
```julia
abs2(q₁ - q₂)
```
If `multiplicative=true` is passed as a keyword, the function returns the
natural measure in the *multiplicative* group:
```julia
abs2(log(q₁ / q₂))
```
[Note that this will return `NaN` if either input is 0.]

If `sqrt=true` is passed, the square-root of the result will be taken (so,
`abs` instead of `abs2`).

Finally, if `rotation=true` is passed, the input quaternions will be
interpreted as quaternions, in which case the result will be the smallest
possible value for any combination of their signs.

See also [`distance_rotation`](@ref).

# Examples
```jldoctest example
julia> distance(imz, -imz)
4
julia> distance(imz, -imz, rotation=true)
0
julia> distance(1, imx, multiplicative=true, sqrt=true)  # π/2
1.5707963267948966
```

"""
@inline function distance(q₁, q₂; multiplicative=false, sqrt=false, rotation=false)
    if multiplicative
        if sqrt
            if rotation
                return distance_multiplicative_sqrt_rotation(q₁, q₂)
            else
                return distance_multiplicative_sqrt(q₁, q₂)
            end
        else
            if rotation
                return distance_multiplicative_rotation(q₁, q₂)
            else
                return distance_multiplicative(q₁, q₂)
            end
        end
    else
        if sqrt
            if rotation
                return distance_additive_sqrt_rotation(q₁, q₂)
            else
                return distance_additive_sqrt(q₁, q₂)
            end
        else
            if rotation
                return distance_additive_rotation(q₁, q₂)
            else
                return distance_additive(q₁, q₂)
            end
        end
    end
end


"""
    distance_rotation(q₁, q₂; sqrt=false)

Return `abs2(log(q₁/q₂))`, but assume that the signs and magnitudes of the
input quaternions do not matter.  If `sqrt=true`, we return the square-root of
that result.

This function is just a simple wrapper calling the [`distance`](@ref) function
as
```julia
distance(q₁, q₂; multiplicative=true, sqrt=sqrt, rotation=true)
```
"""
@inline distance_rotation(q₁, q₂; sqrt=false) = distance(q₁, q₂; multiplicative=true, sqrt=sqrt, rotation=true)


function _abs2_small_log_vec(q::Quaternion{T}) where {T}
    # Like `abs2(log(q))`, but assume magnitude is 1, and don't bother returning -π𝐤 for -1
    absolutevec = absvec(q)
    if absolutevec ≤ eps(q.w)
        return zero(float(T))
    end
    atan(absolutevec, q.w)^2
end


@inline distance_multiplicative(q₁, q₂) = abs2(log(q₁ / q₂))
@inline distance_multiplicative_rotation(q₁, q₂) = _abs2_small_log_vec(q₁ / q₂)
@inline distance_multiplicative_sqrt(q₁, q₂) = √distance_multiplicative(q₁, q₂)
@inline distance_multiplicative_sqrt_rotation(q₁, q₂) = √distance_multiplicative_rotation(q₁, q₂)
@inline distance_additive(q₁, q₂) = abs2(q₁ - q₂)
@inline distance_additive_rotation(q₁, q₂) = min(abs2(q₁ - q₂), abs2(q₁ + q₂))
@inline distance_additive_sqrt(q₁, q₂) = √distance_additive(q₁, q₂)
@inline distance_additive_sqrt_rotation(q₁, q₂) = √distance_additive_rotation(q₁, q₂)
