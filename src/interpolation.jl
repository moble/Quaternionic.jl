# This requires a little more understanding of acting along certain axes of an
# array.  The only useful references I can find for that are
# https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes
# and
# https://discourse.julialang.org/t/crazy-allocations-using-cartesianindices/42262/2

const Rotator = Union{Quaternion, Rotor}


function _unflip!(q, Rpre, R, Rpost)
    @inbounds for Ipost in Rpost
        for i in R
            for Ipre in Rpre
                if q[Ipre, i-1, Ipost] ⋅ q[Ipre, i, Ipost] < 0
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
                if p[Ipre, i-1, Ipost] ⋅ q[Ipre, i, Ipost] < 0
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
4-element Vector{QuatVec{Int64}}:
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


"""
    slerp(q₁, q₂, τ; [unflip=false])

Spherical linear interpolation of quaternions

The result of a "slerp" is given by

        (q₂ / q₁)^τ * q₁

When `τ` is 0, this evaluates to `q₁`; when `τ` is 1, this evaluates to `q₂`;
for any other values the result varies between the two.

Note that applying this to successive pairs of quaternions as in `slerp(q₁, q₂,
τₐ)` and `slerp(q₂, q₃, τᵦ)` will be continuous, but the derivative will be
discontinuous when moving from the first pair to the second.  See
[`squad`](@ref) for a more continuous curve.

If `unflip=true` is passed as a keyword, and the input quaternions are more
anti-parallel than parallel, the sign of `q₂` will be flipped before the result
is computed.

"""
function slerp(q₁::R1, q₂::R2, τ::Real; unflip::Bool=false) where {R1<:Rotator, R2<:Rotator}
    if unflip && q₁⋅q₂ < 0
        return (-q₂ / q₁)^τ * q₁
    end
    (q₂ / q₁)^τ * q₁
end


@doc raw"""
    squad_control_points(R::AbstractVector{Rotor}, t::AbstractVector{<:AbstractFloat}, i::Int)

This is a helper function for the `squad` routines, returning the control
points between one pair of input rotors.

The expressions for ``A`` and ``B`` (assuming all indices are valid) are
```math
\begin{aligned}
A_{i} &= R_{i}\, \exp\left\{
  \frac{1}{4}
  \left[
    \log\left(\bar{R}_{i-1}\, R_i\right) \frac{t_{i+1} - t_{i}} {t_{i} - t_{i-1}}
    - \log\left(\bar{R}_{i}\, R_{i+1}\right)
  \right]
\right\},
\\
B_{i} &= R_{i+1}\, \exp\left\{
  -\frac{1}{4}
  \left[
    \log\left(\bar{R}_{i+1}\, R_{i+2}\right) \frac{t_{i+1} - t_{i}} {t_{i+2} - t_{i+1}}
    - \log\left(\bar{R}_{i}\, R_{i+1}\right)
  \right]
\right\}
\end{aligned}.
```
The indices will be invalid for `A[begin]` and `A[end]`, and for `B[end-1]` and
`B[end]`.  We can simply extend the input `R` values by linearly extrapolating,
which results in the following simplified results:
```math
\begin{aligned}
A_{\mathrm{begin}} &= R_{\mathrm{begin}} \\
A_{\mathrm{end}} &= R_{\mathrm{end}} \\
B_{\mathrm{end-1}} &= R_{\mathrm{end}} \\
B_{\mathrm{end}} &= R_{\mathrm{end}}\, \bar{R}_{\mathrm{end-1}}\, R_{\mathrm{end}} \\
                 &= 2\left(R_{\mathrm{end}}\cdot R_{\mathrm{end-1}}\right)\, R_{\mathrm{end}} - R_{\mathrm{end-1}}.
\end{aligned}
```
"""
function squad_control_points(R::AbstractVector{<:Rotor}, t::AbstractVector{<:AbstractFloat}, i::Int)
    if i==1
        A = R[1]
    elseif i==length(R)
        A = R[end]  # COV_EXCL_LINE
    else
        A = R[i] * exp(
            (
                log(conj(R[i-1]) * R[i]) * ((t[i+1] - t[i]) / (t[i] - t[i-1]))
                - log(conj(R[i]) * R[i+1])
            ) / 4
        )
    end
    if i<length(R)-1
        B = R[i+1] * exp(
            (
                log(conj(R[i+1]) * R[i+2]) * ((t[i+1] - t[i]) / (t[i+2] - t[i+1]))
                - log(conj(R[i]) * R[i+1])
            ) / -4
        )
    elseif i==length(R)-1
        B = R[i+1]
    else # i==length(R)
        B = Rotor{eltype(R)}((2*(R[i] ⋅ R[i-1]) * R[i] - R[i-1]).components...)  # COV_EXCL_LINE
    end
    A, B
end


const unflip_func = unflip  # `unflip` will be used as a local variable in
                            # `squad`, but the function will also be needed

"""
    squad(Rin, tin, tout; [unflip=false], [validate=false])

Perform "Spherical QUADrangle interpolation" on the input `Rotor`s `Rin` with
corresponding times `tin`, to the output times `tout`.

This is a slightly generalized version of [Shoemake's "spherical Bézier
curves"](https://doi.org/10.1145/325165.325242), to allow for time steps of
varying sizes.

The input `Rin` and `tin` must be vectors of the same length.  The output
`tout` may be either a single real number or a vector of real numbers.  Both
`tin` and `tout` are assumed to be sorted, and `tout` is assumed to be
contained entirely within `tin`; no extrapolation will be done.

If `unflip=true` is passed as a keyword, the [`unflip`](@ref) function will be
applied to `Rin`.

If `validate=true` is passed as a keyword, the time ordering of the input `tin`
and `tout` will be tested to ensure that no extrapolation will be done.

See also [`squad_control_points`](@ref).

"""
function squad(Rin::AbstractVector{<:Rotor}, tin::AbstractVector{<:Real}, tout::AbstractVector{<:Real}; unflip=false, validate=false)
    t_begin, t_end = extrema(tin)
    @assert t_begin < t_end  # Proves that there are at least 2 tin
    @assert length(Rin) == length(tin)  # Proves that there are at least 2 Rin
    if validate
        @assert minimum(diff(tin)) > 0
        if length(tout) > 1
            @assert minimum(diff(tout)) > 0
            tout_begin, tout_end = extrema(tout)
            @assert t_begin ≤ tout_begin
            @assert tout_end ≤ t_end
        elseif length(tout) == 1
            @assert t_begin ≤ tout[1] ≤ t_end
        end
    end
    if unflip
        Rin = unflip_func(Rin)
    end
    Rout = similar(Rin, length(tout))
    if length(tout) == 0
        return Rout
    end
    j = 1
    toutj = tout[j]
    if toutj == tin[1]
        Rout[1] = Rin[1]
        j += 1
    end
    while j ≤ length(Rout)
        toutj = tout[j]
        i = searchsortedfirst(tin, toutj) - 1
        if i == length(tin)
            error("Searching for $toutj went out of range [$t_begin, $t_end]")  # COV_EXCL_LINE
        end
        A, B = squad_control_points(Rin, tin, i)
        ta, tb = tin[i], tin[i+1]
        while j ≤ length(Rout) && tb ≥ tout[j]
            τ = (tout[j] - ta) / (tb - ta)
            Rout[j] = slerp(
                slerp(Rin[i], Rin[i+1], τ),
                slerp(A, B, τ),
                2τ*(1-τ)
            )
            j += 1
        end
    end
    Rout
end


function squad(Rin::AbstractVector{<:Rotor}, tin::AbstractVector{<:Real}, tout::Real; unflip=false, validate=false)
    squad(Rin, tin, [tout]; unflip=unflip, validate=validate)[1]
end
