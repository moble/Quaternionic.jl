# This requires a little more understanding of acting along certain axes of an
# array.  The only useful references I can find for that are
# https://julialang.org/blog/2016/02/iteration/#filtering_along_a_specified_dimension_exploiting_multiple_indexes
# and
# https://discourse.julialang.org/t/crazy-allocations-using-cartesianindices/42262/2

const Rotator = Union{Quaternion, Rotor}


# The sign test uses the real part of `componentdot`, so that complex quaternions are also
# supported.  For a Lorentz rotor, the real part of the scalar part of the relative
# transformation is cos(θ/2) cosh(η/2), which selects the hemisphere just as in the real
# case.

# In-place worker: `I` is the axis along the dimension being unflipped, and `Rpre` and
# `Rpost` are the Cartesian indices before and after that dimension.
function unflip!(q::AbstractArray, Rpre::CartesianIndices, I::AbstractUnitRange, Rpost::CartesianIndices)
    @inbounds for Ipost in Rpost
        for i in first(I)+1:last(I)
            for Ipre in Rpre
                if real(componentdot(q[Ipre, i-1, Ipost], q[Ipre, i, Ipost])) < 0
                    # Unary minus is exact and keeps the element type; multiplying a
                    # `Rotor` by -1 would give a `Quaternion`, which is renormalized when
                    # it is stored back.
                    q[Ipre, i, Ipost] = -q[Ipre, i, Ipost]
                end
            end
        end
    end
    q
end


# Copying worker: the unflipped values of `q` are written into `p`, which has the same axes.
function unflip!(p::AbstractArray, q::AbstractArray, Rpre::CartesianIndices, I::AbstractUnitRange, Rpost::CartesianIndices)
    isempty(I) && return p
    @inbounds for Ipost in Rpost
        let i = first(I)
            for Ipre in Rpre
                p[Ipre, i, Ipost] = q[Ipre, i, Ipost]
            end
        end
        for i in first(I)+1:last(I)
            for Ipre in Rpre
                if real(componentdot(p[Ipre, i-1, Ipost], q[Ipre, i, Ipost])) < 0
                    p[Ipre, i, Ipost] = -q[Ipre, i, Ipost]
                else
                    p[Ipre, i, Ipost] = q[Ipre, i, Ipost]
                end
            end
        end
    end
    p
end


"""
    unflip(q; dim=1)
    unflip!(q; dim=1)

Flip the signs of successive quaternions along dimension `dim` so that they are as
continuous as possible.

If `q` represents a series of rotations, the sign of each element is arbitrary.  However,
for certain purposes — such as interpolation and differentiation — the continuity of the
quaternions matters, and so we want the *quaternions* to be as continuous as possible
without changing the *rotations* that they represent.

The first element along `dim` is never changed.  Each subsequent element `q` is negated if
the real part of `p[1]*q[1] + p[2]*q[2] + p[3]*q[3] + p[4]*q[4]` is negative, where `p` is
the (possibly negated) preceding element.  `unflip` returns a new array, while `unflip!`
modifies `q` in place and returns it.

# Examples
```jldoctest
julia> q = [imx, -imx, imx, -imx];

julia> unflip(q)
4-element Vector{QuatVec{Int64}}:
  + 1𝐢 + 0𝐣 + 0𝐤
  + 1𝐢 + 0𝐣 + 0𝐤
  + 1𝐢 + 0𝐣 + 0𝐤
  + 1𝐢 + 0𝐣 + 0𝐤
```
"""
function unflip(q::AbstractArray{<:AbstractQuaternion}; dim::Integer=1)
    Rpre = CartesianIndices(axes(q)[1:dim-1])
    I = axes(q, dim)
    Rpost = CartesianIndices(axes(q)[dim+1:end])
    unflip!(similar(q), q, Rpre, I, Rpost)
end


"""
    unflip!(q; dim=1)

Flip the signs of successive quaternions along dimension `dim` of `q` in place, and return
`q`.  See [`unflip`](@ref) for details.
"""
function unflip!(q::AbstractArray{<:AbstractQuaternion}; dim::Integer=1)
    Rpre = CartesianIndices(axes(q)[1:dim-1])
    I = axes(q, dim)
    Rpost = CartesianIndices(axes(q)[dim+1:end])
    unflip!(q, Rpre, I, Rpost)
end


"""
    slerp(q₁, q₂, τ; [unflip=false])

"Spherical Linear intERPolation" of a pair of quaternions.

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
    if unflip && real(componentdot(q₁, q₂)) < 0
        return (-q₂ / q₁)^τ * q₁
    end
    (q₂ / q₁)^τ * q₁
end


@doc raw"""
    squad_control_points(R::AbstractVector{<:Rotor}, t::AbstractVector{<:Real}, i::Int)

This is a helper function for the `squad` routines, returning the control points between one
pair of input rotors.  Both control points are returned as `Rotor`s of the same type, which
is promoted from the types of `R` and `t`.

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
\right\}.
\end{aligned}
```
These expressions will be invalid for ``A_{\mathrm{begin}}``, ``A_{\mathrm{end}}``,
``B_{\mathrm{end-1}}``, and ``B_{\mathrm{end}}``, because they all involve out-of-bounds indices of
``R_i``.  We can simply extend the input `R` values by linearly extrapolating, which results in the
following simplified results:
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
function squad_control_points(R::AbstractVector{<:Rotor}, t::AbstractVector{<:Real}, i::Int)
    # Both control points are converted to one type, so that the return type does not depend
    # on `i`.
    RT = Rotor{promote_type(basetype(eltype(R)), typeof(float(one(eltype(t)))))}
    n = length(R)
    if i==1
        A = R[1]
    elseif i==n
        A = R[n]  # COV_EXCL_LINE
    else
        A = R[i] * exp(
            (
                log(conj(R[i-1]) * R[i]) * ((t[i+1] - t[i]) / (t[i] - t[i-1]))
                - log(conj(R[i]) * R[i+1])
            ) / 4
        )
    end
    if i<n-1
        B = R[i+1] * exp(
            (
                log(conj(R[i+1]) * R[i+2]) * ((t[i+1] - t[i]) / (t[i+2] - t[i+1]))
                - log(conj(R[i]) * R[i+1])
            ) / -4
        )
    elseif i==n-1
        B = R[n]
    else # i==n
        B = R[n] * conj(R[n-1]) * R[n]  # COV_EXCL_LINE
    end
    RT(A), RT(B)
end


# The value of `squad` at time `t` on the segment from `ta` to `tb`, and its derivative with
# respect to `t`.  With τ = (t - ta) / (tb - ta), the value is
#
#   s = slerp(X, Y, σ) = r^σ X,  where  X = slerp(qᵢ, qᵢ₊₁, τ),  Y = slerp(A, B, τ),
#   σ = 2τ(1 - τ),  and  r = Y / X.
#
# Since slerp(q₁, q₂, τ) = exp(τ L) q₁ with the constant L = log(q₂ / q₁), the derivatives
# with respect to τ (written with a dot) are Ẋ = L₁ X and Ẏ = L₂ Y, where
# L₁ = log(qᵢ₊₁ / qᵢ) and L₂ = log(B / A).  Then ṙ = L₂ r - r L₁, and with l = log(r), the
# derivative of r^σ = exp(σ l) follows from the pushforwards of `log` and `exp`:
#
#   l̇ = log_pushforward(r, ṙ),   (r^σ)˙ = exp_pushforward(σ l, σ̇ l + σ l̇),
#
# so that ṡ = (r^σ)˙ X + r^σ Ẋ, and the derivative with respect to `t` is ṡ / (tb - ta).
# The value is computed exactly as `squad!` computes it when no derivatives are requested.
function squad_with_derivative(qᵢ, A, B, qᵢ₊₁, ta, tb, t)
    # `float` keeps rational times away from `^(::Rotor, ::Rational)`.
    τ = float((t - ta) / (tb - ta))
    X = slerp(qᵢ, qᵢ₊₁, τ)
    Y = slerp(A, B, τ)
    σ = 2τ*(1-τ)
    r = Y / X
    e = r^σ
    s = e * X
    L₁ = log(qᵢ₊₁ / qᵢ)
    L₂ = log(B / A)
    l = quaternion(log(r))
    l̇ = log_pushforward(quaternion(r), L₂ * r - r * L₁)
    ė = exp_pushforward(σ * l, (2 - 4τ) * l + σ * l̇)
    ṡ = ė * X + e * (L₁ * X)
    (s, ṡ / (tb - ta))
end


const unflip_func = unflip  # `unflip` will be used as a local variable in
                            # `squad`, but the function will also be needed

"""
    squad!(Rout, Ω⃗out, Ṙout, Rin, tin, tout; [unflip=false], [validate=false])
    squad!(Rout, Rin, tin, tout; [unflip=false], [validate=false])

In-place evaluation of "Spherical QUADrangle interpolation".  Note that this is intended
mostly as a utility function; [`squad`](@ref) is more user-friendly.  However, for
efficiency, this function may be preferable.

The output arrays `Rout`, `Ω⃗out`, and `Ṙout` will be modified in place, and must have the
same length as `tout`.  Their elements must be `Rotor`, `QuatVec`, and `Quaternion`,
respectively.  Optionally, either or both of `Ω⃗out` and `Ṙout` may be `nothing`, in which
case they will not be computed.  The second form computes only `Rout`.

The times `tout` must lie within the range of `tin`; otherwise an error is thrown.

See also [`squad`](@ref).

"""
function squad!(
    Rout::AbstractVector{<:Rotor}, Ω⃗out::Union{Nothing, AbstractVector{<:QuatVec}}, Ṙout::Union{Nothing, AbstractVector{<:Quaternion}},
    Rin::AbstractVector{<:Rotor}, tin::AbstractVector{<:Real}, tout::AbstractVector{<:Real}; unflip=false, validate=false
)
    if length(tout) == 0
        return
    end
    Base.require_one_based_indexing(Rout, Rin, tin, tout)
    # These checks of the input use `ArgumentError` rather than `@assert`, which may be
    # disabled at some optimization levels.
    t_begin, t_end = extrema(tin)
    t_begin < t_end ||  # Proves that there are at least 2 tin
        throw(ArgumentError("`tin` must contain at least two distinct times"))
    length(Rin) == length(tin) ||  # Proves that there are at least 2 Rin
        throw(ArgumentError("`Rin` and `tin` must have the same length"))
    length(Rout) == length(tout) ||
        throw(ArgumentError("`Rout` and `tout` must have the same length"))
    evaluate_Ω⃗ = (Ω⃗out !== nothing)
    evaluate_Ṙ = (Ṙout !== nothing)
    if evaluate_Ω⃗
        length(Ω⃗out) == length(tout) ||
            throw(ArgumentError("`Ω⃗out` and `tout` must have the same length"))
    end
    if evaluate_Ṙ
        length(Ṙout) == length(tout) ||
            throw(ArgumentError("`Ṙout` and `tout` must have the same length"))
    end
    if validate
        minimum(diff(tin)) > 0 ||
            throw(ArgumentError("`tin` must be strictly increasing"))
        if length(tout) > 1
            minimum(diff(tout)) > 0 ||
                throw(ArgumentError("`tout` must be strictly increasing"))
        end
        tout_begin, tout_end = extrema(tout)
        t_begin ≤ tout_begin && tout_end ≤ t_end ||
            throw(ArgumentError("`tout` must lie within [$t_begin, $t_end]"))
    end
    if unflip
        Rin = unflip_func(Rin)
    end
    n = length(tin)
    j = 1
    while j ≤ length(Rout)
        # The segment is located by value.  AD types such as `ForwardDiff.Dual` may also
        # compare their partials when their values are equal, which would push a time equal
        # to a knot into a neighboring segment, or out of range at either end.
        toutj = value(tout[j])
        i = searchsortedfirst(tin, toutj; by=value) - 1
        if i == 0 && toutj == value(tin[1])
            i = 1
        end
        if i < 1 || i ≥ n
            error("Searching for $(tout[j]) went out of range [$t_begin, $t_end]")
        end
        A, B = squad_control_points(Rin, tin, i)
        ta, tb = tin[i], tin[i+1]
        qᵢ, qᵢ₊₁ = Rin[i], Rin[i+1]
        # Both ends of the segment are checked, so that an out-of-order time goes back to
        # the segment search instead of being extrapolated from the current segment.
        while j ≤ length(Rout) && value(ta) ≤ value(tout[j]) ≤ value(tb)
            if evaluate_Ω⃗ || evaluate_Ṙ
                s, ∂s∂t = squad_with_derivative(qᵢ, A, B, qᵢ₊₁, ta, tb, tout[j])
                Rout[j] = s
                if evaluate_Ω⃗
                    Ω⃗out[j] = 2 * eltype(Ω⃗out)(∂s∂t / s)
                end
                if evaluate_Ṙ
                    Ṙout[j] = ∂s∂t
                end
            else
                # `float` keeps rational times away from `^(::Rotor, ::Rational)`.
                τ = float((tout[j] - ta) / (tb - ta))
                Rout[j] = slerp(
                    slerp(qᵢ, qᵢ₊₁, τ),
                    slerp(A, B, τ),
                    2τ*(1-τ)
                )
            end
            j += 1
        end
    end
end

function squad!(
    Rout::AbstractVector{<:Rotor}, Rin::AbstractVector{<:Rotor}, tin::AbstractVector{<:Real},
    tout::AbstractVector{<:Real}; unflip=false, validate=false
)
    squad!(Rout, nothing, nothing, Rin, tin, tout; unflip=unflip, validate=validate)
end


"""
    squad(Rin, tin, tout; [kwargs...])

"Spherical QUADrangle interpolation" of the input `Rotor`s `Rin` with corresponding times
`tin`, to the output times `tout`.

This is a slightly generalized version of [Shoemake's "spherical Bézier
curves"](https://doi.org/10.1145/325165.325242), to allow for time steps of varying sizes.

The input `Rin` and `tin` must be vectors of the same length.  The output `tout` may be
either a single real number or a vector of real numbers.  The times `tin` are assumed to be
strictly increasing, and `tout` must be contained entirely within the range of `tin`; no
extrapolation will be done.  Sorted `tout` is evaluated most efficiently, but the result is
correct for any order.

See also [`squad!`](@ref) for in-place versions of this function.

# Keyword arguments

If `unflip=true` is passed as a keyword, the [`unflip`](@ref) function will be applied to
`Rin`.

If `validate=true` is passed as a keyword, the time ordering of the input `tin` and `tout`
will be tested to ensure that no extrapolation will be done, and an `ArgumentError` will be
thrown otherwise.  Even without validation, an error is thrown when an element of `tout`
falls outside the range of `tin`.

If `compute_angular_velocity=true` is passed as a keyword, the return value will be a tuple.
The first element of the tuple will be a vector of `Rotor`s as before, but the second
element will be a vector of `QuatVec`s representing the angular velocity.

If `compute_derivative=true` is passed as a keyword, the return value will be a tuple.  The
first element of the tuple will be a vector of `Rotor`s as before, but the last element will
be a vector of `Quaternion`s representing the time-derivative of the rotors.  Note that if
`compute_angular_velocity=true`, this tuple will have three elements.

"""
@inline function squad(
        Rin::AbstractVector{Rotor{T}}, tin::AbstractVector{<:Real}, tout::Union{Real, AbstractVector{<:Real}};
        unflip=false, validate=false, compute_angular_velocity=false, compute_derivative=false
) where {T}
    # The flags are lifted into the type domain, so that the return type can be inferred
    # whenever they are constants, as they are with the default values.  Inlining this
    # method lets the compiler see those constants.
    squad(
        Rin, tin, tout, Val(compute_angular_velocity), Val(compute_derivative);
        unflip=unflip, validate=validate
    )
end

function squad(
        Rin::AbstractVector{Rotor{T}}, tin::AbstractVector{<:Real}, tout::AbstractVector{<:Real},
        ::Val{Ω⃗}, ::Val{Ṙ}; unflip=false, validate=false
) where {T, Ω⃗, Ṙ}
    Rout_eltype = promote_type(eltype(tin), eltype(tout), T)
    Rout = similar(Rin, Rotor{Rout_eltype}, length(tout))
    Ω⃗out = Ω⃗ ? similar(Rin, QuatVec{Rout_eltype}, length(tout)) : nothing
    Ṙout = Ṙ ? similar(Rin, Quaternion{Rout_eltype}, length(tout)) : nothing
    squad!(Rout, Ω⃗out, Ṙout, Rin, tin, tout; unflip=unflip, validate=validate)
    if Ω⃗ && Ṙ
        return (Rout, Ω⃗out, Ṙout)
    elseif Ω⃗
        return (Rout, Ω⃗out)
    elseif Ṙ
        return (Rout, Ṙout)
    end
    return Rout
end

function squad(
        Rin::AbstractVector{Rotor{T}}, tin::AbstractVector{<:Real}, tout::Real,
        Ω⃗::Val, Ṙ::Val; unflip=false, validate=false
) where {T}
    result = squad(Rin, tin, [tout], Ω⃗, Ṙ; unflip=unflip, validate=validate)
    result isa Tuple ? map(first, result) : result[1]
end
