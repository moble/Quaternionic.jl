"""
    randn([rng=GLOBAL_RNG], [T=Quaternion{Float64}], [dims...])

Generate a normally distributed random quaternion of type `T` with mean 0 and standard
deviation of norm 1.  Optionally generate an *array* of such quaternions.  This module
currently provides an implementation for the types `QuaternionF16`, `QuaternionF32`, and
`QuaternionF64` (the default).  The values are drawn from the spherically symmetric
quaternionic normal distribution of variance 1 (corresponding to each component having
independent normal distribution with mean zero and variance 1/4).

See also: [`randn_rotor`](@ref)

# Examples
```julia
julia> randn(QuaternionF64)
0.4336736009756228 - 0.45087190792840853𝐢 - 0.24723937675211696𝐣 - 0.4514571469326208𝐤
julia> randn(QuaternionF16, 2, 2)
2×2 Matrix{QuaternionF16}:
   0.4321 + 1.105𝐢 + 0.2664𝐣 - 0.1359𝐤   0.064 + 0.9263𝐢 - 0.4138𝐣 + 0.05505𝐤
 0.2512 - 0.2585𝐢 - 0.2803𝐣 - 0.00964𝐤  -0.1256 + 0.1848𝐢 + 0.03607𝐣 - 0.752𝐤
```

"""
Base.randn(rng::AbstractRNG, ::Type{Quaternion{T}}) where {T<:AbstractFloat} =
    Quaternion{T}(randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2)

"""
    randn_rotor([rng=GLOBAL_RNG], [T=Quaternion{Float64}], [dims...])

Generate a normally distributed random quaternion of type `T` with mean 0 and norm 1.  (Note that
the *norm* is always precisely 1 with this function, but otherwise the individual components are
randomly distributed.)  The result is spherically symmetric, and gives rise a truly random rotation.

See also: [`randn`](@ref)
"""
function randn_rotor(rng::AbstractRNG, ::Type{T}, dims::Dims) where {T<:AbstractFloat}
    q = randn(rng, Quaternion{T}, dims)
    @. q / abs(q)
end
randn_rotor(rng::AbstractRNG, ::Type{Quaternion{T}}, dims::Dims) where {T<:AbstractFloat} =
    randn_rotor(rng, T, dims)
# Note: The following are more-or-less as given in the original `randn` definition
function randn_rotor(rng::AbstractRNG, ::Type{T}, dim1::Integer, dims::Integer...) where {T<:AbstractFloat}
    q = randn(rng, Quaternion{T}, dim1, dims...)
    @. q / abs(q)
end
randn_rotor(rng::AbstractRNG, ::Type{Quaternion{T}}, dim1::Integer, dims::Integer...) where {T<:AbstractFloat} =
    randn_rotor(rng, T, dim1, dims...)
randn_rotor(rng::AbstractRNG, ::Type{Quaternion{T}}                                 ) where {T<:AbstractFloat} =
    randn_rotor(rng, T, ())
randn_rotor(rng::AbstractRNG, ::Type{T}                                             ) where {T<:AbstractFloat} =
    randn_rotor(rng, T, ())
randn_rotor(                  ::Type{Quaternion{T}}, dims::Dims                     ) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, dims)
randn_rotor(                  ::Type{T},             dims::Dims                     ) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, dims)
randn_rotor(                  ::Type{Quaternion{T}}, dim1::Integer, dims::Integer...) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, dim1, dims...)
randn_rotor(                  ::Type{T},             dim1::Integer, dims::Integer...) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, dim1, dims...)
randn_rotor(                  ::Type{Quaternion{T}}                                 ) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, ())
randn_rotor(                  ::Type{T},                                            ) where {T<:AbstractFloat} =
    randn_rotor(default_rng(), T, ())
randn_rotor(rng::AbstractRNG,                        dims::Dims                     ) =
    randn_rotor(rng, QuaternionF64, dims)
randn_rotor(rng::AbstractRNG,                        dims::Integer...               ) =
    randn_rotor(rng, QuaternionF64, dims...)
randn_rotor(                                         dims::Dims                     ) =
    randn_rotor(default_rng(), QuaternionF64, dims)
randn_rotor(                                         dim1::Integer, dims::Integer...) =
    randn_rotor(default_rng(), QuaternionF64, dim1, dims...)
randn_rotor(                                                                        ) =
    randn_rotor(default_rng(), QuaternionF64, ())
