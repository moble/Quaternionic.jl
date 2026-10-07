"""
    randn([rng=default_rng()], QT, [dims...])

Generate a normally distributed random quaternion of type `QT`, or an *array* of such
quaternions with dimensions `dims`.  The type `QT` must specify the element type, which can
be any `AbstractFloat`, as in `QuaternionF64`, `Rotor{Float32}`, or `QuatVec{BigFloat}`.
(There is no default type; `randn()` with no type argument still returns a `Float64`.)

If `QT` is a `Quaternion` type, the values are drawn from the spherically symmetric
quaternionic normal distribution with mean 0 and variance 1, so that the expected value of
`abs2(q)` is 1.  This corresponds to each of the four components having an independent
normal distribution with mean 0 and variance 1/4.

If `QT` is a `Rotor` type, the result is normalized.  Because the distribution is
spherically symmetric, the result is a uniformly distributed random rotation.

If `QT` is a `QuatVec` type, the result has zero scalar part, and its vector part has mean 0
and variance 1, corresponding to each of the three vector components having an independent
normal distribution with mean 0 and variance 1/3.

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
Base.randn(rng::AbstractRNG, QT::Type{<:AbstractQuaternion{T}}) where {T<:AbstractFloat} =
    QT(randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2)

Base.randn(rng::AbstractRNG, QT::Type{<:Rotor{T}}) where {T<:AbstractFloat} =
    rotor(randn(rng, T), randn(rng, T), randn(rng, T), randn(rng, T))

function Base.randn(rng::AbstractRNG, QT::Type{QuatVec{T}}) where {T<:AbstractFloat}
    SQRT_ONE_THIRD = sqrt(inv(T(3)))
    QT(
        zero(T),
        SQRT_ONE_THIRD * randn(rng, T),
        SQRT_ONE_THIRD * randn(rng, T),
        SQRT_ONE_THIRD * randn(rng, T)
    )
end
