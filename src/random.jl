"""
    randn([rng=GLOBAL_RNG], T=Quaternion{Float64}, [dims...])

Generate a normally distributed random quaternion of type `T` with mean 0 and standard
deviation of norm 1.  Optionally generate an *array* of such quaternions.  This module
currently provides an implementation for the types `QuaternionF16`, `QuaternionF32`, and
`QuaternionF64` (the default).  The values are drawn from the spherically symmetric
quaternionic normal distribution of variance 1 (corresponding to each component having
independent normal distribution with mean zero and variance 1/4).

Note that this function works with `Quaternion{BigFloat}`, even though `Base.randn` does
not work with `BigFloat`; we just use the [Box-Muller
transform](https://en.wikipedia.org/wiki/BoxβMuller_transform) to obtain the desired
result.

If the quaternion type passed in is `Rotor`, the result will be normalized correctly.
Because the distribution is spherically symmetric, the result is a truly random
rotation.

If the quaternion type is `QuatVec`, the result will have a 0 scalar component, and the
vector will have mean 0 standard deviation of norm 1.

# Examples
```julia
julia> randn(QuaternionF64)
0.4336736009756228 - 0.45087190792840853π’ - 0.24723937675211696π£ - 0.4514571469326208π€
julia> randn(QuaternionF16, 2, 2)
2Γ2 Matrix{QuaternionF16}:
   0.4321 + 1.105π’ + 0.2664π£ - 0.1359π€   0.064 + 0.9263π’ - 0.4138π£ + 0.05505π€
 0.2512 - 0.2585π’ - 0.2803π£ - 0.00964π€  -0.1256 + 0.1848π’ + 0.03607π£ - 0.752π€
```

"""
Base.randn(rng::AbstractRNG, QT::Type{<:AbstractQuaternion{T}}) where {T<:AbstractFloat} =
    QT(randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2, randn(rng, T)/2)

function Base.randn(rng::AbstractRNG, QT::Type{<:AbstractQuaternion{BigFloat}})
    # Use the Box-Muller transform to get randn BigFloats from rand BigFloat
    c = rand(rng, BigFloat, 4)
    l1 = β(-log(c[1])/2)
    s2, c2 = sincospi(2*c[2])
    l3 = β(-log(c[3])/2)
    s4, c4 = sincospi(2*c[4])
    QT(l1*c2, l1*s2, l3*c4, l3*s4)
end

_q3v_factor(::Type{T}) where T = inv(βT(3))

Base.randn(rng::AbstractRNG, QT::Type{QuatVec{T}}) where {T<:AbstractFloat} =
    QT(0, randn(rng, T)*_q3v_factor(T), randn(rng, T)*_q3v_factor(T), randn(rng, T)*_q3v_factor(T))

function Base.randn(rng::AbstractRNG, QT::Type{QuatVec{BigFloat}})
    # Use the Box-Muller transform to get randn BigFloats from rand BigFloat
    c = rand(rng, BigFloat, 4)
    l1 = β(-log(c[1])/2)
    s2 = sinpi(2*c[2])
    l3 = β(-log(c[3])/2)
    s4, c4 = sincospi(2*c[4])
    QuatVec(0, l1*s2, l3*c4, l3*s4)
end
