"""
    precessing_nutating_example([Ωₒ], [Ωₚ], [α], [α̇], [ν], [R₀])

Return an exact quaternionic function of time representing nutating
precessional orbital motion, as described in Sec. 6.2 of [this
paper](https://arxiv.org/abs/1604.08139).  This example is useful because it
provides physically realistic but *very* complicated motion, while still being
simple to code up and differentiate analytically.

The output represents the rotation of the line joining the orbiting bodies and
their angular velocity.  The following are the input parameters, along with
their default values and physical interpretations:

  * `Ωₒ=2π/1_000`: The orbital frequency
  * `Ωₚ=2π/10_000`: The precessional frequency
  * `α=π/8`: The opening angle of the precession cone at ``t=0``
  * `α̇=2α/100_000`: The rate of opening of the precession cone
  * `ν=π/80`: The angle of nutation
  * `R₀=exp(-3α*imx/10)`: Overal rotation of the system

The default values are chosen to be typical of a (potentially) real precessing
binary black hole system shortly before merger.

The returned objects are three functions of time: `R`, `ω⃗`, and `Ṙ`, which
return the orientation as a `Rotor`, followed by the angular velocity as a
`QuatVec`, and the time-derivative of `R` as a `Quaternion`.


# Example
```jldoctest
julia> R, ω⃗, Ṙ = precessing_nutating_example();

julia> R(12.34)
rotor(0.9944579779058746 + 0.09804177421238346𝐢 - 0.0008485045352531196𝐣 + 0.03795287510453948𝐤)
julia> ω⃗(345.67)
 + 0.00046343007342866996𝐢 - 0.0007032818419003173𝐣 + 0.006214814810035087𝐤
julia> ϵ = 1e-6; (R(ϵ) - R(-ϵ)) / 2ϵ  # Approximate derivative at t=0
-3.8491432263754177e-7 + (3.9080960689830135e-6)𝐢 - (6.861695854245622e-5)𝐣 + 0.003076329202503836𝐤
julia> Ṙ(0)
-3.849124099815341e-7 + (3.9080812828476746e-6)𝐢 - (6.861695854245653e-5)𝐣 + 0.003076329202503835𝐤
```

"""
function precessing_nutating_example(Ωₒ=2π/1_000, Ωₚ=2π/10_000, α=π/8, α̇=2α/100_000, ν=π/80, R₀=exp(-3α*imx/10))
    R₁(t) = exp(Ωₒ*t*imz/2)
    R₂(t) = exp((α + α̇*t)*imx/2)
    R₃(t) = exp(Ωₚ*t*imz/2)
    R₄(t) = exp(ν*imx/2)
    R(t) = R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * R₁(t)
    Ṙ₀ = 0imz
    Ṙ₁(t) = R₁(t) * (Ωₒ*imz/2)
    Ṙ₂(t) = R₂(t) * (α̇*imx/2)
    Ṙ₃(t) = R₃(t) * (Ωₚ*imz/2)
    Ṙ₄(t) = 0imz
    Ṙ(t) = (
        Ṙ₀ * R₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * R₁(t) +
        R₀ * Ṙ₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * R₁(t) +
        R₀ * R₁(t) * Ṙ₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * R₁(t) +
        -R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * Ṙ₁(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * R₁(t) +
        R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * Ṙ₃(t) * R₂(t) * inv(R₃(t)) * R₁(t) +
        R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * Ṙ₂(t) * inv(R₃(t)) * R₁(t) +
        -R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * Ṙ₃(t) * inv(R₃(t)) * R₁(t) +
        R₀ * R₁(t) * R₄(t) * inv(R₁(t)) * R₃(t) * R₂(t) * inv(R₃(t)) * Ṙ₁(t)
    )
    ω⃗(t) = 2 * quatvec(Ṙ(t) * inv(R(t)))
    (R, ω⃗, Ṙ)
end
