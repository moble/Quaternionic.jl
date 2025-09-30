module QuaternionicChainRulesExt

using Quaternionic
import Quaternionic: _sincu, _cossu
using StaticArrays
isdefined(Base, :get_extension) ?
    (using ChainRules; import ChainRulesCore: ChainRulesCore, rrule, rrule_via_ad, RuleConfig, ProjectTo) :
    (using ..ChainRules; import ...ChainRulesCore: ChainRulesCore, rrule, rrule_via_ad, RuleConfig, ProjectTo)

if length(methods(ChainRulesCore.rrule, (typeof(hypot), Real, Real, Real ))) > 0

    function frule(
        (_, Δx, Δy, Δz, Δxs...),
        ::typeof(hypot),
        x::Union{Real,Complex},
        y::Union{Real,Complex},
        z::Union{Real,Complex},
        xs::Union{Real,Complex}...,
    )
        Ω = hypot(x, y, z, xs...)
        n = ifelse(iszero(Ω), oneunit(Ω), Ω)
        ∂Ω = sum(map(realdot, (x, y, z, xs...), (Δx, Δy, Δz, Δxs...))) / n
        return Ω, ∂Ω
    end

    function rrule(
        ::typeof(hypot),
        x::Union{Real,Complex},
        y::Union{Real,Complex},
        z::Union{Real,Complex},
        xs::Union{Real,Complex}...,
    )
        Ω = hypot(x, y, z, xs...)
        function hypot_pullback(ΔΩ)
            c = real(ΔΩ) / ifelse(iszero(Ω), one(Ω), Ω)
            return (NoTangent(), c * x, c * y, c * z, map(xi -> c * xi, xs)...)
        end
        return (Ω, hypot_pullback)
    end

end

end
