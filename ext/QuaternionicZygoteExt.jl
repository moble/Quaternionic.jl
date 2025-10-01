module QuaternionicZygoteExt

import Quaternionic: Quaternionic, components, AbstractQuaternion
isdefined(Base, :get_extension) ?
    (import Zygote; import ZygoteRules) :
    (import ...Zygote; import ...ZygoteRules)

import StaticArraysCore: SVector
import ChainRulesCore: unthunk

# Lift any Δq form into a tuple cotangent
function Δtuple(Δq)
    Δ = unthunk(Δq)
    if Δ isa AbstractQuaternion
        return components(Δ)
    elseif Δ isa SVector{4}
        return Tuple(Δ)
    elseif Δ isa Tuple && length(Δ) == 4
        return Δ
    elseif Δ isa Real
        return (Δ, zero(Δ), zero(Δ), zero(Δ))
    else
        error("Unsupported cotangent type: $(typeof(Δ))")
    end
end

Zygote.@adjoint function Base.:+(t::Real, q::AbstractQuaternion)
    q′ = t + q
    function pullback(Δq′)
        Δw, Δx, Δy, Δz = Δtuple(Δq′)
        return (Δw, typeof(q)(Δw, Δx, Δy, Δz))
    end
    return q′, pullback
end

Zygote.@adjoint function Base.:+(q::AbstractQuaternion, t::Real)
    q′ = q + t
    function pullback(Δq′)
        Δw, Δx, Δy, Δz = Δtuple(Δq′)
        return (typeof(q)(Δw, Δx, Δy, Δz), Δw)
    end
    return q′, pullback
end

Zygote.@adjoint function Base.:-(t::Real, q::AbstractQuaternion)
    q′ = t - q
    function pullback(Δq′)
        Δw, Δx, Δy, Δz = Δtuple(Δq′)
        return (Δw, typeof(q)(-Δw, -Δx, -Δy, -Δz))
    end
    return q′, pullback
end

Zygote.@adjoint function Base.:-(q::AbstractQuaternion, t::Real)
    q′ = q - t
    function pullback(Δq′)
        Δw, Δx, Δy, Δz = Δtuple(Δq′)
        return (typeof(q)(Δw, Δx, Δy, Δz), -Δw)
    end
    return q′, pullback
end


end # module
