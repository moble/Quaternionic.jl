module Quaternionic

using StaticArrays, LinearAlgebra, PrecompileTools
import LaTeXStrings
import Random: AbstractRNG, default_rng, randn!
using TestItems: @testitem

export AbstractQuaternion
export Quaternion, quaternion,
    QuaternionF64, QuaternionF32, QuaternionF16,
    imx, imy, imz, 𝐢, 𝐣, 𝐤
export Rotor, rotor, RotorF64, RotorF32, RotorF16
export QuatVec, quatvec, QuatVecF64, QuatVecF32, QuatVecF16
export components, basetype
export (⋅), (×), (×̂), normalize
export abs2vec, absvec
export from_float_array, to_float_array,
    from_euler_angles, to_euler_angles,
    from_euler_phases, to_euler_phases!, to_euler_phases,
    from_spherical_coordinates, to_spherical_coordinates,
    from_rotation_matrix, to_rotation_matrix
export distance, distance2
export align
export unflip, unflip!, slerp, squad
export ∂log, log∂log, ∂exp, exp∂exp, slerp∂slerp, slerp∂slerp∂τ, squad∂squad∂t
export precessing_nutating_example

abstract type AbstractQuaternion{T<:Number} <: Number end


include("utilities.jl")
include("quaternion.jl")
include("base.jl")
include("algebra.jl")
include("math.jl")
include("gradients_exp_log.jl")
include("random.jl")
include("conversion.jl")
include("distance.jl")
include("alignment.jl")
include("interpolation.jl")
include("gradients_interpolation.jl")
include("examples.jl")

include("precompilation.jl")

# This symbol is only defined on Julia versions that support extensions
if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    # COV_EXCL_START

    function __init__()
        @require ChainRules="082447d4-558c-5d27-93f4-14fc19e9eca2" include("../ext/QuaternionicChainRulesExt.jl")
        @require ChainRulesCore="d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4" include("../ext/QuaternionicChainRulesCoreExt.jl")
        @require FastDifferentiation="eb9bf01b-bf85-4b60-bf87-ee5de06c00be" include("../ext/QuaternionicFastDifferentiationExt.jl")
        @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" include("../ext/QuaternionicForwardDiffExt.jl")
        @require Symbolics="0c5d862f-8b57-4792-8d23-62f2024744c7" include("../ext/QuaternionicSymbolicsExt.jl")
        @require Latexify="23fbe1c1-3f47-55db-b15f-69d7ec21a316" include("../ext/QuaternionicLatexifyExt.jl")
    end

    # COV_EXCL_STOP
end

end  # module
