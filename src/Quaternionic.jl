module Quaternionic

using StaticArrays, LinearAlgebra, PrecompileTools
using Latexify, LaTeXStrings
import Random: AbstractRNG, default_rng, randn!
import Symbolics

export AbstractQuaternion
export Quaternion, QuaternionF64, QuaternionF32, QuaternionF16, imx, imy, imz, 𝐢, 𝐣, 𝐤
export Rotor, RotorF64, RotorF32, RotorF16
export QuatVec, QuatVecF64, QuatVecF32, QuatVecF16
export components
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


include("quaternion.jl")
include("base.jl")
include("algebra.jl")
include("math.jl")
include("random.jl")
include("conversion.jl")
include("distance.jl")
include("alignment.jl")
include("interpolation.jl")
include("gradients.jl")
include("examples.jl")

include("precompilation.jl")

# This symbol is only defined on Julia versions that support extensions
if !isdefined(Base, :get_extension)
    using Requires
end

@static if !isdefined(Base, :get_extension)
    function __init__()
        @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" include("../ext/QuaternionicForwardDiffExt.jl")
    end
end

end  # module
