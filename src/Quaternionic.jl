module Quaternionic

using StaticArrays, Latexify, LaTeXStrings, LinearAlgebra, Requires
import Random: AbstractRNG, default_rng, randn!
import Symbolics

export AbstractQuaternion
export Quaternion, QuaternionF64, QuaternionF32, QuaternionF16, imx, imy, imz, 𝐢, 𝐣, 𝐤
export Rotor, RotorF64, RotorF32, RotorF16
export QuatVec, QuatVecF64, QuatVecF32, QuatVecF16
export (⋅), (×), (×̂), normalize
export abs2vec, absvec
export from_float_array, to_float_array,
    from_euler_angles, to_euler_angles,
    from_euler_phases, to_euler_phases!, to_euler_phases,
    from_spherical_coordinates, to_spherical_coordinates,
    from_rotation_matrix, to_rotation_matrix
export distance, distance2
export unflip, unflip!, slerp, squad
export ∂log, log∂log, ∂exp, exp∂exp, slerp∂slerp, slerp∂slerp∂τ, squad∂squad∂t
export precessing_nutating_example

abstract type AbstractQuaternion{T<:Real} <: Number end


include("quaternion.jl")
include("base.jl")
include("algebra.jl")
include("math.jl")
include("random.jl")
include("conversion.jl")
include("distance.jl")
include("interpolation.jl")
include("gradients.jl")
include("examples.jl")


function __init__()
    @require ForwardDiff="f6369f11-7733-5829-9624-2563aa707210" begin
        # Let ForwardDiff act naturally on quaternions.
        @inline function ForwardDiff.extract_derivative(::Type{T}, y::AbstractQuaternion) where {T}
            Quaternion(ForwardDiff.extract_derivative(T, y.components))
        end
    end

    # @require DifferentialEquations="0c46a032-eb83-5123-abaf-570d42b7fbaa" begin
    # end
end

end  # module
