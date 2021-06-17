module Quaternionic

using StaticArrays, Latexify, LaTeXStrings, LinearAlgebra
import Random: AbstractRNG, default_rng, randn!
import Symbolics

export AbstractQuaternion
export Quaternion, QuaternionF64, QuaternionF32, QuaternionF16, imx, imy, imz, 𝐢, 𝐣, 𝐤
export Rotor, RotorF64, RotorF32, RotorF16
export Quat3Vec, Quat3VecF64, Quat3VecF32, Quat3VecF16
export abs2vec, absvec
export randn_rotor
export from_float_array, to_float_array, from_euler_angles, to_euler_angles,
    from_euler_phases, to_euler_phases!, to_euler_phases,
    from_spherical_coordinates, to_spherical_coordinates,
    from_rotation_matrix, to_rotation_matrix
export distance, distance_rotation
export unflip, unflip!

abstract type AbstractQuaternion{T<:Real} <: Number end

include("quaternion.jl")
include("algebra.jl")
include("math.jl")
include("random.jl")
include("conversion.jl")
include("distance.jl")
include("interpolation.jl")

end  # module
