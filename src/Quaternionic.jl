module Quaternionic

import LinearAlgebra: LinearAlgebra, Symmetric, eigen, norm, normalize, (⋅), (×)
import GenericLinearAlgebra
import PrecompileTools: PrecompileTools, @compile_workload, @setup_workload
import StaticArrays: StaticArrays, @SMatrix, @SVector, SA, SMatrix, SVector
import Random: AbstractRNG, default_rng

# The `public` keyword is a syntax error before Julia 1.11, so we wrap it in a macro that
# emits `public` where it exists and does nothing otherwise.
macro public(names)
    @static if VERSION ≥ v"1.11.0-DEV.469"
        syms = names isa Symbol ? (names,) : names.args
        esc(Expr(:public, syms...))
    else
        nothing
    end
end

export AbstractQuaternion
export Quaternion, quaternion,
    QuaternionF64, QuaternionF32, QuaternionF16,
    imx, imy, imz, 𝐢, 𝐣, 𝐤
export Rotor, rotor, RotorF64, RotorF32, RotorF16
export QuatVec, quatvec, QuatVecF64, QuatVecF32, QuatVecF16
export components, basetype
@public value, iszerovalue
@public dominant_eigenvector
@public ℂconj, ℂreal, ℂimag, ℂreim
@public RB, BR, Rv, vR, KAN
export (⋅), (×), (×̂), normalize, norm
export abs2vec, absvec
export from_float_array, to_float_array,
    from_euler_angles, to_euler_angles,
    from_euler_phases, to_euler_phases!, to_euler_phases,
    from_spherical_coordinates, to_spherical_coordinates,
    from_rotation_matrix, to_rotation_matrix
export distance, distance2
export align
export Lorentz, Boost, ga_components
export unflip, unflip!, slerp, squad, squad!
export precessing_nutating_example

abstract type AbstractQuaternion{T<:Number} <: Number end


include("utilities.jl")
include("quaternion.jl")
include("base.jl")
include("algebra.jl")
include("math.jl")
include("exp_log_derivatives.jl")
include("random.jl")
include("conversion.jl")
include("distance.jl")
include("alignment.jl")
include("interpolation.jl")
include("examples.jl")
include("Lorentz.jl")

include("precompilation.jl")

end  # module
