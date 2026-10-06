# Derivative rules (`rrule`s and `frule`s) and projections for ChainRules consumers, such as
# Zygote.  This file is loaded as a package extension on Julia 1.9 and later, and by Requires
# (included from `src/Quaternionic.jl`, so nested inside `Quaternionic`) on earlier versions.
#
# The conventions for tangents and cotangents are documented at the top of `projection.jl`,
# and the helper functions that every rule file uses are listed at the top of `helpers.jl`.
# Rules and opt-outs live in separate files, because an `@opt_out` with exactly the signature
# of a supplied rule silently deletes that rule.
module QuaternionicChainRulesCoreExt

# Under Requires, this file is included into `Quaternionic` itself, so both `Quaternionic`
# and `ChainRulesCore` (which Requires binds there) are reached as relative imports.
isdefined(Base, :get_extension) ?
    (using Quaternionic; import Quaternionic: basetype, value, iszerovalue) :
    (using ..Quaternionic; import ..Quaternionic: basetype, value, iszerovalue)
isdefined(Base, :get_extension) ?
    (using ChainRulesCore; import ChainRulesCore: rrule, frule, ProjectTo, @opt_out) :
    (using ..ChainRulesCore; import ..ChainRulesCore: rrule, frule, ProjectTo, @opt_out)
using StaticArrays: SVector
using LinearAlgebra: LinearAlgebra, dot, norm, normalize

include("helpers.jl")
include("projection.jl")
include("structural.jl")
include("arithmetic.jl")
include("elementary.jl")
include("geometry.jl")
include("optouts.jl")

end  # module
