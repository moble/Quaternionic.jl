# Derivative rules (`rrule`s and `frule`s) and projections for ChainRules consumers, such as
# Zygote.
#
# The conventions for tangents and cotangents are documented at the top of `projection.jl`,
# and the helper functions that every rule file uses are listed at the top of `helpers.jl`.
# Rules and opt-outs live in separate files, because an `@opt_out` with exactly the signature
# of a supplied rule silently deletes that rule.
module QuaternionicChainRulesCoreExt

using Quaternionic
import Quaternionic: basetype, value, iszerovalue
using ChainRulesCore
import ChainRulesCore: rrule, frule, ProjectTo, @opt_out
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
