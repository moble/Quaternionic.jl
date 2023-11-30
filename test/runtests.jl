"""
Run this script (from this directory) as

    time julia --code-coverage=tracefile-%p.info --code-coverage=user --project=. ./runtests.jl

Then, if you have `lcov` installed, you should also have `genhtml`, and you can run this

    genhtml tracefile-<your_PID>.info --output-directory coverage/ && open coverage/index.html

to view the coverage locally as HTML.  I find that this sometimes requires removing files
that aren't really there from the .info file.

It's a well-hidden fact that you can turn coverage on and off by adding certain comments
around the code you don't want to measure:

    # COV_EXCL_START
    untested_code_that_wont_show_up_in_coverage()
    # COV_EXCL_STOP

"""

# This is hopefully a temporary clause, for checking the GC behavior on Julia pre-1.10 nightlies,
# which are failing for unknown reasons, but in ways that resemble behavior other people have found
# and traced back to the new GC.  This block's code is copied from AbstractAlgebra.jl (under the
# Simplified "2-clause" BSD License): https://github.com/Nemocas/AbstractAlgebra.jl/pull/1432/files
if VERSION >= v"1.8.0"
  GC.enable_logging(true)

  # print gc settings
  jlmax = @ccall jl_gc_get_max_memory()::UInt64
  totalmem = @ccall uv_get_total_memory()::UInt64
  constrmem = @ccall uv_get_constrained_memory()::UInt64
  println("Memory:   max: ", Base.format_bytes(jlmax))
  println("        total: ", Base.format_bytes(totalmem))
  println("       constr: ", Base.format_bytes(constrmem))

#= FIXME/TODO: in the future we may wish to experiment with limiting the GC heap here
  if VERSION >= v"1.10.0-"
    # adjust heap size hint
    memenv = parse(Int, get(ENV, "OSCARCI_MAX_MEM_GB", "5")) * 2^30
    println("Setting heap size hint to ", Base.format_bytes(memenv))
    @ccall jl_gc_set_max_memory(memenv::UInt64)::Cvoid
  end
=#
end



using Quaternionic
using Test
using Random, Symbolics, StaticArrays, ForwardDiff, GenericLinearAlgebra,
    ChainRulesTestUtils, Zygote, ChainRulesTestUtils

import LinearAlgebra

using ChainRulesCore
ChainRulesCore.debug_mode() = true

@variables w x y z a b c d e  # Symbolic variables

# NOTE: `FloatTypes` and `IntTypes` must be in descending order of width
FloatTypes = [BigFloat, Float64, Float32, Float16]
IntTypes = [BigInt, Int128, Int64, Int32, Int16, Int8]
SymbolicTypes = [Num]
Types = [FloatTypes...; IntTypes...; SymbolicTypes...]
PrimitiveTypes = [T for T in Types if isbitstype(T)]

QTypes = [Quaternion, Rotor, QuatVec]

# Handy assignments for now
Base.eps(::Quaternion{T}) where {T} = eps(T)
Base.eps(T::Type{<:Integer}) = zero(T)
Base.eps(n::Num) = zero(n)
Base.:≈(a::Num, b::Num; kwargs...) = iszero(Symbolics.simplify(a-b; expand=true))

enabled_tests = lowercase.(ARGS)

help = ("help" ∈ enabled_tests || "--help" ∈ enabled_tests)
helptests = []

# This block is cribbed from StaticArrays.jl/test/runtests.jl
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if help
        push!(helptests, key)
    else
        if isempty(enabled_tests) || key in enabled_tests
            println("Running $key.jl")
            Random.seed!(42)
            include(fname)
        end
    end
end

@testset verbose=true "All tests" begin
    addtests("quaternion.jl")
    addtests("basis.jl")
    addtests("base.jl")
    addtests("algebra.jl")
    addtests("math.jl")
    addtests("random.jl")
    addtests("conversion.jl")
    addtests("distance.jl")
    addtests("alignment.jl")
    addtests("interpolation.jl")
    addtests("gradients.jl")
    addtests("auto_differentiation.jl")
    addtests("doctests.jl")
end

if help
    println()
    println("Pass no args to run all tests, or select one or more of the following:")
    for helptest in helptests
        println("    ", helptest)
    end
end
