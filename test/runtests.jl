"""
Run this script (from this directory) as

    time julia --code-coverage=tracefile-%p.info --code-coverage=user --project=. ./runtests.jl

Then, if you have `lcov` installed, you should also have `genhtml`, and you can run this

    genhtml tracefile-<your_PID>.info --output-directory coverage/ && open coverage/index.html

to view the coverage locally as HTML.  I find that this sometimes requires
removing files that aren't really there from the .info file.

"""

using Quaternionic
using Test, Random, Symbolics

@variables w x y z a b c d  # Symbolic variables

FloatTypes = [Float64, Float32, Float16, BigFloat]
IntTypes = [Int128, Int64, Int32, Int16, Int8, BigInt]
SymbolicTypes = [Num]
Types = [FloatTypes...; IntTypes...; SymbolicTypes...]

# Handy assignments for now
Base.eps(::Quaternion{T}) where {T} = eps(T)
Base.eps(T::Type{<:Integer}) = zero(T)
Base.eps(n::Num) = zero(n)
Base.:≈(a::Num, b::Num; kwargs...) = Symbolics.simplify(a-b; expand=true) == 0

# This block is cribbed from StaticArrays.jl/test/runtests.jl
#
# Hook into Pkg.test so that tests from a single file can be run.  For example,
# to run only the MVector and SVector tests, use:
#
#   Pkg.test("StaticArrays", test_args=["MVector", "SVector"])
#
enabled_tests = lowercase.(ARGS)
function addtests(fname)
    key = lowercase(splitext(fname)[1])
    if isempty(enabled_tests) || key in enabled_tests
        Random.seed!(42)
        include(fname)
    end
end


addtests("basis.jl")
addtests("fundamentals.jl")
addtests("math.jl")
addtests("doctests.jl")
