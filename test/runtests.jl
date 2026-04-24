#!/usr/bin/env julia

"""
Run this file from the root of the package as

    julia --project=test test/runtests.jl

To filter to specific test files (substring match on the filename):

    julia --project=test test/runtests.jl --file algebra
    julia --project=test test/runtests.jl --file lorentz --verbose

Filters for @testitem-based tests (e.g., differentiation_interface.jl):

    julia --project=test test/runtests.jl --tags unit,fast
    julia --project=test test/runtests.jl --name "some test name"
    julia --project=test test/runtests.jl --pattern "some pattern"
    julia --project=test test/runtests.jl --exclude slow

For code coverage, run from the root as

    julia --project=test --code-coverage=tracefile-%p.info --code-coverage=user test/runtests.jl

Then, if you have `lcov` installed, you should also have `genhtml`, and you can run

    genhtml tracefile-<your_PID>.info --output-directory coverage/ && open coverage/index.html

to view the coverage locally as HTML.  I find that this sometimes requires removing files
that aren't really there from the .info file.

It's a well-hidden fact that you can turn coverage on and off by adding certain comments
around the code you don't want to measure:

    # COV_EXCL_START
    untested_code_that_wont_show_up_in_coverage()
    # COV_EXCL_STOP

"""

using Quaternionic
using Test
using TestItemRunner
using Random, StaticArrays, ForwardDiff, GenericLinearAlgebra,
    ChainRulesTestUtils, Zygote, ChainRulesTestUtils, Aqua
import Symbolics, FastDifferentiation, Latexify
import LinearAlgebra
using ChainRulesCore


Symbolics.@variables w x y z a b c d e  # Symbolic variables

# NOTE: `FloatTypes` and `IntTypes` must be in descending order of width
FloatTypes = [BigFloat, Float64, Float32, Float16]
IntTypes = [BigInt, Int128, Int64, Int32, Int16, Int8]
SymbolicTypes = [Symbolics.Num]
Types = [FloatTypes...; IntTypes...; SymbolicTypes...]
PrimitiveTypes = [T for T in Types if isbitstype(T)]

QTypes = [Quaternion, Rotor, QuatVec]

# Handy assignments for now
Base.eps(::Quaternion{T}) where {T} = eps(T)
Base.eps(T::Type{<:Integer}) = zero(T)
Base.eps(n::Symbolics.Num) = zero(n)
Base.:≈(a::Symbolics.Num, b::Symbolics.Num; kwargs...) =
    iszero(Symbolics.simplify(a-b; expand=true))
function _sym_iszero(diff)
    d = Symbolics.simplify(diff; expand=true)
    iszero(d) || iszero(Symbolics.simplify(d^2; expand=true))
end
Base.:≈(a::Symbolics.Num, b::Number; kwargs...) = _sym_iszero(a - b)
Base.:≈(a::Number, b::Symbolics.Num; kwargs...) = _sym_iszero(a - b)
Base.:≈(a::AbstractQuaternion{Symbolics.Num}, b::AbstractQuaternion{Symbolics.Num}; kwargs...) =
    all(iszero(Symbolics.simplify(x - y; expand=true)) for (x, y) in zip(components(a), components(b)))

# ─── CLI ──────────────────────────────────────────────────────────────────────

function _print_help()
    @info """Command-line runner for Quaternionic.jl tests

    Basic usage (from the root of the package):
        julia --project=test test/runtests.jl                               # Run all tests
        julia --project=test test/runtests.jl --help                        # Show this help
        julia --project=test test/runtests.jl --verbose                     # Verbose output
        julia --project=test test/runtests.jl --list-tags                   # List tags

    File filter (substring match; applies to include()-based test modules):
        julia --project=test test/runtests.jl --file algebra                # Run algebra.jl
        julia --project=test test/runtests.jl --file differentiation_interface

    Filters for @testitem-based tests (e.g., differentiation_interface.jl):
        julia --project=test test/runtests.jl --tags unit,fast              # All tags must match
        julia --project=test test/runtests.jl --name "Some test name"       # Name substring
        julia --project=test test/runtests.jl --pattern "Some pattern"      # Name or filename
        julia --project=test test/runtests.jl --exclude slow                # Exclude any tag

    Multiple filters can be combined (all must be satisfied):
        julia --project=test test/runtests.jl --file algebra --verbose
    """
    return
end

function _list_available_tags()
    @info "Available tags for @testitem tests:"
    println()
    for (tag, desc) in TAGS_DATA
        if isempty(desc)
            println("  $tag")
        else
            println("  $tag - $desc")
        end
    end
    return
end

const TAGS_DATA = Dict(
    # Test Type (What kind of test?)
    :integration => "End-to-end tests with real datasets and full workflows",
    :unit => "Single component or function tests",
    :validation => "Tests verifying expected values, behavior, or mathematical correctness",

    # Complexity (How resource-intensive?)
    :fast => "Quick tests suitable for frequent execution",
    :slow => "Resource-intensive tests requiring significant time or memory",
)

"""
    _parse_argument_with_value(flag, transform = identity)

Parse a command-line argument that expects a value following the flag.
If multiple occurrences exist, uses the last one (warns about duplicates).
"""
function _parse_argument_with_value(flag, transform = identity)
    occurrences = findall(x -> x == flag, ARGS)
    isempty(occurrences) && return nothing

    if length(occurrences) > 1
        @warn "Duplicate argument '$flag' found, using last occurrence"
    end

    idx = last(occurrences)
    if idx == length(ARGS)
        error("Missing argument for '$flag'")
    end
    try
        return transform(ARGS[idx+1])
    catch
        error("Invalid value for flag '$flag': $(ARGS[idx + 1])")
    end
end

"""
    _validate_arguments()

Validate that all command-line arguments are recognized.
"""
function _validate_arguments()
    valid_flags = Set([
        "--verbose", "-v",
        "--help", "-h",
        "--list-tags", "-l",
        "--file",
        "--tags",
        "--exclude",
        "--name",
        "--pattern",
    ])

    i = 1
    while i <= length(ARGS)
        arg = ARGS[i]
        if startswith(arg, "-")
            if !(arg in valid_flags)
                error("Unknown argument: $arg")
            end
            if arg in ["--file", "--tags", "--exclude", "--name", "--pattern"]
                i += 1  # Skip the value
            end
        end
        i += 1
    end
end

"""
    parse_arguments()

Parse command-line arguments for the test runner.  Returns a named tuple with:
- `verbose`: Enable verbose @testitem output
- `help`: Show help and exit
- `list`: List available tags and exit
- `file`: Run only test files whose name contains this substring
- `tags`: Run @testitems that have ALL of these tags (AND logic)
- `exclude`: Skip @testitems that have ANY of these tags (OR logic)
- `name`: Run @testitems whose name contains this substring
- `pattern`: Run @testitems whose name OR filename contains this substring
"""
function parse_arguments()
    _validate_arguments()

    verbose = "--verbose" in ARGS || "-v" in ARGS
    help = "--help" in ARGS || "-h" in ARGS
    list = "--list-tags" in ARGS || "-l" in ARGS

    file_filter = _parse_argument_with_value("--file")
    name_filter = _parse_argument_with_value("--name")
    pattern_filter = _parse_argument_with_value("--pattern")

    function ensure_tag_existence(tag)
        if !haskey(TAGS_DATA, tag)
            error(
                "Tag '$tag' is not a valid tag. Update `TAGS_DATA` in `test/runtests.jl` if necessary",
            )
        end
        return tag
    end

    tag_transform(list_of_tags) =
        map(split(list_of_tags, ",")) do tag
            ensure_tag_existence(Symbol(tag))
        end

    tags_filter = _parse_argument_with_value("--tags", tag_transform)
    exclude_filter = _parse_argument_with_value("--exclude", tag_transform)

    return (
        verbose = verbose,
        help = help,
        list = list,
        file = file_filter,
        tags = tags_filter,
        exclude = exclude_filter,
        name = name_filter,
        pattern = pattern_filter,
    )
end

function _create_filter(args)
    filters = []

    if !isnothing(args.file)
        push!(filters, test_item -> contains(test_item.filename, args.file))
    end
    if !isnothing(args.tags)
        push!(filters, test_item -> all(tag in test_item.tags for tag in args.tags))
    end
    if !isnothing(args.exclude)
        push!(filters, test_item -> !(any(tag in test_item.tags for tag in args.exclude)))
    end
    if !isnothing(args.name)
        push!(filters, test_item -> contains(test_item.name, args.name))
    end
    if !isnothing(args.pattern)
        push!(
            filters,
            test_item ->
                contains(test_item.name, args.pattern) ||
                contains(test_item.filename, args.pattern),
        )
    end

    isempty(filters) && return nothing
    return test_item -> all(f(test_item) for f in filters)
end


# ─── Test execution ───────────────────────────────────────────────────────────

args = parse_arguments()

if args.help
    _print_help()
elseif args.list
    _list_available_tags()
else
    filter_func = _create_filter(args)

    # Run a file-based test module if its name matches the --file filter
    function addtests(fname)
        key = lowercase(splitext(fname)[1])
        if isnothing(args.file) || contains(key, args.file)
            println("Running $fname")
            Random.seed!(42)
            include(joinpath(@__DIR__, fname))
        end
    end

    @testset verbose=true "All tests" begin

        # @testitem-based tests (differentiation_interface.jl and any @testitem blocks
        # embedded in the package source).  These support the full filter set.
        if isnothing(args.file) || contains("differentiation_interface", args.file)
            println("Running differentiation_interface.jl")
            if isnothing(filter_func)
                @run_package_tests verbose = args.verbose
            else
                @run_package_tests verbose = args.verbose filter = filter_func
            end
        end

        addtests("aqua.jl")
        addtests("quaternion.jl")
        addtests("basis.jl")
        addtests("base.jl")
        addtests("algebra.jl")
        addtests("math.jl")
        addtests("lorentz.jl")
        addtests("lorentz_group.jl")
        addtests("random.jl")
        addtests("conversion.jl")
        addtests("distance.jl")
        addtests("alignment.jl")
        addtests("interpolation.jl")
        addtests("gradients.jl")
        addtests("auto_differentiation.jl")
        addtests("doctests.jl")
    end
end
