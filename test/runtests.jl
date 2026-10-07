# Run every test item under `Pkg.test`.
#
# Day-to-day test runs are better made with the test-item tools (the VS Code test explorer,
# `juliati`, or the JuliaMCP server), which can select items by name, file, or tag, and
# which keep their worker processes alive between runs.  This file exists for `Pkg.test`,
# which CI uses.  It accepts two filters on the tags of the test items, which are passed as
# `test_args` to `Pkg.test`, or on the command line when this file is run directly:
#
#     julia --project=test test/runtests.jl --tags slow          # items with all these tags
#     julia --project=test test/runtests.jl --exclude slow,ad    # items with none of them
#
# The tags in use are
#
#   :unit         single component or function tests
#   :validation   tests of expected values, behavior, or mathematical correctness
#   :fast         quick tests suitable for frequent execution
#   :slow         resource-intensive tests requiring significant time or memory
#   :ad           automatic-differentiation tests, along with one tag for each set of rules
#                 or backend: :chainrules, :zygote, :forwarddiff, :enzyme, :mooncake, and
#                 :reversediff
#
# For code coverage, run from the root of the package as
#
#     julia --project=test --code-coverage=tracefile-%p.info --code-coverage=user test/runtests.jl
#
# Then, if you have `lcov` installed, you should also have `genhtml`, and you can run
#
#     genhtml tracefile-<your_PID>.info --output-directory coverage/ && open coverage/index.html
#
# to view the coverage locally as HTML.  This sometimes requires removing files that aren't
# really there from the .info file.  Code can be excluded from the coverage measurement by
# surrounding it with these comments:
#
#     # COV_EXCL_START
#     untested_code_that_wont_show_up_in_coverage()
#     # COV_EXCL_STOP

using TestItemRunner

function tag_filter(args)
    flags = ("--tags", "--exclude")
    tags = Dict(flag => Symbol[] for flag in flags)
    i = 1
    while i ≤ length(args)
        flag = args[i]
        (flag ∈ flags && i < length(args)) ||
            error("Unrecognized test arguments $(args[i:end]); only `--tags` and `--exclude` are accepted")
        append!(tags[flag], Symbol.(split(args[i+1], ',')))
        i += 2
    end
    testitem -> all(∈(testitem.tags), tags["--tags"]) && !any(∈(testitem.tags), tags["--exclude"])
end

@run_package_tests filter=tag_filter(ARGS)
