### These benchmarks can be run from the top-level directory of this repo with
###
###     julia -e 'using PkgBenchmark; results=benchmarkpkg("Quaternionic"); export_markdown("benchmark/results.md", results)'
###
### This runs the benchmarks (possibly tuning them automatically first), and writes the results to a nice markdown file.


using BenchmarkTools
using Random
using Quaternionic

Random.seed!(1234)

const SUITE = BenchmarkGroup()

const N = 100_000

SUITE["algebra"] = BenchmarkGroup()
for T in [Float64, Float32, Float16]
    a = collect(randn(Quaternion{T}, N))
    b = collect(randn(Quaternion{T}, N))
    for op in [:+, :-, :*, :/]
        SUITE["algebra"][op, T] = @benchmarkable $op.($a, $b)
    end
end

SUITE["math"] = BenchmarkGroup()
for T in [Float64, Float32, Float16]
    a = collect(randn(Quaternion{T}, N))
    for f in [abs, abs2, absvec, abs2vec, inv, log, exp, sqrt, angle]
        SUITE["algebra"][f, T] = @benchmarkable $f.($a)
    end
end
