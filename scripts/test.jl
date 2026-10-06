import Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate(".")

# Coverage is written even when the tests fail, but the failure is reported and the script
# then exits with a nonzero status.
testerror = nothing
try
    Δt = @elapsed Pkg.test("Quaternionic"; coverage=true, test_args=ARGS)
    println("Running tests took $Δt seconds.")
catch e
    e isa InterruptException && rethrow()
    global testerror = e
    showerror(stderr, e)
    println(stderr, "\nTests failed; proceeding to coverage")
end

Pkg.activate("test")  # Coverage is a dependency of the test environment.
using Coverage
cd((@__DIR__) * "/..")
coverage = vcat(Coverage.process_folder("src"), Coverage.process_folder("ext"))
Coverage.writefile("lcov.info", coverage)
# Pkg.test writes .cov files only next to the package's own sources.  Cleaning just these
# directories avoids following symbolic links out of the repository.
foreach(Coverage.clean_folder, ("src", "ext", "test"))

isnothing(testerror) || exit(1)
