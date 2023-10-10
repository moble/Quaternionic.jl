import Dates
println("Running tests starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate(".")

try
    Δt = @elapsed Pkg.test("Quaternionic"; coverage=true, test_args=ARGS)
    println("Running tests took $Δt seconds.")
catch e
    println("Tests failed; proceeding to coverage")
end

Pkg.activate()  # Activate Julia's base (home) directory
using Coverage
cd((@__DIR__) * "/..")
coverage = vcat(Coverage.process_folder("src"), Coverage.process_folder("ext"))
Coverage.writefile("lcov.info", coverage)
Coverage.clean_folder(".")
