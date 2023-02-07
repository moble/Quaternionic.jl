using Pkg
try
    Pkg.test("Quaternionic"; coverage=true, test_args=ARGS)
catch e
    println("Tests failed; proceeding to coverage")
end

Pkg.activate()
using Coverage
cd(@__DIR__)
coverage = Coverage.process_folder()
Coverage.writefile("lcov.info", coverage)
Coverage.clean_folder(".")
