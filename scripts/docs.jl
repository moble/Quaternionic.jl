# Run this script from any directory as
#
#   julia -t 4 scripts/docs.jl
#
# The docs will build and the browser should open automatically.  `LiveServer`
# will monitor the docs for any changes, then rebuild them and refresh the browser
# until this script is stopped.

import Dates
println("Building docs starting at ", Dates.format(Dates.now(), "HH:MM:SS"), ".")

using Pkg
cd((@__DIR__) * "/..")
Pkg.activate("docs")
# The docs share the workspace's root Manifest.toml with the tests, so this only
# instantiates it; updating it would also change the packages that the tests use.
Pkg.instantiate()

using LiveServer
servedocs(launch_browser=true)
