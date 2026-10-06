# This simply uses `make.jl` in this directory to build the docs, then serves them locally.
# Run `julia --project=docs docs/serve.jl` from the top directory to execute this script.
#
# LiveServer rebuilds the docs whenever a file in `docs/src` or in the package's `src`
# changes.  Docstring edits in `src` appear in the rebuilt docs only if Revise is available
# (for example, from the default global environment), because otherwise the package is not
# reloaded between builds.

# Activate the docs environment even if this script was started with another project (for
# example, `--project=.`), because the doctests need packages such as Symbolics that only
# the docs environment provides.
import Pkg
Pkg.activate(@__DIR__)

try
    using Revise
catch
    @warn "Revise is not available, so docstring edits will not appear until the server is restarted."
else
    using Quaternionic
end

import LiveServer: servedocs

servedocs(;
    include_dirs=["src/", realpath(joinpath("docs", "src", "local_notes"))],
    launch_browser=true,
)
