# This simply uses `make.jl` in this directory to build the docs, then serves them locally.
# Run `julia --project=docs docs/serve.jl` from the top directory to execute this script.

#using Quaternionic
import LiveServer: servedocs

servedocs(; launch_browser=true)
