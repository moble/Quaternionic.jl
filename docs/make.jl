using Documenter, Quaternionic

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)

makedocs(
    sitename="Quaternionic.jl",
    modules = [Quaternionic]
)
