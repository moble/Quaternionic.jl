using Quaternionic
using Documenter

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)

makedocs(;
    modules=[Quaternionic],
    authors="Michael Boyle <michael.oliver.boyle@gmail.com>",
    repo="https://github.com/moble/Quaternionic.jl/blob/{commit}{path}#{line}",
    sitename="Quaternionic.jl",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://moble.github.io/Quaternionic.jl/stable/",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Basics" => "manual.md",
        "Functions of time" => "functions_of_time.md",
        "Differentiating by quaternions" => "differentiation.md",
    ],
    # doctest = false
)

deploydocs(;
    repo="github.com/moble/Quaternionic.jl",
    devbranch="main",
    push_preview=true,
)
