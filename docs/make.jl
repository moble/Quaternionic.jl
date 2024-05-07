# Run with
#   time julia --project=. make.jl && julia --project=. -e 'using LiveServer; serve(dir="build")'
# assuming you are in this `docs` directory (otherwise point the project argument here)

using Quaternionic
using Documenter
using DocumenterCitations

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    #style=:authoryear,
)

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)

makedocs(;
    plugins=[bib],
    sitename="Quaternionic.jl",
    modules=[Quaternionic],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        canonical="https://moble.github.io/Quaternionic.jl/stable/",
        assets = String["assets/citations.css"],
    ),
    authors="Michael Boyle <michael.oliver.boyle@gmail.com>",
    repo=Remotes.GitHub("moble", "Quaternionic.jl"),
    pages=[
        "Introduction" => "index.md",
        "Basics" => "manual.md",
        "Functions of time" => "functions_of_time.md",
        "Differentiating by quaternions" => "differentiation.md",
        "All functions" => "functions.md"
    ],
    # doctest = false
)

deploydocs(;
    repo="github.com/moble/Quaternionic.jl",
    devbranch="main",
    push_preview=true,
)
