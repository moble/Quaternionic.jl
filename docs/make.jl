# Run with
#   time julia --project=. make.jl && julia --project=. -e 'using LiveServer; serve(dir="build")'
# assuming you are in this `docs` directory (otherwise point the project argument here)

using Quaternionic
using Documenter
using DocumenterCitations
using LinearAlgebra

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "references.bib");
    #style=:authoryear,
)

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)

include("local_notes.jl")
(notes_pages, notes_remotes) = local_notes()

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
    remotes=notes_remotes,
    pages=[
        "Introduction" => "index.md",
        "Manual" => [
            "Types and construction" => "types.md",
            "Algebra and mathematical functions" => "math.md",
            "Rotations and conversions" => "conversions.md",
            "Distances and alignment" => "distances.md",
            "Lorentz transformations" => "lorentz.md",
        ],
        "Functions of time" => "functions_of_time.md",
        "Differentiating by quaternionic arguments" => "differentiation.md",
        "Theory" => [
            "Geometric algebra" => "geometric_algebra.md",
            "Spacetime algebra" => "spacetime_algebra.md",
        ],
        "References" => "references.md",
        notes_pages...,
        "API index" => "api.md",
    ],
    # Every exported or public name must be documented.  Internal helpers, such as
    # `squad_control_points`, may have docstrings that are deliberately left out.
    checkdocs=:public,
    # doctest = false,
    doctestfilters = [
        # Drop any digit after the 12th digit after a decimal, throughout the docs
        r"(?<=\d\.\d{12})\d+",
        # Ignore any warning involving Symbolics or SymbolicUtils
        r"WARNING: .* Symbolic.*",
    ],
)

deploydocs(;
    repo="github.com/moble/Quaternionic.jl",
    devbranch="main",
    push_preview=true,
)
