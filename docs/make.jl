using Documenter, Quaternionic

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)

makedocs(
    sitename="Quaternionic.jl",
    modules = [Quaternionic],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),  # Use clean URLs, unless built as a "local" build
        edit_link = "main",  # Link out to "main" branch on github
        #= UNCOMMENT ONCE STABLE:
        canonical = "https://moble.github.io/Quaternionic.jl/stable/",
        warn_outdated = true,
        :UNCOMMENT ONCE STABLE =#
    ),
    pages = [
        "Introduction" => "index.md",
        "Manual" => "manual.md",
    ],
)

deploydocs(
    repo="github.com/moble/Quaternionic.jl.git",
    devbranch="main"
)
