@testitem "Doctests" tags=[:unit] begin
    using Documenter
    # Documenter evaluates the manual's `@meta` blocks, such as `CurrentModule =
    # Quaternionic`, in `Main`, so the package must be loaded there.
    @eval Main using Quaternionic
    DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)
    doctest(Quaternionic)
end
