using Test, Documenter, Quaternionic

DocMeta.setdocmeta!(Quaternionic, :DocTestSetup, :(using Quaternionic); recursive=true)
doctest(Quaternionic)
