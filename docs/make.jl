using BetaEnsemblesSampling
using Documenter

DocMeta.setdocmeta!(BetaEnsemblesSampling, :DocTestSetup, :(using BetaEnsemblesSampling); recursive=true)

makedocs(;
    modules=[BetaEnsemblesSampling],
    authors="Dan Betea",
    repo="https://github.com/danbetea/BetaEnsemblesSampling.jl/blob/{commit}{path}#{line}",
    sitename="BetaEnsemblesSampling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://danbetea.github.io/BetaEnsemblesSampling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/danbetea/BetaEnsemblesSampling.jl",
    devbranch="main",
)
