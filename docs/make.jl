using Documenter, AutoHist, StatsBase

makedocs(
    sitename="AutoHist.jl",
    modules = [AutoHist],
    format=Documenter.HTML(),
    pages = [
        "Introdution" => "index.md",
        "Supported Methods" => "methods.md",
        "Examples" => [
            "examples/density_estimation.md",
            "examples/algorithm_choice.md"
        ],
        "API" => "api.md",
        "Algorithms" => "algorithms.md"
    ],
    checkdocs=:none
)

deploydocs(;
    repo="git@github.com/oskarhs/AutoHist.jl",
    push_preview = true
)