using Documenter, AutoHist, StatsBase

makedocs(
    sitename="AutoHist.jl",
    modules = [AutoHist],
    format=Documenter.HTML(; assets=["assets/favicon.ico"]),
    pages = [
        "Introduction" => "index.md",
        "Supported Methods" => "methods.md",
        "Examples" => [
            "examples/plotting.md",
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