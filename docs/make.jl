using Documenter, AutoHist

makedocs(
    sitename="AutoHist.jl",
    modules = [AutoHist],
    format=Documenter.HTML(),
    pages = [
        "Introdution" => "index.md",
        "Supported Methods" => "methods.md",
        "API" => "api.md"
    ]
)

deploydocs(;
    repo="git@github.com/oskarhs/AutoHist.jl",
    push_preview = true
)