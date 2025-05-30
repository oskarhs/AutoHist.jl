using Documenter, AutoHist

makedocs(
    sitename="AutoHist.jl",
    modules = [AutoHist],
    format=Documenter.HTML(),
    pages = [
        "Introdution" => "index.md",
        "API" => "api.md"
    ]
)