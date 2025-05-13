include(joinpath(@__DIR__, "..", "src", "AutoHist.jl"))
using .AutoHist
using StatsBase, Random, Distributions, Plots

import BenchmarkTools: @benchmark, @btime

function benchmark_autohist()
    rng = Xoshiro(1812)
    n = 10^6
    x = rand(rng, Normal(), n)

    @btime histogram_regular($x, rule="aic")
    H = histogram_regular(x, rule="aic")

    @btime histogram_irregular($x; grid="data")
    H2 = histogram_irregular(x; grid="data")

    p = plot(H)
    savefig(p, "example.pdf")

    p2 = plot(H2)
    savefig(p2, "example2.pdf")
    return nothing
end

benchmark_autohist()