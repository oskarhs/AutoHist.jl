include(joinpath(@__DIR__, "..", "src", "AutoHist.jl"))

using AutoHist, StatsBase, Random, Distributions

import BenchmarkTools: @benchmark

function benchmark_autohist()
    rng = Xoshiro(1812)
    n = 25000
    x = rand(rng, Normal(), n)

    @btime histogram_regular($x)

    @btime histogram_irregular($x; grid="data")
    return nothing
end

benchmark_autohist()