using AutoHist, StatsBase, Random, Distributions

import BenchmarkTools: @benchmark

function benchmark_autohist()
    rng = Xoshiro(1812)
    n = 25000
    x = rand(rng, Normal(), n)

    b_reg = @benchmark histogram_regular(x)
    @show b_reg

    b_irr = @benchmark histogram_irregular(x; grid="data")
    @show b_irr
    return nothing
end

benchmark_autohist()