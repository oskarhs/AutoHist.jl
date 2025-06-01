using AutoHist, Distributions
using Test

import StatsBase: Histogram, fit

@testset "return type" begin
    x = collect(LinRange(0,1,11))

    for rule in ["bayes", "aic", "bic", "br", "mdl", "nml", "l2cv", "klcv"]
        H = histogram_regular(x; rule=rule)
        @test typeof(H) <: Histogram
    end
    for rule in ["pena", "penb", "penr", "bayes", "klcv", "l2cv", "nml"]
        H = histogram_irregular(x; rule=rule)
        @test typeof(H) <: Histogram
    end
    for grid in ["regular","data", "quantile"]
        H = histogram_irregular(x; grid=grid)
        @test typeof(H) <: Histogram
    end
end

@testset "left open and right open intervals" begin
    x = collect(LinRange(0,1,11))
    H1 = histogram_regular(x; right=true)
    H2 = histogram_regular(x; right=false)
    H3 = histogram_irregular(x; right=true)
    H4 = histogram_irregular(x; right=false)

    @test H1.closed == :right
    @test H2.closed == :left
    @test H3.closed == :right
    @test H4.closed == :left
end

@testset "estimated support" begin
    n = 100
    x = [-5.0, 4.5]
    H1 = histogram_regular(x)
    edges1 = H1.edges[1]
    H2 = histogram_irregular(x)
    edges2 = H2.edges[1]

    @test isapprox(edges1[1], minimum(x); atol=1e-10)
    @test isapprox(edges1[end], maximum(x); atol=1e-10)
    @test isapprox(edges2[1], minimum(x); atol=1e-10)
    @test isapprox(edges2[end], maximum(x); atol=1e-10)
end

@testset "given support" begin
    n = 100
    x = rand(n)
    H1 = histogram_regular(x; support=(0.0, 1.0))
    edges1 = H1.edges[1]
    H2 = histogram_irregular(x; support=(0.0, 1.0))
    edges2 = H2.edges[1]

    @test isapprox(edges1[1], 0.0; atol=1e-10)
    @test isapprox(edges1[end], 1.0; atol=1e-10)
    @test isapprox(edges2[1], 0.0; atol=1e-10)
    @test isapprox(edges2[end], 1.0; atol=1e-10)
end

@testset "is density" begin
    n = 100
    x = randn(n)

    H1 = histogram_regular(x)
    edges1 = H1.edges[1]
    dens1 = H1.weights
    H2 = histogram_irregular(x)
    edges2 = H2.edges[1]
    dens2 = H2.weights

    @test sum(dens1 .* (edges1[2:end] - edges1[1:end-1])) ≈ 1.0
    @test sum(dens2 .* (edges2[2:end] - edges2[1:end-1])) ≈ 1.0
end

@testset "erroneous to default" begin
    # The arguments are of the right type, but are set to nonsensical values.
    # In this case, the function should use the default parameter values instead.
    x = randn(10^3)
    H1 = histogram_regular(x; rule="nonsense", a=-1.0, maxbins=-100)
    H2 = histogram_irregular(x; rule="nonsense", grid="nonsense", a=-1.0)

    @test H1 == histogram_regular(x)
    @test H2 == histogram_irregular(x)
end

@testset "a as function" begin
    x = randn(10^3)
    H1 = histogram_regular(x; a = k->8.0)
    H2 = histogram_regular(x; a = k->-1.0) # should return a = 1.0 for all k where a(k) is negative

    @test H1 == histogram_regular(x; a=8.0)
    @test H2 == histogram_regular(x; a=1.0)
end