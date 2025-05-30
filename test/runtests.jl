using AutoHist, Distributions
using Test

import StatsBase: Histogram, fit

@testset "bin regular" begin
    x = rand(10^3)
    k = 10
    H = fit(Histogram, x, LinRange(0.0, 1.0, k+1), closed=:right)
    N = bin_regular(x, 0.0, 1.0, k, true)
    @test all(N .≈ H.weights)
end

@testset "bin irregular" begin
    x = rand(10^3)
    edges = [0.0, 0.2, 0.5, 0.85, 1.0]
    H = fit(Histogram, x, edges, closed=:right)
    N = bin_irregular(x, edges, true)
    @test all(N .≈ H.weights)
end

@testset "return type" begin
    x = collect(LinRange(0,1,11))
    H1 = histogram_regular(x)
    H2 = histogram_irregular(x)

    @test typeof(H1) <: Histogram
    @test typeof(H2) <: Histogram
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