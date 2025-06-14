using AutoHist, Distributions
using Test

import StatsBase: Histogram, fit

@testset "return type" begin
    x = collect(LinRange(0,1,11))

    for rule in [:bayes, :aic, :bic, :br, :mdl, :nml, :l2cv, :klcv,
                 :sturges, :fd, :scott]
        for closed in [:left, :right]
            H = histogram_regular(x; rule=rule, closed=closed)
            @test typeof(H) <: Histogram
        end
    end
    for scalest in [:minim, :iqr, :stdev]
        for level in [0,1,2,3,4,5]
            for closed in [:left, :right]
                H = histogram_regular(x; rule=:wand, closed=closed, scalest=scalest, level=level)
                @test typeof(H) <: Histogram
            end
        end
    end
    H = histogram_regular(x; rule=:wand)
    @test typeof(H) <: Histogram

    for rule in [:pena, :penb, :penr, :bayes, :klcv, :l2cv, :nml]
        H = histogram_irregular(x; rule=rule)
        @test typeof(H) <: Histogram
    end
    for grid in [:regular, :data, :quantile] # test grid, right-left open interval combinations
        for closed in [:left, :right]
            H = histogram_irregular(x; grid=grid, closed=closed)
            @test typeof(H) <: Histogram
        end
    end
    H = histogram_irregular(x; greedy=false) # check that greedy != false works
    @test typeof(H) <: Histogram
end

@testset "left open and right open intervals" begin
    x = collect(LinRange(0,1,11))
    H1 = histogram_regular(x; closed=:right)
    H2 = histogram_regular(x; closed=:left)
    H3 = histogram_irregular(x; closed=:right)
    H4 = histogram_irregular(x; closed=:left)

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
    # In this case, the function should throw an appropriate error
    x = randn(10^3)

    @test_throws ArgumentError histogram_regular(x; rule=:nonsense, maxbins=-100)
    @test_throws DomainError histogram_regular(x; rule=:bayes, a=-1.0)
    @test_throws DomainError histogram_regular(x; rule=:bayes, a=k->-2.0*k)
    @test_throws DomainError histogram_irregular(x; a=-1.0)
    @test_throws ArgumentError histogram_irregular(x; rule=:nonsense)
    @test_throws ArgumentError histogram_irregular(x; grid=:nonsense)
    @test_throws ArgumentError histogram_regular(x; rule=:wand, scalest=:nonsense)
    @test_throws ArgumentError histogram_irregular(x; closed=:nonsense)
    @test_throws ArgumentError histogram_regular(x; closed=:nonsense)

end

@testset "a as function" begin
    x = randn(10^3)

    @test histogram_regular(x; a = k->8.0) == histogram_regular(x; a=8.0)
end

@testset "min length" begin
    n = 10^3
    x = randn(n)
    H1 = histogram_irregular(x; rule=:klcv, use_min_length=true)
    H2 = histogram_irregular(x; rule=:l2cv, use_min_length=true)

    min_length = (maximum(x) - minimum(x))*log(n)^(1.5)/n

    @test minimum(H1.edges[1][2:end] - H1.edges[1][1:end-1]) ≥ min_length
    @test minimum(H2.edges[1][2:end] - H2.edges[1][1:end-1]) ≥ min_length
end

@testset "throws error misspecified support" begin
    x = [-1.0, 1.0]
    @test_throws DomainError histogram_irregular(x; support=(-0.5, Inf))
    @test_throws DomainError histogram_irregular(x; support=(-Inf, 0.5))
end