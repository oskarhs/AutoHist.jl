using Plots; gr()
using AutoHist, Distributions
using Test

import StatsBase: Histogram, fit

@testset "return type" begin
    x = collect(LinRange(0,1,11))

    for rule in [:bayes, :aic, :bic, :br, :mdl, :nml, :l2cv, :klcv,
                 :sturges, :fd, :scott]
        for closed in [:left, :right]
            H = histogram_regular(x; rule=rule, closed=closed)
            @test typeof(H) <: AutomaticHistogram
        end
    end
    for scalest in [:minim, :iqr, :stdev]
        for level in [0,1,2,3,4,5]
            for closed in [:left, :right]
                H = histogram_regular(x; rule=:wand, closed=closed, scalest=scalest, level=level)
                @test typeof(H) <: AutomaticHistogram
            end
        end
    end
    H = histogram_regular(x; rule=:wand)
    @test typeof(H) <: AutomaticHistogram

    for rule in [:pena, :penb, :penr, :bayes, :klcv, :l2cv, :nml]
        H = histogram_irregular(x; rule=rule)
        @test typeof(H) <: AutomaticHistogram
    end
    for grid in [:regular, :data, :quantile] # test grid, right-left open interval combinations
        for closed in [:left, :right]
            H = histogram_irregular(x; grid=grid, closed=closed)
            @test typeof(H) <: AutomaticHistogram
        end
    end
    H = histogram_irregular(x; greedy=false) # check that greedy != false works
    @test typeof(H) <: AutomaticHistogram
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
    #= edges1 = H1.edges[1]
    H2 = histogram_irregular(x)
    edges2 = H2.edges[1]
    @test isapprox(edges1[1], minimum(x); atol=1e-10)
    @test isapprox(edges1[end], maximum(x); atol=1e-10)
    @test isapprox(edges2[1], minimum(x); atol=1e-10)
    @test isapprox(edges2[end], maximum(x); atol=1e-10) =#

    breaks1 = H1.breaks
    H2 = histogram_irregular(x)
    breaks2 = H2.breaks
    @test isapprox(breaks1[1], minimum(x); atol=1e-10)
    @test isapprox(breaks1[end], maximum(x); atol=1e-10)
    @test isapprox(breaks2[1], minimum(x); atol=1e-10)
    @test isapprox(breaks2[end], maximum(x); atol=1e-10)
end

@testset "given support" begin
    n = 100
    x = rand(n)
    #= H1 = histogram_regular(x; support=(0.0, 1.0))
    edges1 = H1.edges[1]
    H2 = histogram_irregular(x; support=(0.0, 1.0))
    edges2 = H2.edges[1]
    @test isapprox(edges1[1], 0.0; atol=1e-10)
    @test isapprox(edges1[end], 1.0; atol=1e-10)
    @test isapprox(edges2[1], 0.0; atol=1e-10)
    @test isapprox(edges2[end], 1.0; atol=1e-10) =#

    h1 = histogram_regular(x; support=(0.0, 1.0))
    breaks1 = h1.breaks
    h2 = histogram_irregular(x; support=(0.0, 1.0))
    breaks2 = h2.breaks
    @test isapprox(breaks1[1], 0.0; atol=1e-10)
    @test isapprox(breaks1[end], 1.0; atol=1e-10)
    @test isapprox(breaks2[1], 0.0; atol=1e-10)
    @test isapprox(breaks2[end], 1.0; atol=1e-10)
end

@testset "is density" begin
    n = 100
    x = randn(n)
#=     H1 = histogram_regular(x)
    edges1 = H1.edges[1]
    dens1 = H1.weights
    H2 = histogram_irregular(x)
    edges2 = H2.edges[1]
    dens2 = H2.weights
    @test sum(dens1 .* (edges1[2:end] - edges1[1:end-1])) ≈ 1.0
    @test sum(dens2 .* (edges2[2:end] - edges2[1:end-1])) ≈ 1.0 =#

    h1 = histogram_regular(x)
    breaks1 = h1.breaks
    dens1 = h1.density
    h2 = histogram_irregular(x)
    breaks2 = h2.breaks
    dens2 = h2.density
    @test sum(dens1 .* (breaks1[2:end] - breaks1[1:end-1])) ≈ 1.0
    @test sum(dens2 .* (breaks2[2:end] - breaks2[1:end-1])) ≈ 1.0
end

@testset "erroneous kwargs throw error" begin
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
    @test_throws ArgumentError histogram_regular(x; rule=:wand, level=100)
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
    h1 = histogram_irregular(x; rule=:klcv, use_min_length=true)
    h2 = histogram_irregular(x; rule=:l2cv, use_min_length=true)

    min_length = (maximum(x) - minimum(x))*log(n)^(1.5)/n

    #= @test minimum(H1.edges[1][2:end] - H1.edges[1][1:end-1]) ≥ min_length
    @test minimum(H2.edges[1][2:end] - H2.edges[1][1:end-1]) ≥ min_length =#
    @test minimum(h1.breaks[2:end] - h1.breaks[1:end-1]) ≥ min_length
    @test minimum(h2.breaks[2:end] - h2.breaks[1:end-1]) ≥ min_length
end

@testset "throws error misspecified support" begin
    x = [-1.0, 1.0]
    @test_throws DomainError histogram_irregular(x; support=(-0.5, Inf))
    @test_throws DomainError histogram_irregular(x; support=(-Inf, 0.5))
end

@testset "AutomaticHistogram plot and string" begin
    breaks = [0.0, 0.4, 0.6, 1.0]
    counts = [2, 5, 10]
    density = counts ./ ((breaks[2:end] - breaks[1:end-1])*sum(counts))

    h = AutomaticHistogram(breaks, density, counts, :irregular, :right, 1.0)
    @test typeof(plot(h)) == Plots.Plot{Plots.GRBackend}                            # check that Plots extension works

    io = IOBuffer() # just checks that we can call the show method
    show(io, h)
    output = String(take!(io))
    @test typeof(output) == String
end

@testset "AutomaticHistogram fit" begin
    x = randn(10^3)

    @test fit(AutomaticHistogram, x) == histogram_irregular(x)                  # check defaults, irregular
    @test fit(AutomaticHistogram, x; type=:regular) == histogram_regular(x)     # check defaults, regular

    kwargs1 = Dict(:rule => :penb, :grid => :quantile)
    kwargs2 = Dict(:rule => :wand, :scalest => :iqr, :level => 4)

    @test fit(AutomaticHistogram, x; kwargs1...) == histogram_irregular(x; kwargs1...)
    @test fit(AutomaticHistogram, x; kwargs2...) == histogram_regular(x; kwargs2...)

    @test_throws ArgumentError fit(AutomaticHistogram, x; rule=:nonsense)               # test error handling
    @test_throws ArgumentError fit(AutomaticHistogram, x; rule=:l2cv, type=:nonsense)
end

@testset "AutomaticHistogram loglik, logmarglik" begin
    breaks = [0.0, 0.4, 0.6, 1.0]
    counts = [2, 5, 10]
    density = counts ./ ((breaks[2:end] - breaks[1:end-1])*sum(counts))

    @test_throws ArgumentError logmarginallikelihood(AutomaticHistogram(breaks, density, counts, :irregular, :right))
    @test logmarginallikelihood(AutomaticHistogram(breaks, density, counts, :irregular, :right, 1.0)) == logmarginallikelihood(AutomaticHistogram(breaks, density, counts, :irregular, :right), 1.0)
    @test !isnan(loglikelihood(AutomaticHistogram(breaks, density, counts, :irregular, :right)))
end

@testset "AutomaticHistogram extrema" begin
    x = rand(10^3)
    support = (0.0, 1.0)
    h = fit(AutomaticHistogram, x; rule=:aic, support=support) 

    @test isapprox(minimum(h), support[1]; atol=1e-10)
    @test isapprox(maximum(h), support[2]; atol=1e-10)
    @test isapprox(extrema(h)[1], support[1]; atol=1e-10) && isapprox(extrema(h)[2], support[2]; atol=1e-10)
end

@testset "AutomaticHistogram modes" begin
    breaks1 = LinRange(0, 1, 11)
    density1 = [108, 105, 105, 114, 91, 88, 80, 105, 101, 103] / 100.0
    counts1 = [108, 105, 105, 114, 91, 88, 80, 105, 101, 103]
    true_modes1 = [0.05, 0.35, 0.75, 0.95]

    breaks2 = LinRange(0, 1, 11)
    density2 = [95, 84, 101, 101, 102, 106, 100, 104, 104, 103] / 100.0
    counts2 = [95, 84, 101, 101, 102, 106, 100, 104, 104, 103]
    true_modes2 = [0.05, 0.55, 0.8]

    breaks3 = LinRange(0, 1, 11)
    density3 = fill(100, 10) / 100.0
    counts3 = fill(100, 10)
    true_modes3 = [0.5]

    breaks4 = LinRange(0, 1, 11)
    density4 = [1.0]
    counts4 = [1]
    true_modes4 = [0.5]

    @test modes(AutomaticHistogram(breaks1, density1, counts1, :irregular, :right)) == true_modes1
    @test modes(AutomaticHistogram(breaks2, density2, counts2, :irregular, :right)) == true_modes2
    @test modes(AutomaticHistogram(breaks3, density3, counts3, :irregular, :right)) == true_modes3
    @test modes(AutomaticHistogram(breaks4, density4, counts4, :irregular, :right)) == true_modes4
end