import Plots; Plots.gr() # avoid namespace pollution
using AutoHist, Distributions
using Test
import StatsBase: Histogram, fit

@testset "return type regular" begin
    x = collect(LinRange(0,1,11))

    for rule in [:bayes, :knuth, :aic, :bic, :br, :mdl, :nml, :l2cv, :klcv,
                 :sturges, :fd, :scott]
        for closed in [:left, :right]
            h = histogram_regular(x; rule=rule, closed=closed)
            @test typeof(h) <: AutomaticHistogram
        end
    end
    for scalest in [:minim, :iqr, :stdev]
        for level in [0,1,2,3,4,5]
            for closed in [:left, :right]
                h = histogram_regular(x; rule=:wand, closed=closed, scalest=scalest, level=level)
                @test typeof(h) <: AutomaticHistogram
            end
        end
    end
    h = histogram_regular(x; rule=:wand)
    @test typeof(h) <: AutomaticHistogram
end

@testset "return type irregular" begin
    x = collect(LinRange(0,1,11))

    for rule in [:pena, :penb, :penr, :bayes, :klcv, :l2cv, :nml]
        h = histogram_irregular(x; rule=rule)
        @test typeof(h) <: AutomaticHistogram
    end
    for grid in [:regular, :data, :quantile] # test grid, right-left open interval combinations
        for closed in [:left, :right]
            h = histogram_irregular(x; grid=grid, closed=closed)
            @test typeof(h) <: AutomaticHistogram
        end
    end
    x = collect(LinRange(0,1,11))
    h = histogram_irregular(x; alg=DP(greedy=false)) # check that greedy == false works
    @test typeof(h) <: AutomaticHistogram
end

@testset "left open and right open intervals" begin
    x = collect(LinRange(0,1,11))
    h1 = histogram_regular(x; closed=:right)
    h2 = histogram_regular(x; closed=:left)
    h3 = histogram_irregular(x; closed=:right)
    h4 = histogram_irregular(x; closed=:left)

    @test h1.closed == :right
    @test h2.closed == :left
    @test h3.closed == :right
    @test h4.closed == :left
end

@testset "algorithms" begin
    x = collect(LinRange(0,1,11))

    for alg in [DP(), GPDP()]
        for rule in [:penr, :bayes]
            h = histogram_irregular(x; rule=rule, alg=alg)
            @test typeof(h) <: AutomaticHistogram
        end
    end
    for rule in [:klcv, :l2cv]
        h = histogram_irregular(x; rule=rule, alg=DP())
        @test typeof(h) <: AutomaticHistogram
    end

    for alg in [DP(gr_maxbins=5), GPDP(gr_maxbins=5, max_cand=10)]
        h = histogram_irregular(x; alg=alg)
        @test typeof(h) <: AutomaticHistogram
    end

    # Error handling
    @test_throws ArgumentError DP(gr_maxbins = :nonsense)
    @test_throws ArgumentError GPDP(gr_maxbins = :nonsense)
    @test_throws DomainError DP(gr_maxbins = -1)
    @test_throws DomainError GPDP(gr_maxbins = -1)
    @test_throws ArgumentError GPDP(max_cand = :nonsense)
    @test_throws DomainError GPDP(max_cand = -1)

    @test_throws ArgumentError histogram_irregular(x; rule=:l2cv, alg=GPDP())
end

@testset "estimated support" begin
    n = 100
    x = [-5.0, 4.5]
    h1 = histogram_regular(x)
    breaks1 = h1.breaks
    h2 = histogram_irregular(x)
    breaks2 = h2.breaks
    @test isapprox(breaks1[1], minimum(x); atol=1e-10)
    @test isapprox(breaks1[end], maximum(x); atol=1e-10)
    @test isapprox(breaks2[1], minimum(x); atol=1e-10)
    @test isapprox(breaks2[end], maximum(x); atol=1e-10)
end

@testset "given support" begin
    n = 100
    x = rand(n)

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
    @test_throws DomainError histogram_regular(x; maxbins=-100)
end

@testset "a as function" begin
    x = randn(10^3)

    @test histogram_regular(x; a = k->8.0) == histogram_regular(x; a=8.0)
    @test histogram_regular(x; a = k->8.0) ≈ histogram_regular(x; a=8.0)
end

@testset "min length" begin
    n = 10^3
    x = randn(n)
    h1 = histogram_irregular(x; rule=:klcv, use_min_length=true)
    h2 = histogram_irregular(x; rule=:l2cv, use_min_length=true)
    min_length = (maximum(x) - minimum(x))*log(n)^(1.5)/n

    @test minimum(h1.breaks[2:end] - h1.breaks[1:end-1]) ≥ min_length
    @test minimum(h2.breaks[2:end] - h2.breaks[1:end-1]) ≥ min_length
end

@testset "throws error misspecified support" begin
    x = [-1.0, 1.0]
    @test_throws DomainError histogram_irregular(x; support=(-0.5, Inf))
    @test_throws DomainError histogram_irregular(x; support=(-Inf, 0.5))
    @test_throws DomainError histogram_regular(x; support=(-0.5, Inf))
    @test_throws DomainError histogram_regular(x; support=(-Inf, 0.5))
end

@testset "AutomaticHistogram Plots and string" begin
    breaks = [0.0, 0.4, 0.6, 1.0]
    counts = [2, 5, 10]
    density = counts ./ ((breaks[2:end] - breaks[1:end-1])*sum(counts))

    h = AutomaticHistogram(breaks, density, counts, :irregular, :right, 1.0)
    @test typeof(Plots.plot(h)) == Plots.Plot{Plots.GRBackend}    # check that Plots extension works
    
    @static if VERSION ≥ v"1.10"
        import Makie, CairoMakie

        @test typeof(Makie.plot(h)) == Makie.FigureAxisPlot # check that Makie extension works
        @test typeof(Makie.plot!(h)) == Makie.PlotList{Tuple{Makie.PlotSpec}}
        @test typeof(Makie.barplot(h)) == Makie.FigureAxisPlot
        @test typeof(Makie.hist(h)) == Makie.FigureAxisPlot
        @test typeof(Makie.stephist(h)) == Makie.FigureAxisPlot

        F = Makie.Figure(); ax = Makie.Axis(F[1,1])
        Makie.plot!(ax, h)
        @test length(ax.scene.plots) == 1 && typeof(ax.scene.plots[1]) == Makie.PlotList{Tuple{Makie.PlotSpec}}
    end

    io = IOBuffer() # just checks that we can call the show method
    show(io, h)
    output = String(take!(io))
    @test typeof(output) == String
end

@testset "Makie" begin
    breaks = [0.0, 0.4, 0.6, 1.0]
    counts = [2, 5, 10]
    density = counts ./ ((breaks[2:end] - breaks[1:end-1])*sum(counts))
    h = AutomaticHistogram(breaks, density, counts, :irregular, :right, 1.0)
    
    @static if VERSION ≥ v"1.10"
        import Makie, CairoMakie

        @test typeof(Makie.plot(h)) == Makie.FigureAxisPlot # check that Makie extension works
        @test typeof(Makie.plot!(h)) == Makie.PlotList{Tuple{Makie.PlotSpec}}
        @test typeof(Makie.barplot(h)) == Makie.FigureAxisPlot
        @test typeof(Makie.hist(h)) == Makie.FigureAxisPlot
        @test typeof(Makie.stephist(h)) == Makie.FigureAxisPlot

        F = Makie.Figure(); ax = Makie.Axis(F[1,1])
        Makie.plot!(ax, h)
        @test length(ax.scene.plots) == 1 && typeof(ax.scene.plots[1]) == Makie.PlotList{Tuple{Makie.PlotSpec}}
    end
end

@testset "AutomaticHistogram approx different dims" begin
    breaks1 = LinRange(0, 1, 11)
    density1 = [95, 84, 101, 101, 102, 106, 100, 104, 104, 103] / 100.0
    counts1 = [95, 84, 101, 101, 102, 106, 100, 104, 104, 103]

    breaks2 = LinRange(0, 1, 6)
    density2 = fill(100, 5) / 100.0
    counts2 = fill(100, 5)

    @test !(AutomaticHistogram(breaks1, density1, counts1, :irregular, :right) ≈ AutomaticHistogram(breaks2, density2, counts2, :irregular, :right))
    @test !(AutomaticHistogram(breaks1, density1, counts1, :irregular, :right, 1.0) ≈ AutomaticHistogram(breaks2, density2, counts2, :irregular, :right))
    @test !(AutomaticHistogram(breaks1, density1, counts1, :irregular, :right) ≈ AutomaticHistogram(breaks2, density2, counts2, :irregular, :right, 1.0))
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
    @test_throws ArgumentError fit(AutomaticHistogram, x; type=:nonsense)
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
    @test AutoHist.insupport(h, 0.5)
    @test !AutoHist.insupport(h, -0.5)
    @test !AutoHist.insupport(h, 1.5)
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

    @test peaks(AutomaticHistogram(breaks1, density1, counts1, :regular, :right)) == true_modes1
    @test peaks(AutomaticHistogram(breaks2, density2, counts2, :regular, :right)) == true_modes2
    @test peaks(AutomaticHistogram(breaks3, density3, counts3, :regular, :right)) == true_modes3
    @test peaks(AutomaticHistogram(breaks4, density4, counts4, :regular, :right)) == true_modes4
end

@testset "AutomaticHistogram pdf" begin
    breaks1 = LinRange(0, 1, 11)
    density1 = [1.08, 1.05, 1.05, 1.14, 0.91, 0.88, 0.80, 1.05, 1.01, 1.03]
    counts1 = [108, 105, 105, 114, 91, 88, 80, 105, 101, 103]

    breaks2 = [0.0, 0.2, 0.8, 0.9, 1.0]
    density2 = [0.945, 1.045, 0.880, 0.960]
    counts2 = [189, 627, 88, 96]

    @test AutoHist.pdf(AutomaticHistogram(breaks1, density1, counts1, :regular, :right), -0.1) == 0.0    # test values outside of support
    @test AutoHist.pdf(AutomaticHistogram(breaks2, density2, counts2, :irregular, :right), 1.2) == 0.0

    for j in 1:9    # test at each boundary that closed=:right and closed=:left behave as expected
        @test AutoHist.pdf(AutomaticHistogram(breaks1, density1, counts1, :regular, :right), breaks1[j+1]) == density1[j]
        @test AutoHist.pdf(AutomaticHistogram(breaks1, density1, counts1, :regular, :left), breaks1[j+1]) == density1[j+1]
    end

    for j in 1:3    # test at each boundary that closed=:right and closed=:left behave as expected
        @test AutoHist.pdf(AutomaticHistogram(breaks2, density2, counts2, :irregular, :right), breaks2[j+1]) == density2[j]
        @test AutoHist.pdf(AutomaticHistogram(breaks2, density2, counts2, :irregular, :left), breaks2[j+1]) == density2[j+1]
    end

    # Test the broadcasted versions as well
    @test AutoHist.pdf.(AutomaticHistogram(breaks1, density1, counts1, :regular, :right), [0.1]) == [1.08]
    @test AutoHist.pdf.(AutomaticHistogram(breaks2, density2, counts2, :irregular, :right), [0.2]) == [0.945]
end

@testset "AutomaticHistogram distance" begin
    breaks = LinRange(0, 1, 11)
    density = [1.08, 1.05, 1.05, 1.14, 0.91, 0.88, 0.80, 1.05, 1.01, 1.03]
    counts = [108, 105, 105, 114, 91, 88, 80, 105, 101, 103]

    h = AutomaticHistogram(breaks, density, counts, :regular, :right)

    for dist in [:iae, :ise, :hell, :sup, :kl]
        @test distance(h, h, dist) == 0.0
    end
    @test distance(h, h, :lp; p=3.0) == 0.0
    @test distance(h, h, :lp; p=Inf) == 0.0

    @test_throws ArgumentError distance(h, h, :nonsense) # error handling
    @test_throws DomainError distance(h, h, :lp; p=-1.0)

    @test distance(
        AutomaticHistogram([-1.0, 0.0], [1.0], [1], :regular, :right),
        AutomaticHistogram([0.0, 1.0], [1.0], [1], :regular, :right),
        :iae
    ) == distance(
        AutomaticHistogram([0.0, 1.0], [1.0], [1], :regular, :right),
        AutomaticHistogram([-1.0, 0.0], [1.0], [1], :regular, :right),
        :iae
    ) == 2.0
    @test Inf == distance(
        AutomaticHistogram([-1.0, 0.0], [1.0], [1], :regular, :right),
        AutomaticHistogram([0.0, 1.0], [1.0], [1], :regular, :right),
        :kl
    )
end