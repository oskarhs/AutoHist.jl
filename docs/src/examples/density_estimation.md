# Density estimation

The following document illustrates the use of AutoHist.jl through examples from the world of density estimation. In particular, we showcase some of the relative advantages and disadvantages of regular and irregular histogram procedures.

## Estimating the LogNormal probability density

We start by considering an example with some simulated data from the [LogNormal-distribution](https://en.wikipedia.org/wiki/Log-normal_distribution). To start, we fit a regular histogram to the data, using the approach of [Birgé and Rozenholc (2006)](https://doi.org/10.1016/j.csda.2010.04.021), which corresponds to `rule=BR()`.
```@example LogNormal; continued=true
using AutoHist, Random, Distributions
x = rand(Xoshiro(1812), LogNormal(), 10^4)
h1 = fit(AutomaticHistogram, x, BR())
```
Alternatively, since the standard LogNormal pdf has known support ``[0,\infty)``, we can incorporate this knowledge in our histogram estimate through the `support` keyword.
```@example LogNormal; continued = true
h2 = fit(AutomaticHistogram, x, BR(); support=(0.0, Inf))
```
To quantify the difference of using the correct, known support in this case, we compute the integrated absolute error between the two densities in the following code snippet, which is given by ``\int |f(x) - g(x)|\text{d}x``.
```@example LogNormal
distance(h1, h2, :iae)
```
The resulting ``L_1`` distance of ``0.042`` indicates that the new bin origin at ``0`` has a moderate effect on the resulting density estimate.

The standard LogNormal is quite challenging to estimate well using a regular histogram procedure due to its heavy tails. These two factors make irregular methods an appealing alternative in this case. Here, we use the penalized log-likelihood approach from [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) with penalty R and a data-based grid to construct the histogram, (implemented in AutoHist.jl via `rule = RMG_penR()`).
```@example LogNormal; continued = true
h3 = fit(AutomaticHistogram, x, RMG_penR(grid=:data))
```

To compare the two approaches, we can plot the resulting histograms along with the true density:
```@example LogNormal
using Distributions, Plots; gr() # hide
t = LinRange(0.0, 8.5, 1000) # hide
p = plot(t, pdf(LogNormal(), t), xlabel="x", label="True density", color="blue", lw=2.0, linestyle=:dash) # hide
p1 = plot(p, h2, ylabel="Density", label="Regular", alpha=0.4, color="red") # hide
xlims!(p1, -0.5, 8.5) # hide
p2 = plot(p, h3, label="Irregular", alpha=0.4, color="black") # hide
xlims!(p2, -0.5, 8.5) # hide
plot(p1, p2, layout=(1, 2), size=(670, 320)) # hide
```
The irregular procedure selects smaller bin widths near the origin, reflecting the fact that the LogNormal density is rapidly changing in this area. On the other hand, the bin widths are larger in the flatter region in the right tail of the density. Both histogram procedures provide quite reasonable estimates of the density, owing to the fairly large sample size.

## Mode hunting
In an exploratory data analysis setting, identifying key features in a dataset such as modes is frequently of great interest to statisticians and practicioners alike. Unfortunately, most popular regular histogram have been designed with good performance in terms of statistical risk with respect to classical, integral-based loss functions, which typically results in a large amount of spurious histogram modes in regions where the true density is flat and in the tails of the density [(Scott, 1992)](https://doi.org/10.1002/9780470316849). In practice, this means that a data-analyst must use subjective judgement to infer whether a regular histogram estimate is indicative of a mode being present or not. If the presence of a mode is deemed likely, subjective visual smoothing is typically required to get a rough idea of its location. In constrast, some irregular histogram procedures have been shown empirically to perform quite well with regard to automatic mode detection in cases where the true density has a smallish amount of well-separated modes, see e.g. [Davies et al. (2009)](https://doi.org/10.1051/ps:2008005); [Li et al. (2020)](https://doi.org/10.1093/biomet/asz081); [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).

To illustrate the advantage of irregular histograms when it comes to mode identification, we will consider an example where we attempt to locate the modes of the Harp density of [Li et al. (2020)](https://doi.org/10.1093/biomet/asz081), plotted below.

```@example Harp; continued=true
module TestDistributions # hide

using Distributions # hide
using Random # hide
using Plots # hide

import Distributions: pdf, cdf # hide
import Random: rand # hide

export pdf, cdf, rand, peaks, plot_test_distribution, Harp # hide

struct Harp <: ContinuousUnivariateDistribution end # hide
function pdf(d::Harp, x::Real) # hide
    means = [0.0, 5.0, 15.0, 30.0, 60.0] # hide
    sds = [0.5, 1.0, 2.0, 4.0, 8.0] # hide
    dens = 0.0 # hide
    for j = 1:5 # hide
        dens = dens + 0.2 * pdf(Normal(means[j], sds[j]), x) # hide
    end # hide
    return dens # hide
end # hide
function cdf(d::Harp, x::Real) # hide
    means = [0.0, 5.0, 15.0, 30.0, 60.0] # hide
    sds = [0.5, 1.0, 2.0, 4.0, 8.0] # hide
    cum = 0.0 # hide
    for j = 1:5 # hide
        cum += 0.2 * cdf(Normal(means[j], sds[j]), x) # hide
    end # hide
    return cum # hide
end # hide
function rand(rng::AbstractRNG, d::Harp) # hide
    means = [0.0, 5.0, 15.0, 30.0, 60.0] # hide
    sds = [0.5, 1.0, 2.0, 4.0, 8.0] # hide
    j = rand(rng, DiscreteUniform(1, 5)) # hide
    return rand(rng, Normal(means[j], sds[j])) # hide
end # hide
function peaks(d::Harp) # hide
    return Float64[0.0, 5.0, 15.001837158203125, 30.003319148997427, 60.000555419921874] # hide
end # hide
function plot_test_distribution(d::Harp) # hide
    if d != Harp() # hide
        throw(ArgumentError("Supplied distribution $d is not supported.")) # hide
    end # hide
    dom = LinRange(-5.0, 80.0, 5000) # hide
    x_lims = [-5.0, 80.0] # hide
    p = plot() # hide
    t = dom # hide
    plot!(p, t, pdf.(d, t), label="", ylims=(0.0, Inf), titlefontsize=10, color="blue", linestyle=:dash) # hide
    xlims!(p, x_lims...) # hide
    return p # hide
end # hide
end # hide
```

```@example Harp
using Plots; gr() # hide
import .TestDistributions as TD # hide
p = TD.plot_test_distribution(TD.Harp()) # hide
plot!(p, title="", xlab="x", ylab="Density") # hide
plot(p, size=(650, 310)) # hide
```
The Harp density has 5 peaks with varying degrees of sharpness, with the location of each mode becoming gradually more difficult to locate in absolute terms as we move rightward on the ``x``-axis. In the numerical experiment to follow, we generate a random sample of size ``n = 5000`` from the Harp density, and fit both an irregular and a regular histogram to the data. Motivated by the results of the simulation studies in [Davies et al. (2009)](https://doi.org/10.1051/ps:2008005) and [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034), we have used the BIC criterion to draw the regular histogram, and the random irregular histogram method (`rule=RIH()`), as both of these have been shown to perform relatively well compared to their respective regular/irregular counterparts in automatic mode identification. To access the Harp distribution, we use the TestDistributions package, which can be found in [the following github repository](https://github.com/oskarhs/Random-Histograms---Paper/tree/main/TestDistributions).

```@example Harp; continued = true
import .TestDistributions as TD
using AutoHist, Distributions, Random
x = rand(Xoshiro(1812), TD.Harp(), 5000)

h_reg = fit(AutomaticHistogram, x, BIC())
h_irr = fit(AutomaticHistogram, x, RIH())
```

We can now add the two fitted histograms along with the location of their modes to the plot of the Harp density above:
```@example Harp
using Plots; gr() # hide
p = TD.plot_test_distribution(TD.Harp()) # hide
plot!(p, title="", xlab="x", ylab="Density") # hide
p1 = plot(p, h_reg, title="Regular", color="red", alpha=0.4, label="", ylims=[-0.015, 0.17]) # hide
scatter!(p1, peaks(h_reg), fill(-0.008, length(h_reg)), color="red", label="", ms=3) # hide
p2 = plot(p, h_irr, title="Irregular", color="black", alpha=0.4, label="", ylims=[-0.015, 0.17]) # hide
scatter!(p2, peaks(h_irr), fill(-0.008, length(h_irr)), color="black", label="", ms=3) # hide
plot(p1, p2, layout=(2, 1), size=(670, 650)) # hide
```

The spatial inhomogeneity of the Harp density causes some trouble for the regular histogram, as the use of a global bin width leads to undersmoothing near the sharpest mode and oversmoothing near the flattest mode and in the tails. In particular, the regular estimate introduces many spurious modes near the rightmost mode of the true density, making it difficult to exactly infer the location and the number of modes in the region surrounding it. In contrast, the irregular histogram estimate is able to correctly infer the number of modes in the density automatically, with no subjective judgements required. The good mode-finding behavior is in this case enabled by the fact that the bin widths are able adapt to the local smoothness of the true density, resulting in more smoothing near modes. To inspect the accuracy of the inferred modes, we print the locations of the true modes and the peaks of the irregular histogram.

```@example Harp
println("True modes: ", TD.Harp() |> TD.peaks |> x -> round.(x; digits=2))
println("AutoHist modes: ", h_irr |> peaks |> x -> round.(x; digits=2))
```
We observe that all of the histogram peaks are quite close to a corresponding true mode, especially when taking the sharpness of each individual peak into account.