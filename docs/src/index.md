# AutoHist.jl Documentation

Fast automatic histogram construction. Supports a plethora of regular and irregular histogram procedures.

## Introduction
Despite being the oldest nonparametric density estimator, the histogram remains widespread in use even to this day. Regrettably, the quality of a histogram density estimate is rather sensitive to the choice of partition used to draw the histogram, which has lead to the development of automatic histogram methods that select the partition based on the sample itself. Unfortunately, most default histogram plotting software only support a few regular automatic histogram procedures, where all the bins are of equal length, and use very simple plug-in rules by default to compute the the number of bins, frequently leading to poor density estimates for non-normal data. Moreover, fast and fully automatic irregular histogram methods are rarely supported by default plotting software, which has prevented their adaptation by practitioners.

The AutoHist.jl package makes it easy to construct both regular and irregular histograms automatically based on a given one-dimensional sample. It currently supports 7 different methods for irregular histograms and 12 criteria for regular histograms from the statistical literature. In addition, the package provides a number of convenience functions for automatic histograms, such as methods for evaluating the histogram probability density function or identifying the location of modes.

## Quick Start
The main functions exported by this package are [`fit`](@ref) and [`autohist`](@ref), which allows the user to fit a histogram to 1-dimensional data with an automatic and data-based choice of bins. The following short example shows how the `fit` method is used in practice

```@example index; continued=true
using AutoHist, Random, Distributions
x = rand(Xoshiro(1812), Normal(), 10^6) # simulate some data
h_irr = fit(AutomaticHistogram, x)      # compute an automatic irregular histogram (default method)
```
The third argument to `fit` controls the rule used to select the histogram partition, with the default being the [`RIH`](@ref) method. To fit an automatic histogram with a specific rule, e.g. Knuth's rule, all we have to do is change the value of the `rule` argument.
```@example index; continued=true
h_reg = fit(AutomaticHistogram, x, Knuth()) # compute an automatic regular histogram
```
The above calls to `fit` returns an object of type [`AutomaticHistogram`](@ref), with weights normalized so that the resulting histograms are probability densities. This type represents the histogram in a similar fashion to [StatsBase.Histogram](https://juliastats.org/StatsBase.jl/stable/empirical/#Histograms), but has more fields to enable the use of several convenience functions.
```@example index
h_irr
```
Alternatively, an automatic histogram can be fitted to data through the `autohist` method, which served as an alias for `fit(AutomaticHistogram, x, rule; kwargs...)`.
```@example index; continued=true
h_reg = autohist(x, Knuth())     # equivalent to fit(AutomaticHistogram, x, Knuth())
```

AutomaticHistogram objects are compatible with [Plots.jl](https://github.com/JuliaPlots/Plots.jl), which allows us to easily plot the two histograms resulting from the above code snippet via e.g. `Plots.plot(h_irr)`. To show both histograms side by side, we can create a plot as follows:

```@example index
import Plots; Plots.gr()
p_irr = Plots.plot(h_irr, xlabel="x", ylabel="Density", title="Irregular", alpha=0.4, color="black", label="")
p_reg = Plots.plot(h_reg, xlabel="x", title="Regular", alpha=0.4, color="red", label="")
Plots.plot(p_irr, p_reg, layout=(1, 2), size=(670, 320))
```

Alternatively, [Makie.jl](https://github.com/MakieOrg/Makie.jl) can also be used to make graphical displays of the fitted histograms via e.g. `Makie.plot(h_irr)`. To produce a plot similar to the above display, we may for instance do the following:
```@example index
import CairoMakie, Makie # using the CairoMakie backend
fig = Makie.Figure(size=(670, 320))
ax1 = Makie.Axis(fig[1, 1], title="Irregular", xlabel="x", ylabel="Density")
ax2 = Makie.Axis(fig[1, 2], title="Regular", xlabel="x")
p_irr = Makie.plot!(ax1, h_irr, alpha=0.4, color="black")
p_reg = Makie.plot!(ax2, h_reg, alpha=0.4, color="red")
fig
```

## Supported methods
Both the regular and the irregular procedure support a large number of criteria to select the histogram partition. The keyword argument `rule` controls the criterion used to choose the best partition, and includes the following options:

- Regular Histograms:
    - Knuth's rule, Random regular histogram, [`RRH`](@ref)
    - L2 cross-validation, [`L2CV_R`](@ref)
    - Kullback-Leibler cross-validation, [`KLCV_R`](@ref)
    - Akaike's information criterion, [`AIC`](@ref)
    - The Bayesian information criterion, [`BIC`](@ref)
    - Birgé and Rozenholc's criterion, [`BR`](@ref)
    - Normalized Maximum Likelihood, [`NML_R`](@ref)
    - Minimum Description Length. [`MDL`](@ref)
    - Sturges' rule, [`Sturges`](@ref)
    - Freedman and Diaconis' rule, [`FD`](@ref)
    - Scott's rule, [`Scott`](@ref)
    - Wand's rule, [`Wand`](@ref)
- Irregular Histograms:
    - Random irregular histogram, [`RIH`](@ref)
    - L2 cross-validation, [`L2CV_I`](@ref)
    - Kullback-Leibler cross-validation, [`KLCV_I`](@ref)
    - Rozenholc et al. penalty A, [`RMG_penA`](@ref)
    - Rozenholc et al. penalty B, [`RMG_penB`](@ref)
    - Rozenholc et al. penalty R, [`RMG_penR`](@ref)
    - Normalized Maximum Likelihood, [`NML_I`](@ref)
    - Bayesian Blocks [`BayesBlocks`](@ref)

A description of each method along with references for each method can be found on the [methods page](methods.md).