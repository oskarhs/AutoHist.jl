# Examples

The following document illustrates the use of AutoHist.jl through examples. In particular, we showcase some of the relative advantages and disadvantages of regular and irregular histogram procedures.

## Estimating the Chi-square probability density

We start by considering an example with some simulated data from the [LogNormal-distribution](https://en.wikipedia.org/wiki/Log-normal_distribution). To start, we fit a regular histogram to the data, using the approach of [Birgé and Rozenholc (2006)](https://doi.org/10.1016/j.csda.2010.04.021), which corresponds to `rule=:br`.
```@example LogNormal; continued=true
using AutoHist, Random, Distributions
x = rand(Xoshiro(1812), LogNormal(), 10^4)
h1 = fit(AutomaticHistogram, x; rule=:br)
```
Alternatively, since the standard LogNormal pdf has known support ``[0,\infty)``, we can incorporate this knowledge in our histogram estimate through the `support` keyword.
```@example LogNormal; continued = true
h2 = fit(AutomaticHistogram, x; rule=:br, support=(0.0, Inf))
```
In this case, specifying the support makes little difference since the sample size is quite large, and the standard lognormal distribution assigns sufficient mass near ``0``.

The standard LogNormal is quite challenging to estimate well using a regular histogram procedure due to its heavy tails. These two factors make irregular methods an appealing alternative in this case. Here, we use the penalized log-likelihood approach from [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) with penalty `:penr` and a data-based grid to construct the histogram.
```@example LogNormal; continued = true
h3 = fit(AutomaticHistogram, x; rule=:penr, grid=:data, support=(0.0, Inf))
```

To compare the two approaches, we can plot the resulting histograms along with the true density:
```@example LogNormal
using Distributions, Plots; gr() # hide
t = LinRange(0.0, 8.5, 1000) # hide
p = plot(t, pdf(LogNormal(), t), xlabel="x", label="True density", color="black", lwd=1.5) # hide
p1 = plot(p, h2, ylabel="Density", label="Regular", alpha=0.4, color="red") # hide
xlims!(p1, -0.5, 8.5) # hide
p2 = plot(p, h3, label="Irregular", alpha=0.4, color="blue") # hide
xlims!(p2, -0.5, 8.5) # hide
plot(p1, p2, layout=(1, 2), size=(600, 300)) # hide
```

## Mode hunting
