# AutoHist.jl Documentation

Fast automatic histogram construction. Supports a plethora of regular and irregular histogram procedures.

## Introduction
Despite being the oldest nonparametric density estimator, the histogram remains widespread in use even to this day. Regrettably, the quality of a histogram density estimate is rather sensitive to the choice of partition used to draw the histogram, which has lead to the development of automatic histogram methods that select the partition based on the sample itself. Unfortunately, most default histogram plotting software only support a few regular automatic histogram procedures, where all the bins are of equal length, and use very simple plug-in rules by default to compute the the number of bins, frequently leading to poor density estimates for non-normal data. Moreover, fast and fully automatic irregular histogram methods are rarely supported by default plotting software, which has prevented their adaptation by practitioners.

The AutoHist.jl package makes it easy to construct both regular and irregular histograms automatically based on a given one-dimensional sample. It currently supports 7 different methods for irregular histograms and 12 criteria for regular histograms from the statistical literature. In addition, the package provides a number of convenience functions for automatic histograms, such as methods for evaluating the histogram probability density function or identifying the location of modes.

## Quick Start
The three main functions exported by this package are [`fit`](@ref), [`histogram_irregular`](@ref) and [`histogram_regular`](@ref), which constructs an irregular or regular histogram with automatic selection of the number of bins based on the sample. The following example shows how to compute and display an irregular and a regular histogram, with an automatic selection of the number of bins.

```@example index; continued=true
using AutoHist, Random, Distributions
x = rand(Xoshiro(1812), Normal(), 10^6)     # simulate some data
h_irr = histogram_irregular(x)              # compute an automatic irregular histogram
h_reg = histogram_regular(x)                # compute an automatic regular histogram
```
Alternatively, irregular and regular automatic histograms can be fitted to data using the `fit` method by controlling the `type` keyword argument.
```@example index; continued=true
h_irr = fit(AutomaticHistogram, x; type=:irregular)  # equivalent to h_irr = histogram_irregular(x)
h_reg = fit(AutomaticHistogram, x; type=:regular)    # equivalent to h_reg = histogram_regular(x)
```

All of the above functions return an object of type [`AutomaticHistogram`](@ref), with weights normalized so that the resulting histograms are probability densities. This type represents the histogram in a similar fashion to [StatsBase.Histogram](https://juliastats.org/StatsBase.jl/stable/empirical/#Histograms), but has more fields to enable the use of several convenience functions.
```@example index
h_irr
```

AutomaticHistogram objects are compatible with Plots.jl, which allows us to easily plot the two histograms resulting from the above code snippet:

```@example index
using Plots; gr()
# Plot the resulting histograms
p_irr = plot(h_irr, xlabel="x", ylabel="Density", label="Irregular", alpha=0.4, color="red")
p_reg = plot(h_reg, xlabel="x", label="Regular", alpha=0.4, color="blue")
plot(p_irr, p_reg, layout=(1, 2), size=(600, 300))
```

## Supported methods
Both the regular and the irregular procedure support a large number of criteria to select the histogram partition. The keyword argument `rule` controls the criterion used to choose the best partition, and includes the following options:

- Regular Histograms:
    - Random regular histogram, :bayes (default)
    - L2 cross-validation, :l2cv
    - Kullback-Leibler cross-validation: :klcv
    - AIC, :aic
    - BIC, :bic
    - Birgé and Rozenholc's criterion, :br
    - Normalized Maximum Likelihood, :nml
    - Minimum Description Length, :mdl
    - Sturges' rule, :sturges
    - Freedman and Diaconis' rule, :fd
    - Scott's rule, :scott
    - Wand's rule, :wand
- Irregular Histograms:
    - Random irregular histogram, :bayes (default)
    - L2 cross-validation, :l2cv
    - Kullback-Leibler cross-validation: :klcv
    - Rozenholc et al. penalty R: :penr
    - Rozenholc et al. penalty B: :penb
    - Rozenholc et al. penalty A: :pena
    - Normalized Maximum Likelihood: :nml

A more detailed description along with references for each method can be found on the [methods page](methods.md).

Example usage with different rules:
```julia
histogram_irregular(x; rule=:penr)
histogram_regular(x; rule=:aic)
```