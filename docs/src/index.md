# AutoHist.jl Documentation

Fast automatic histogram construction. Supports a plethora of regular and irregular histogram procedures.

## Quick Start
The two main functions exported by this package are `histogram_irregular` and `histogram_regular`, which constructs an irregular or regular histogram with automatic selection of the number of bins based on the sample. The following example shows how to compute and display a regular and an irregular histogram, with an automatic selection of the number of bins.

```@example index; continued=true
using AutoHist, Random, Distributions
x = rand(Xoshiro(1812), Normal(), 10^6)     # simulate some data
h_irr = histogram_irregular(x)              # compute an automatic irregular histogram
h_reg = histogram_regular(x)                # compute an automatic regular histogram
```

Both `histogram_irregular` and `histogram_regular` return a [StatsBase.Histogram](https://juliastats.org/StatsBase.jl/stable/empirical/#StatsBase.Histogram), with weights normalized so that the resulting histograms are probability densities. This allows us to easily plot the two histograms resulting from the above code snippet:

```@example index
using Plots; gr()
# Plot the resulting histograms
p1 = plot(h_irr, xlabel="x", ylabel="Density", label="Irregular", alpha=0.4, color="red")
p2 = plot(h_reg, xlabel="x", label="Regular", alpha=0.4, color="blue")
plot(p1, p2, layout=(1, 2), size=(600, 300))
```


## Supported methods
Both the regular and the irregular procedure support a large number of criteria to select the histogram partition. The keyword argument `rule` controls the criterion used to choose the best partition, and includes the following criteria:

- Regular Histograms:
    - Regular random histogram, "bayes" (default)
    - L2 cross-validation, "l2cv"
    - Kullback-Leibler cross-validation: "klcv"
    - AIC, "aic"
    - BIC, "bic"
    - Birgé and Rozenholc's criterion, "br"
    - Normalized Maximum Likelihood, "nml"
    - Minimum Description Length, "mdl"
    - Sturges' rule, "sturges"
    - Freedman and Diaconis' rule, "fd"
    - Scott's rule, "scott"
    - Wand's rule, "wand"
- Irregular Histograms:
    - Irregular random histogram, "bayes" (default)
    - L2 cross-validation, "l2cv"
    - Kullback-Leibler cross-validation: "klcv"
    - Rozenholc et al. penalty R: "penR"
    - Rozenholc et al. penalty B: "penB"
    - Normalized Maximum Likelihood: "nml"

A more detailed description along with references for each method can be found on the [methods page](methods.md).

Example usage with different rules:
```julia
histogram_irregular(x; rule="penr")
histogram_regular(x; rule="aic")
```

## Features 
In addition to providing automatic histogram construction, this library will at a later point in time include several convenience functions for histograms. These include functions to determine the number and the location of the modes of a histogram, and functions to compute numerical estimation error made with piecewise continuous densities in mind.