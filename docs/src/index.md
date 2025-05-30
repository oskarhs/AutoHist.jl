# AutoHist.jl Documentation

Fast automatic histogram construction. Supports a plethora of regular and irregular histogram procedures.

## Quick Start
The two main functions exported by this package are `histogram_irregular` and `histogram_regular`, which constructs an irregular or regular histogram with automatic selection of the number of bins based on the sample.

```julia
julia> using AutoHist, Plots
julia> x = randn(10^6);
julia> H1 = histogram_regular(x);
julia> plot(H1)

julia> H2 = histogram_irregular(x);
julia> plot(H2)
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
- Irregular Histograms:
    - Irregular random histogram, "bayes" (default)
    - L2 cross-validation, "l2cv"
    - Kullback-Leibler cross-validation: "klcv"
    - Rozenholc et al. penalty R: "penR"
    - Rozenholc et al. penalty B: "penB"
    - Normalized Maximum Likelihood: "nml"

```julia
julia> histogram_irregular(x; rule="penR")
julia> histogram_regular(x; rule="aic")
```

WIP: Add another page where each method is described in greater detail.

## Features 
In addition to providing automatic histogram construction, this library will at a later point in time include several convenience functions for histograms. These include functions to determine the number and the location of a histogram, and functions to compute numerical estimation error made with piecewise continuous densities in mind.