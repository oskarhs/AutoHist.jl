# AutoHist.jl

[![Build Status](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/oskarhs/AutoHist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/oskarhs/AutoHist.jl)


A pure Julia implementation of automatic regular and irregular histogram methods based on maximizing a goodness-of-fit criterion.

## Introduction
Most default histogram plotting software only support a few regular automatic histogram procedures and use very simple plug-in rules by default to compute the the number of bins, frequently leading to poor density estimates for non-normal data \[cf. [Birgé and Rozenholc (2006)](#birge2006bins), [Simensen et al. (2025)](#simensen2025random)\]. The purpose of this software package is to offer the user a fast and simple-to-use implementation of more sophisticated regular and irregular histogram procedures. Our package supports a variety of methods including those based on asymptotic risk minimization, leave-one-out cross-validiation, penalized maximum likelihood and fully Bayesian approaches.

## Installation
Installing the package is most easily done via Julia's builtin package manager `Pkg`. This package is part of the Julia general registry, so the installation can be done via the two following lines of code:
```julia
using Pkg
Pkg.add("AutoHist")
```

## Quick start

To get started, we illustrate the basic use of the package by fitting histograms to a normal random sample with an automatic selection of the histogram partition.

```julia
using AutoHist
x = randn(10^6)
h1 = fit(AutomaticHistogram, x)
```

Several rules are available to select the histogram partition, and can be controlled through the use of the `rule` keyword argument. For instance, a regular histogram based on maximizing the AIC can be fit as follows:
```julia
h2 = fit(AutomaticHistogram, x; rule=:aic)
```

Alternatively, we can fit regular and irregular histograms by calling the functions `histogram_regular` or `histogram_regular`.
```julia
h3 = histogram_regular(x; rule=:wand, scalest=:iqr)
h4 = histogram_irregular(x; rule=:penr)
```

A detailed exposition of the keyword arguments passed to each of these functions can be found by typing `?fit`, `?histogram_irregular` and `?histogram_regular` in the repl or in the [API documentation](https://oskarhs.github.io/AutoHist.jl/stable/api/).

The fitted histograms can be displayed through the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) or [Makie.jl](https://github.com/MakieOrg/Makie.jl) packages as follows:

```julia
import Plots, CairoMakie
Plots.plot(h)
Makie.plot(h)
```

## Supported criteria

The keyword argument `rule` determines the method used to construct the histogram for both of the histogram functions. The rule used to construct the histogram can be changed by setting `rule` equal to a symbol indicating the method to be used, e.g. `rule=:aic`.

The default method is the Bayesian approach of [Simensen et al. (2025)](#simensen2025random), corresponding to keyword `rule=:bayes`.
A detailed description of the supported methods can be found in the [methods documentation](https://oskarhs.github.io/AutoHist.jl/stable/methods/).

## References
<a name="simensen2025random"></a> Simensen, O. H., Christensen, D. & Hjort, N. L. (2025). Random Irregular Histograms. _arXiv preprint_. doi: [10.48550/ARXIV.2505.22034](https://doi.org/10.48550/ARXIV.2505.22034)

<a name="rozenholc2010combining"></a> Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. _Computational Statistics & Data Analysis_, **54**, 3313–3323. doi: [10.1016/j.csda.2010.04.021](https://doi.org/10.1016/j.csda.2010.04.021)

<a name="kanazawa1988optimal"></a> Kanazawa, Y. (1988). An optimal variable cell histogram. _Communications in Statistics-Theory and Methods_, **17**, 1401–1422. doi: [10.1080/03610928808829688](https://doi.org/10.1080/03610928808829688)

<a name="birge2006bins"></a> Birgé, L., & Rozenholc, Y. (2006). How many bins should be put in a regular histogram. _ESAIM: Probability and Statistics_, **10**, 24–45. doi: [10.1051/ps:2006001](https://doi.org/10.1051/ps:2006001)