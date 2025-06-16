# AutoHist.jl

[![Build Status](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/oskarhs/AutoHist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/oskarhs/AutoHist.jl)


A pure Julia implementation of automatic regular and irregular histogram methods based on maximizing a goodness-of-fit criterion.

## Introduction
Most default histogram plotting software only support regular automatic histogram procedures and use very simple plug-in rules to compute the the number of bins, frequently leading to poor density estimates for non-normal data \[cf. [Birgé and Rozenholc (2006)](#birge2006bins), [Simensen et al. (2025)](#simensen2025random)\]. The purpose of this software package is to offer the user a fast implementation of more sophisticated regular and irregular histogram procedures. Our package supports a variety of methods including those based on asymptotic risk minimization, leave-one-out cross-validiation, penalized maximum likelihood and fully Bayesian approaches.

This module exports three functions that can be used to fit an automatic histogram to a one-dimensional sample, namely `fit`, `histogram_irregular` and `histogram_regular`. A detailed exposition of these functions can be found by typing `?fit`, `?histogram_irregular` and `?histogram_regular` in the repl or in the [API documentation](https://oskarhs.github.io/AutoHist.jl/stable/api/).

## Installation
Installing the package is most easily done via Julia's builtin package manager `Pkg`. This package is part of the Julia general registry, so the installation can be done via the two following lines of code:
```julia
using Pkg
Pkg.add("AutoHist")
```

## Example usage

The following code snippet shows an example where an automatic regular histogram and an automatic irregular histogram are fitted to a normal random sample, and the resulting histograms are plotted.

```julia
julia> using AutoHist, Plots
julia> x = randn(10^6);
julia> h1 = histogram_regular(x);
julia> plot(h1)

julia> h2 = histogram_irregular(x);
julia> plot(h2)

julia> h3 = fit(AutomaticHistogram, x; rule=:wand, scalest=:stdev);
julia> plot(h3)
```

## Supported criteria

The keyword argument `rule` determines the method used to construct the histogram for both of the histogram functions. The rule used to construct the histogram can be changed by setting `rule` equal to a symbol indicating the method to be used, e.g. `rule=:aic`.

The default method is the Bayesian approach of [Simensen et al. (2025)](#simensen2025random), corresponding to keyword `rule=:bayes`.
A detailed description of the supported methods can be found in the [methods documentation](https://oskarhs.github.io/AutoHist.jl/stable/methods/).

## Implementation

### Irregular histograms
Our implementation uses the dynamical programming algorithm of [Kanazawa (1988)](#kanazawa1988optimal) together with the greedy search heuristic of [Rozenholc et al. (2010)](#rozenholc2010combining) to build a histogram quickly, making this package an excellent option for histogram construction for large data sets.

### Regular histograms
For regular histograms, we provide a multithreaded implementation, which will be automatically used when Julia is launched with more than one thread.

## To do
- Add a new type for automatically constructed histograms that allows for the evaluation of loglikelihoods and similar quantities of interest.

## References
<a name="simensen2025random"></a> Simensen, O. H., Christensen, D. & Hjort, N. L. (2025). Random Irregular Histograms. _arXiv preprint_. doi: [10.48550/ARXIV.2505.22034](https://doi.org/10.48550/ARXIV.2505.22034)

<a name="rozenholc2010combining"></a> Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. _Computational Statistics & Data Analysis_, **54**, 3313–3323. doi: [10.1016/j.csda.2010.04.021](https://doi.org/10.1016/j.csda.2010.04.021)

<a name="kanazawa1988optimal"></a> Kanazawa, Y. (1988). An optimal variable cell histogram. _Communications in Statistics-Theory and Methods_, **17**, 1401–1422. doi: [10.1080/03610928808829688](https://doi.org/10.1080/03610928808829688)

<a name="birge2006bins"></a> Birgé, L., & Rozenholc, Y. (2006). How many bins should be put in a regular histogram. _ESAIM: Probability and Statistics_, **10**, 24–45. doi: [10.1051/ps:2006001](https://doi.org/10.1051/ps:2006001)