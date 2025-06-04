# AutoHist.jl

[![Build Status](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/oskarhs/AutoHist.jl/actions/workflows/CI.yml?query=branch%3Amain) [![codecov](https://codecov.io/gh/oskarhs/AutoHist.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/oskarhs/AutoHist.jl)


A pure Julia implementation of automatic regular and irregular histogram methods based on maximizing a goodness-of-fit criterion.

## Introduction
Most default histogram plotting software only support regular automatic histogram procedures and use very simple plug-in rules to compute the the number of bins, frequently leading to poor density estimates for non-normal data \[cf. [Birgé and Rozenholc (2006)](#birge2006bins), [Simensen et al. (2025)](#simensen2025random)\]. The purpose of this software package is to offer the user a fast implementation of more sophisticated regular and irregular histogram procedures. Our package supports a variety of methods including those based on leave-one-out cross-validiation, penalized maximum likelihood and fully Bayesian approaches.

This module exports the two functions `histogram_regular` and `histogram_irregular`, offering automatic histogram construction for 1-dimensional samples. A detailed exposition of all keyword arguments can be found by typing `?histogram_regular` and `?histogram_irregular` in the repl.

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
julia> H1 = histogram_regular(x);
julia> plot(H1)

julia> H2 = histogram_irregular(x);
julia> plot(H2)
```

## Supported criteria

The keyword argument `rule` determines the method used to construct the histogram for both of the histogram functions. The rule used to construct the histogram can be changed by setting `rule` equal to a string indicating the method to be used.

The default method is the Bayesian approach of [Simensen et al. (2025)](#simensen2025random), corresponding to keyword `rule="bayes"`.

A list of the supported methods, along with their corresponding keywords can be found below. A detailed description of the supported methods can be found in the [methods documentation](https://oskarhs.github.io/AutoHist.jl/dev/methods/).

- Regular Histograms:
    - Regular random histogram, "bayes"
    - L2 cross-validation, "l2cv"
    - Kullback-Leibler cross-validation: "klcv"
    - AIC, "aic"
    - BIC, "bic"
    - Birgé and Rozenholc's criterion, "br"
    - Normalized Maximum Likelihood, "nml"
    - Minimum Description Length, "mdl"
    - Sturges rule, "sturges"
    - Freedman and Diaconis' rule, "fd"
    - Scott's rule, "scott"
    - Wand's rule, "wand"
- Irregular Histograms:
    - Irregular random histogram, "bayes"
    - L2 cross-validation, "l2cv"
    - Kullback-Leibler cross-validation: "klcv"
    - Rozenholc et al. penalty R: "penR"
    - Rozenholc et al. penalty B: "penB"
    - Normalized Maximum Likelihood: "nml"

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