---
title: 'AutoHist.jl: A Julia package for fast and automatic histogram construction'
tags:
  - Julia
  - statistics
  - histograms
  - density estimation
authors:
  - name: Oskar Høgberg Simensen
    orcid: 0009-0009-5056-7324
    affiliation: 1
affiliations:
 - name: Karlsruhe Institute of Technology
   index: 1
output:
  github_document:
    html_preview: true
bibliography: paper.bib
date: 16 July 2025
---

# Summary

`AutoHist.jl` is a `Julia` [@bezanson2017julia] package for fitting histograms to univariate data, with an automatic and data-based selection of the histogram partition.
It currently supports 7 irregular and 12 regular automatic histogram procedures from the statistical literature. `AutoHist.jl` provides extensions for `Plots.jl` [@christ2023plots] and `Makie.jl` [@danisch2021makie],
allowing for simple and intuitive visualization of the resulting histograms.

# Statement of need

To this day, histograms remain one of the most widely used nonparametric density estimators. Their popularity is in no doubt due to their simplicity and interpretability,
and they are routinely used by practitioners with limited mathematical backgrounds. 
Although contructing a histogram for a given interval partition is a simple task, the quality of the resulting density estimate is very sensitive to the choice of bins.
As a result, the task of designing automatic histogram procedures, where the number of bins and their location are chosen automatically based on the sample has recieved considerable interest in the statistics community.
Despite advances in our understanding of different bin selection rules, many popular software libraries only include a limited number of simple rules for selecting a regular histogram partition,
and typically do not provide any support for automatic irregular histogram construction.
`AutoHist.jl` fills this gap by providing a fast implementation of state-of-the-art regular and irregular bin selection algorithms from the statistics literature.

We note that some automatic histogram selection rules have been implemented in `Julia`, typically as part of plotting libraries.
The `Plots.jl` package provides an implementation of some plug-in rules for regular histograms, while the `StatsPlots` package offers an implementation of equal-area histograms.
The `R` package `histogram` [@mildenberger2019histogram] supports some regular and irregular histogram methods, but their implementation covers fewer criteria.
Our implementation has the additional advantage of avoiding the awkward conversion of histogram objects from `R` to an appropriate type in `Julia`.

# Examples

To illustrate the basic use of the software, we fit an irregular histogram based on the Bayesian criterion and a regular histogram bas on AIC to a standard normal random sample of size $n = 10^6$.

```julia
using AutoHist, Random, Distributions

n = 10^6
x = rand(Xoshiro(1812), Normal(), n)            # synthetic data

h_irr = fit(AutomaticHistogram, x; rule=:bayes) # fit an irregular histogram
h_reg = fit(AutomaticHistogram, x; rule=:aic)   # fit a regular histogram
```

The call to the `fit` method returns an object of type `AutomaticHistogram`, with fields recording the chosen histogram partition, estimated density and bin counts. 

```julia
using Makie, CairoMakie
```


# References