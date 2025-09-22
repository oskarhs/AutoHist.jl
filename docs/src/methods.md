# Supported Methods
This page provides background on each histogram method supported through the `rule` argument. Our presentation is intended to be rather brief, and we therefore do not cover the theoretical underpinnings of each method in great detail. For some further background on automatic histogram procedures and the theory behind them, we recommend the excellent reviews contained in the articles of [Birgé and Rozenholc (2006)](https://doi.org/10.1016/j.csda.2010.04.021) and [Davies et al. (2009)](https://doi.org/10.1051/ps:2008005).

For ease of exposition, we present all methods covered here in the context of estimating the density of a sample ``\boldsymbol{x} = (x_1, x_2, \ldots, x_n)`` on the unit interval, but note that extending the procedures presented here to other compact intervals is possible through a suitable affine transformation. In particular, if a density estimate with support ``[a,b]`` is desired, we can scale the data to the unit interval through ``z_i = (x_i - a)/(b-a)``, and apply the methods on this transformed sample and rescale the resulting density estimate to ``[a,b]``. In cases where the support of the density is unknown, a natural choice is ``a = x_{(1)}`` and ``b = x_{(n)}``. Cases where only the lower or upper bound is known can be handled similarly. The transformation used to construct the histogram can be controlled through the `support` keyword, where the default argument `support=(-Inf, Inf)` uses the order statistics-based approach described above.

#### Notation
Before we describe the methods included here in more detail, we introduce some notation. We let ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` denote a partition of ``[0,1]`` into ``k`` intervals and write ``|\mathcal{I}_j|`` for the length of interval ``\mathcal{I}_j``. The intervals in the partition ``\mathcal{I}`` can be either right- or left-closed. Whether a left- or right-closed partition is used to draw the histogram is controlled by the keyword argument `closed`, with options `:left` and `:right` (default). This choice is somewhat arbitrary, but is unlikely to matter much in practical applications.

Based on a partition ``\mathcal{I}``, we can write down the corresponding histogram density estimate by

```math
\widehat{f}(x) = \sum_{j=1}^k \frac{\widehat{\theta}_j}{|\mathcal{I}_j|}\mathbf{1}_{\mathcal{I}_j}(x), \quad x\in [0,1],
```
where ``\mathbf{1}_{\mathcal{I}_j}`` is the indicator function, ``\widehat{\theta}_j \geq 0`` for all ``j`` and ``\sum_{j=1}^k \widehat{\theta}_j = 1``. 

For most of the methods considered here, the estimated bin probabilities are the maximum likelihood estimates ``\widehat{\theta}_j = N_j/n``, where ``N_j = \sum_{i=1}^n \mathbb{1}_{\mathcal{I}_j}(x_i)`` is the number of observations landing in interval ``\mathcal{I}_j`` . The exception to this rule are the two Bayesian approaches [`RIH`](@ref) and [`RRH`](@ref), which uses the Bayes estimator ``\widehat{\theta}_j = (a_j + N_j)/(a+n)`` for ``(a_1, \ldots, a_k) \in (0,\infty)^k`` and ``a = \sum_{j=1}^k a_j`` instead.

The goal of an automatic histogram procedure is to find a partition ``\mathcal{I}`` based on the sample alone which produces a reasonable density estimate. Regular histogram procedures only consider regular partitions, where all intervals in the partition are of equal length, so that one only needs to determine the number ``k`` of bins. Irregular histograms allow for partitions with intervals of unequal length, and try to determine both the number of bins and the locations of the cutpoints between the intervals.

#### A short note on using different rules
In order to fit a histogram using a specific `rule`, we call `fit(AutomaticHistogram, x, rule)`, where `x` is the data vector. For many of the rules discussed below, the user can specify additional rule-specific keywords to `rule`, providing additional control over the supplied method when desired. We also provide a set of default values for these parameters, so that the user may for instance call `fit(AutomaticHistogram, x, AIC())` to fit a regular histogram using the AIC criterion without having to worry about explicitly passing any keyword arguments.

## Irregular histograms
The following section provides a description of all the irregular histogram rules that have been implemented in AutoHist.jl. In each case, the best partition is selected among the subset of interval partitions of the unit interval that have cut points belonging to a discrete set of cardinality ``k_n+1``[^1]. In all the irregular procedures covered here, we attempt to find best partition according to a goodness-of-fit criterion among all partitions with endpoints belonging to a given discrete mesh ``\{\tau_{j}\colon 0\leq j \leq k_n\}``.

[^1]: Note that the endpoints ``0`` and ``1`` are included in every candidate partition.

#### Random irregular histogram
```@docs
RIH
```

#### Rozenholc, Mildenberger & Gather penalty A
```@docs
RMG_penA
```

#### Rozenholc, Mildenberger & Gather penalty B
```@docs
RMG_penB
```

#### Rozenholc, Mildenberger & Gather penalty R
```@docs
RMG_penR
```

#### Irregular ``L_2`` leave-one-out cross-validation (L2CV_I)
```@docs
L2CV_I
```

#### Irregular Kullback-Leibler leave-one-out cross-validation (KLCV_I)
```@docs
KLCV_I
```

#### Normalized maximum likelihood, regular (NML_R)
```@docs
NML_I
```

#### Bayesian Blocks (BayesBlocks)
```@docs
BayesBlocks
```

## Regular histograms
The following section details how each value of the `rule` argument selects the number ``k`` of bins to draw a regular histogram automatically based on a random sample. In the following, ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` is the corresponding partition of ``[0,1]`` consisting of ``k`` equal-length bins. In cases where the value of the number of bins is computed by maximizing an expression, we look for the best regular partition among all regular partitions consisting of no more than ``k_n`` bins. For rules falling under this umbrella, ``k_n`` can be controlled through the `maxbins` keyword, as detailed below.

#### Random regular histogram (RRH), Knuth

```@docs
RRH
```

#### AIC

```@docs
AIC
```

#### BIC
```@docs
BIC
```

#### Birgé-Rozenholc (BR)
```@docs
BR
```

#### Regular ``L_2`` leave-one-out cross-validation (L2CV_R)
```@docs
L2CV_R
```

#### Regular Kullback-Leibler leave-one-out cross-validation (KLCV_R)
```@docs
KLCV_R
```

#### Minimum description length (MDL)
```@docs
MDL
```

#### Normalized maximum likelihood, regular (NML_R)
```@docs
NML_R
```

#### Sturges' rule
```@docs
Sturges
```

#### Freedman and Diaconis' rule
```@docs
FD
```

#### Scott's rule
```@docs
Scott
```

#### Wand's rule
```@docs
Wand
```

## Index

```@index
Pages = ["methods.md"]
```