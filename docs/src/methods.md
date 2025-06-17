# Supported Methods
This page provides background on each histogram method supported through the `rule` argument. Our presentation is intended to be rather brief, and we do as such not cover the theoretical underpinnings of each method in great detail. For some further background on automatic histogram procedures and the theory behind them, we recommend the excellent reviews contained in the articles of [Birgé and Rozenholc (2006)]((https://doi.org/10.1016/j.csda.2010.04.021)) and [Davies et al. (2009)](https://doi.org/10.1051/ps:2008005).

For ease of exposition, we present all methods covered here in the context of estimating the density of a sample ``\boldsymbol{x} = (x_1, x_2, \ldots, x_n)`` on the unit interval, but note that extending the procedures presented here to other compact intervals is possible through a suitable affine transformation. In particular, if a density estimate with support ``[a,b]`` is desired, we can scale the data to the unit interval through ``z_i = (x_i - a)/(b-a)``, and apply the methods on this transformed sample and rescale the resulting density estimate to ``[a,b]``. In cases where the support of the density is unknown, a natural choice is ``a = x_{(1)}`` and ``b = x_{(n)}``. Cases where only the lower or upper bound is known can be handled similarly. The transformation used to construct the histogram can be controlled through the `support` keyword, where the default argument `support=(-Inf, Inf)` uses the order statistics-based approach described above.

Before we describe the methods included here in more detail, we introduce some notation. We let ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` denote a partition of ``[0,1]`` into ``k`` intervals and write ``|\mathcal{I}_j|`` for the length of interval ``\mathcal{I}_j``. The intervals in the partition ``\mathcal{I}`` can be either right- or left-closed. Whether a left- or right-closed partition is used to draw the histogram is controlled by the keyword argument `closed`, with options `:left` and `:right` (default). This choice is somewhat arbitrary, but is unlikely to matter much in practical applications.

Based on a partition ``\mathcal{I}``, we can write down the corresponding histogram density estimate by

```math
\widehat{f}(x) = \sum_{j=1}^k \frac{\widehat{\theta}_j}{|\mathcal{I}_j|}\mathbf{1}_{\mathcal{I}_j}(x), \quad x\in [0,1],
```
where ``\mathbf{1}_{\mathcal{I}_j}`` is the indicator function, ``\widehat{\theta}_j \geq 0`` for all ``j`` and ``\sum_{j=1}^k \widehat{\theta}_j = 1``. 

For most of the methods considered here, the estimated bin probabilities are the maximum likelihood estimates ``\widehat{\theta}_j = N_j/n``, where ``N_j = \sum_{i=1}^n \mathbb{1}_{\mathcal{I}_j}(x_i)`` is number of observations landing in interval ``\mathcal{I}_j`` . The exception to this rule is are the two Bayesian approaches, which uses the Bayes estimator ``\widehat{\theta}_j = (a_j + N_j)/(a+n)`` for ``(a_1, \ldots, a_k) \in (0,\infty)^k`` and ``a = \sum_{j=1}^k a_j`` instead.

The goal of an automatic histogram procedure is to find a partition ``\mathcal{I}`` based on the sample alone which produces a reasonable density estimate. Regular histogram procedures only consider regular partitions, where all intervals in the partition are of equal length, so that one only needs to determine the number ``k`` of bins. Irregular histograms allow for partitions with intervals of unequal length, and try to determine both the number of bins and the locations of the cutpoints between the intervals. In all the irregular procedures covered here, we attempt to find best partition according to a criterion among all partitions with endpoints belonging to a given discrete mesh.

## Irregular histograms
The following section describes how each value of the `rule` keyword supported by the `histogram_irregular` function selects the optimal histogram partition. In each case, the best partition is selected among the subset of interval partitions of the unit interval that have cut points belonging to a discrete set of cardinality ``k_n-1``. The construction of the candidate cut point set can be controlled through the `grid` keyword argument, with options `:regular` (default), `:data` and `:quantile`.

#### bayes:
Consists of maximizing the log-marginal likelihood conditional on the partition ``\mathcal{I} = (\mathcal{I}_1, \ldots, \mathcal{I}_k)``,
```math
    \sum_{j=1}^k \big\{\log \Gamma(a_j + N_j) - \log \Gamma(a_j) - N_j\log|\mathcal{I}_j|\big\} + \log p_n(k) - \log \binom{k_n-1}{k-1}
```
Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins, which can be controlled by supplying a function to the `logprior` keyword argument. The default value is ``p_n(k) \propto 1``. Here, ``a_j = a/k``, for a scalar ``a > 0`` which can be controlled by the user through the keyword argument `a` (default `a=5.0`).

This approach to irregular histograms was pioneered by [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).


#### pena:
Consists of maximizing a penalized log-likelihood,
```math
    \sum_{j=1}^k N_j \log (N_j/|\mathcal{I}_j|) - \log \binom{k_n-1}{k-1} - k - 2\log(k) - \sqrt{2(k-1)\Big[\log \binom{k_n-1}{k-1}+ 2\log(k)\Big]}*.
```
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).

#### penb:
Consists of maximizing a penalized log-likelihood,
```math
    \sum_{j=1}^k N_j \log (N_j/|\mathcal{I}_j|) - \log \binom{k_n-1}{k-1} - k - \log^{2.5}(k).
```
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
#### penr:
Consists of maximizing a penalized log-likelihood,
```math
    \sum_{j=1}^k N_j \log (N_j/|\mathcal{I}_j|) - \frac{1}{2n}\sum_{j=1}^k \frac{N_j}{|\mathcal{I}_j|} - \log \binom{k_n-1}{k-1} - \log^{2.5}(k).
```
This criterion was also suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
#### l2cv:
Consists of maximization of an L2 leave-one-out cross-validation criterion,
```math
    \frac{n+1}{n}\sum_{j=1}^k \frac{N_j^2}{|\mathcal{I}_j|} - 2\sum_{j=1}^k \frac{N_j}{|\mathcal{I}_j|}.
```
This approach dates back to [Rudemo (1982)](https://www.jstor.org/stable/4615859).
#### klcv:
Consists of maximization of a Kullback-Leibler cross-validation criterion,
```math
    \sum_{j=1}^k N_j\log(N_j-1) - \sum_{j=1}^k N_j\log |I_j|,
```
where the maximmization is over all partitions with ``N_j \geq 2`` for all ``j``.
This approach was, to our knowledge, first pursued by Simensen et al. (2025).
#### nml:
Consists of maximization of a penalized likelihood,
```math
\begin{aligned}
    &\sum_{j=1}^k N_j\log \frac{N_j}{|\mathcal{I}_j|} - \frac{k-1}{2}\log(n/2) - \log\frac{\sqrt{\pi}}{\Gamma(k/2)} - n^{-1/2}\frac{\sqrt{2}k\Gamma(k/2)}{3\Gamma(k/2-1/2)} \\
    &- n^{-1}\left(\frac{3+k(k-2)(2k+1)}{36} - \frac{\Gamma(k/2)^2 k^2}{9\Gamma(k/2-1/2)^2} \right)  - \log \binom{k_n-1}{k-1}
\end{aligned}
```
This a variant of this criterion first suggested by [Kontkanen and Myllymäki (2007)](https://proceedings.mlr.press/v2/kontkanen07a.html). The above criterion uses an asymptotic expansion of their proposed penalty term, as their proposed penalty can be quite expensive to evaluate.

## Regular histograms
The following section details how each value of the `rule` keyword supported by the `histogram_regular` function selects the number ``k`` of bins to draw a histogram automatically based on a random sample. In the following, ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` is the corresponding partition of ``[0,1]`` consisting of ``k`` equal-length bins. In cases where the value of the number of bins is computed by maximizing an expression, we look for the best regular partition among all regular partitions consisting of no more than ``k_n`` bins.

#### bayes:
Consists of maximizing the log-marginal likelihood for given ``k``,
```math
   n\log (k) + \sum_{j=1}^k \big\{\log \Gamma(a_j + N_j) - \log \Gamma(a_j)\big\} + \log p_n(k).
```
Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins, which can be controlled by supplying a function to the `logprior` keyword argument. The default value is ``p_n(k) \propto 1``. Here, ``a_j = a/k``, for a scalar ``a > 0``, possibly depending on ``k``. The value of ``a`` can be set by supplying a fixed, positive scalar or a function ``a(k)`` to the keyword argument `a`.

The particular choices ``a_j = 0.5`` and ``p_n(k)\propto 1`` were suggested by [Knuth (2019)](https://doi.org/10.1016/j.dsp.2019.102581).

#### aic:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (k) + \sum_{j=1}^k N_j \log (N_j/n) - k.
```
The aic criterion was proposed by [Taylor (1987)](https://doi.org/10.1093/biomet/74.3.636) for histograms.
#### bic:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (k) + \sum_{j=1}^k N_j \log (N_j/n) - \frac{k}{2}\log(n).
```
#### br:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (k) + \sum_{j=1}^k N_j \log (N_j/n) - k - \log^{2.5}(k).
```
This criterion was proposed by [Birgé and Rozenholc (2006)](https://doi.org/10.1051/ps:2006001).
#### l2cv:
Consists of maximizing a L2 leave-one-out cross-validation criterion,
```math
    -2k + k\frac{n+1}{n^2}\sum_{j=1}^k N_j^2.
```
This approach to histogram density estimation was first considered by [Rudemo (1982)](https://www.jstor.org/stable/4615859
).
#### klcv:
Consists of maximizing a Kullback-Leibler leave-one-out cross-validation criterion,
```math
    n\log(k) + \sum_{j=1}^k N_j\log (N_j-1),
```
where the maximmization is over all regular partitions with ``N_j \geq 2`` for all ``j``.
This approach was first studied by [Hall (1990)]((https://doi.org/10.1007/BF01203164)).
#### mdl:
Consists of finding the model providing the shortest encoding of the data, which is equivalent to maximization of
```math
    n\log(k) + \sum_{j=1}^k \big(N_j-\frac{1}{2}\big)\log\big(N_j-\frac{1}{2}\big) - \big(n-\frac{k}{2}\big)\log\big(n-\frac{k}{2}\big) - \frac{k}{2}\log(n),
```
where the maximmization is over all regular partitions with ``N_j \geq 1`` for all ``j``.
The minimum description length principle was first applied to histogram estimation by [Hall and Hannan (1988)](https://doi.org/10.1093/biomet/75.4.705).

#### nml:
Consists of maximization of a penalized likelihood,
```math
\begin{aligned}
    &\sum_{j=1}^k N_j\log \frac{N_j}{|\mathcal{I}_j|} - \frac{k-1}{2}\log(n/2) - \log\frac{\sqrt{\pi}}{\Gamma(k/2)} - n^{-1/2}\frac{\sqrt{2}k\Gamma(k/2)}{3\Gamma(k/2-1/2)} \\
    &- n^{-1}\left(\frac{3+k(k-2)(2k+1)}{36} - \frac{\Gamma(k/2)^2 k^2}{9\Gamma(k/2-1/2)^2} \right).
\end{aligned}
```
This is a regular variant of the normalized maximum likelihood criterion considered by [Kontkanen and Myllymäki (2007)](https://proceedings.mlr.press/v2/kontkanen07a.html).

#### Sturges' rule:
The number ``k`` of bins is computed according to the formula
```math
    k = \lceil \log_2(n) \rceil + 1.
```
This classical rule, due to [Sturges (1926)](https://doi.org/10.1080/01621459.1926.10502161), is the default for determining the number of bins in R.

#### Freedman and Diaconis' rule:
The number ``k`` of bins is computed according to the formula
```math
    k = \big\lceil\frac{n^{1/3}}{2\mathrm{IQR}(\boldsymbol{x})}\big\rceil,
```
where ``\mathrm{IQR}(\boldsymbol{x})`` is the sample interquartile range. This rule dates back to [Freedman and Diaconis (1982)](https://doi.org/10.1007/BF01025868) and is the default bin selection rule used by the `histogram()` function from Plots.jl.

#### Scott's rule:
The number ``k`` of bins is computed according to the formula
```math
    k = \big\lceil \hat{\sigma}^{-1}(24\sqrt{\pi})^{-1/3}n^{1/3}\big\rceil,
```
where ``\hat{\sigma}`` is the sample standard deviation. Scott's normal reference rule was first proposed by [Scott (1979)](https://doi.org/10.1093/biomet/66.3.605).

#### Wand's rule
A more sophisticated version of Scott's rule, Wand's rule proceeds by determining the bin width ``h`` as
```math
    h = \Big(\frac{6}{\hat{C}(f_0) n}\Big)^{1/3},
```
where ``\hat{C}(f_0)`` is an estimate of the functional ``C(f_0) = \\int \\big\\{f_0'(x)\\big\\}^2\\mspace{2mu}\\mathrm{d}x``. The corresponding number of bins ``k = \lceil h^{-1}\rceil``. The full details on this method are given in [Wand (1997)](https://doi.org/10.2307/2684697).
The density estimate is computed based on a scale estimate, which can be controlled through the `scale` keyword argument. Possible choices are `:stdev`, `:iqr` which uses an estimate based on the sample standard deviation or the sample interquartile range as a scale estimate. The default choice `:minim` uses the minimum of the above estimates.
The `level` keyword controls the number of stages of functional estimation used to compute ``\hat{C}``, and can take values `0, 1, 2, 3, 4, 5`, with the default value being `level=2`. The choice `level=0` corresponds to Scott's rule under the chosen scale estimate.


## References
Simensen, O. H., Christensen, D. & Hjort, N. L. (2025). Random Irregular Histograms. _arXiv preprint_. doi: [10.48550/ARXIV.2505.22034](https://doi.org/10.48550/ARXIV.2505.22034)

Taylor, C. C. (1987). Akaike’s information criterion and the histogram. _Biometrika_, **74**, 636–639.
doi: [10.1093/biomet/74.3.636](https://doi.org/10.1093/biomet/74.3.636)

Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. _Computational Statistics & Data Analysis_, **54**, 3313–3323. doi: [10.1016/j.csda.2010.04.021](https://doi.org/10.1016/j.csda.2010.04.021)

Birgé, L., & Rozenholc, Y. (2006). How many bins should be put in a regular histogram. _ESAIM: Probability and Statistics_, **10**, 24–45. doi: [10.1051/ps:2006001](https://doi.org/10.1051/ps:2006001)

Rudemo, M. (1982). Empirical choice of histograms and kernel density estimators. _Scandinavian Journal of Statistics_, **9**, 65-78

Hall, P. (1990). Akaike’s information criterion and Kullback–Leibler loss for histogram density estimation.
_Probability Theory and Related Fields_, **85**, 449–467. doi: [10.1007/BF01203164](https://doi.org/10.1007/BF01203164)

Hall, P. and Hannan, E. J. (1988). On stochastic complexity and nonparametric density estimation.
_Biometrika_, **75**, 705–714. doi: [10.1093/biomet/75.4.705](https://doi.org/10.1093/biomet/75.4.705)

Knuth, K. H. (2019). Optimal data-based binning for histograms and histogram-based probability density models.
_Digital Signal Processing_, **95**, doi: [10.1016/j.dsp.2019.102581](https://doi.org/10.1016/j.dsp.2019.102581)

Kontkanen, P. and Myllymäki, P. (2007). Mdl histogram density estimation. _Proceedings of the Eleventh International Conference on Artificial Intelligence and Statistics_, **2**, 219–226

Sturges, H. A. (1926). The choice of a class interval. _Journal of the American Statistical Association_, **21**, 65–66. doi: [10.1080/01621459.1926.10502161](https://doi.org/10.1080/01621459.1926.10502161).

Freedman, D. and Diaconis, P. (1981) On the histogram as a density estimator: L2 theory.
_Zeitschrift für Wahrscheinlichkeitstheorie und verwandte Gebiete_, **57**, 453–476.
doi: [10.1007/BF01025868](https://doi.org/10.1007/BF01025868).

Scott, D. W. (1979). On optimal and data-based histograms. _Biometrika_, **66**, 605–610,
doi: [10.1093/biomet/66.3.605](https://doi.org/10.1093/biomet/66.3.605).

Wand, M. P. (1997). Data-based choice of histogram bin width. _The American Statistician_, **51**, 59–64.
doi: [10.2307/2684697](https://doi.org/10.2307/2684697)

Davies, P. L., Gather, U., Nordman, D., and Weinert, H. (2009). A comparison of automatic histogram constructions. _ESAIM: Probability and Statistics_, **13**, 181–196. doi: [10.1051/ps:2008005](https://doi.org/10.1051/ps:2008005).