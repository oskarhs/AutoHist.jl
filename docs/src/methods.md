# Supported Methods

Before we describe the methods included here in more detail, we introduce some notation. For ease of exposition, we present all methods covered here in the context of estimating the density of a sample ``x_1, x_2, \ldots, x_n`` on the unit interval, but note that extending the procedures presented here to other compact intervals. For some further background on histograms, we reccomend the excellent review by Birgé and Rozenholc (2006).

We let ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` denote a partition of ``[0,1]`` into ``k`` intervals and write ``|\mathcal{C}|`` for the length of an interval ``\mathcal{C}``. We can then write a histogram density estimate by

```math
\widehat{f}(x) = \sum_{j=1}^k \frac{\widehat{\theta}_j}{|\mathcal{I}_j|}\mathbf{1}_{\mathcal{I}_j}(x), \quad x\in [0,1],
```
where ``\mathbf{1}_{\mathcal{I}_j}`` is the indicator function, ``\widehat{\theta}_j \geq 0`` for all ``j`` and ``\sum_{j=1}^k \widehat{\theta}_j = 1``.

For most of the methods considered here, the estimated bin probabilities are the maximum likelihood estimates ``\widehat{\theta}_j = N_j/n``, where ``N_j = \sum_{i=1}^n \mathbb{1}_{\mathcal{I}_j}(x_i)``. The exception to this rule is the Bayesian approach of Simensen et al. (2025), which uses the Bayes estimator ``\widehat{\theta}_j = (5/k + N_j)/(5+n)`` instead.

The goal of an automatic histogram procedure is to find a partition ``\mathcal{I}`` based on the sample alone which produces a reasonable density estimate. Regular histogram procedures only consider regular partitions, where all intervals in the partition are of equal length, so that one only needs to determine the number ``k`` of bins. Irregular histograms allow for partitions with intervals of unequal length, and try to determine both the number of bins and the locations of the cutpoints between the intervals.

In the sequel, we let ``\log L_n(k) = \sum_{j=1}^k N_j \log (N_j/[n |\mathcal{I}_j|])`` be the maximized log-likelihood.

## Regular histograms
The following presents the selection criteria maximized by each keyword supported by the `histogram_regular` function. In each case the chosen value of ``k`` is understood to be the maximizer of the expression in question, and ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` is the corresponding partition of ``[0,1]`` consisting of ``k`` equal-length bins .
- bayes [Simensen et al. (2025)]:

```math
    n\log (k) + \sum_{j=1}^k \big\{\log \Gamma(a_j + N_j) - \log \Gamma(a_j + N_j)\big\} + \log p_n(k)
```
- 
    Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins.
- aic:
    Maximize ``\log L_n(k) - k``
- bic:
    Maximize ``\log L_n(k) - 0.5k\log(n)``
- br [Birgé and Rozenholc (2006)]:
    Maximize ``\log L_n(k) - k - \log^{2.5}(k)``
- l2cv [Rudemo (1982)]:
    Maximize ````



## References
Simensen, O. H., Christensen, D. & Hjort, N. L. (2025). Random Irregular Histograms. _arXiv preprint_. doi: [10.48550/ARXIV.2505.22034](https://doi.org/10.48550/ARXIV.2505.22034)

Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. _Computational Statistics & Data Analysis_, 54, 3313–3323. doi: [10.1016/j.csda.2010.04.021](https://doi.org/10.1016/j.csda.2010.04.021)

Birgé, L., & Rozenholc, Y. (2006). How many bins should be put in a regular histogram. _ESAIM: Probability and Statistics_, 10, 24–45. doi: [10.1051/ps:2006001](https://doi.org/10.1051/ps:2006001)

Rudemo, M. (1982). Empirical choice of histograms and kernel density estimators. _Scandinavian Journal of Statistics_, 9, 65-78