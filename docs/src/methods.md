# Supported Methods

Before we describe the methods included here in more detail, we introduce some notation. For ease of exposition, we present all methods covered here in the context of estimating the density of a sample ``x_1, x_2, \ldots, x_n`` on the unit interval, but note that extending the procedures presented here to other compact intervals. For some further background on histograms, we reccomend the excellent review by Birgé and Rozenholc (2006).

We let ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` denote a partition of ``[0,1]`` into ``k`` intervals and write ``|\mathcal{C}|`` for the length of an interval ``\mathcal{C}``. We can then write a histogram density estimate by

```math
\widehat{f}(x) = \sum_{j=1}^k \frac{\widehat{\theta}_j}{|\mathcal{I}_j|}\mathbf{1}_{\mathcal{I}_j}(x), \quad x\in [0,1],
```
where ``\mathbf{1}_{\mathcal{I}_j}`` is the indicator function, ``\widehat{\theta}_j \geq 0`` for all ``j`` and ``\sum_{j=1}^k \widehat{\theta}_j = 1``.

For most of the methods considered here, the estimated bin probabilities are the maximum likelihood estimates ``\widehat{\theta}_j = N_j/n``, where ``N_j = \sum_{i=1}^n \mathbb{1}_{\mathcal{I}_j}(x_i)`` is number of observations landing in interval ``\mathcal{I}_j`` . The exception to this rule is the Bayesian approach of Simensen et al. (2025), which uses the Bayes estimator ``\widehat{\theta}_j = (5/k + N_j)/(5+n)`` instead.

The goal of an automatic histogram procedure is to find a partition ``\mathcal{I}`` based on the sample alone which produces a reasonable density estimate. Regular histogram procedures only consider regular partitions, where all intervals in the partition are of equal length, so that one only needs to determine the number ``k`` of bins. Irregular histograms allow for partitions with intervals of unequal length, and try to determine both the number of bins and the locations of the cutpoints between the intervals. In all the irregular procedures covered here, we attempt to find best partition according to a criterion among all partitions with endpoints belonging to a given discrete mesh.

## Regular histograms
The following presents the selection criteria maximized by each keyword supported by the `histogram_regular` function. In each case the chosen value of ``k`` is understood to be the maximizer of the expression in question, and ``\mathcal{I} = (\mathcal{I}_1, \mathcal{I}_2, \ldots, \mathcal{I}_k)`` is the corresponding partition of ``[0,1]`` consisting of ``k`` equal-length bins.

#### bayes:
Consists of maximizing the log-marginal likelihood for given ``k``,
```math
   \quad n\log (k) + \sum_{j=1}^k \big\{\log \Gamma(a_j + N_j) - \log \Gamma(a_j)\big\} + \log p_n(k).
```
Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins.

The particular choices ``a_j = 0.5`` and ``p_n(k)\propto 1`` were suggested by Knuth (2019).

#### aic:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (n) + \sum_{j=1}^k N_j \log (N_j/n) - k.
```
The aic criterion was proposed by Taylor (1987) for histograms.
#### bic:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (n) + \sum_{j=1}^k N_j \log (N_j/n) - \frac{k}{2}\log(n).
```
#### br:
Consists of maximizing a penalized log-likelihood,
```math
    n\log (n) + \sum_{j=1}^k N_j \log (N_j/n) - k - \log^{2.5}(k).
```
This criterion was proposed by Birgé and Rozenholc (2006).
#### l2cv:
Consists of maximizing a L2 leave-one-out cross-validation criterion,
```math
    -2k + k\frac{n+1}{n^2}\sum_{j=1}^k N_j^2.
```
This approach to histogram density estimation was first considered by Rudemo (1982).
#### klcv:
Consists of maximizing a Kullback-Leibler leave-one-out cross-validation criterion,
```math
    n\log(k) + \sum_{j=1}^k N_j\log (N_j-1).
```
This approach was first studied by Hall (1990).
#### mdl:
Consists of finding the model providing the shortest encoding of the data, which is equivalent to maximization of
```math
    n\log(k) + \sum_{j=1}^k \big(N_j-\frac{1}{2}\big)\log\big(N_j-\frac{1}{2}\big) - \big(n-\frac{k}{2}\big)\log\big(n-\frac{k}{2}\big) - \frac{k}{2}\log(n).
```
The minimum description length principle was first applied to histogram estimation by Hall and Hannan (1988).

## References
Simensen, O. H., Christensen, D. & Hjort, N. L. (2025). Random Irregular Histograms. _arXiv preprint_. doi: [10.48550/ARXIV.2505.22034](https://doi.org/10.48550/ARXIV.2505.22034)

Taylor, C. C. (1987). Akaike’s information criterion and the histogram. _Biometrika_. 74, 636–639.
doi: [10.1093/biomet/74.3.636](https://doi.org/10.1093/biomet/74.3.636)

Rozenholc, Y., Mildenberger, T., & Gather, U. (2010). Combining regular and irregular histograms by penalized likelihood. _Computational Statistics & Data Analysis_. 54, 3313–3323. doi: [10.1016/j.csda.2010.04.021](https://doi.org/10.1016/j.csda.2010.04.021)

Birgé, L., & Rozenholc, Y. (2006). How many bins should be put in a regular histogram. _ESAIM: Probability and Statistics_. 10, 24–45. doi: [10.1051/ps:2006001](https://doi.org/10.1051/ps:2006001)

Rudemo, M. (1982). Empirical choice of histograms and kernel density estimators. _Scandinavian Journal of Statistics_. 9, 65-78

Hall, P. (1990). Akaike’s information criterion and Kullback–Leibler loss for histogram density estimation.
_Probability Theory and Related Fields_. 85, 449–467. doi: [10.1007/BF01203164](https://doi.org/10.1007/BF01203164)

Hall, P. and Hannan, E. J. (1988). On stochastic complexity and nonparametric density estimation. _Biometrika_.
75, 705–714. doi: [10.1093/biomet/75.4.705](https://doi.org/10.1093/biomet/75.4.705)

Knuth, K. H. (2019). Optimal data-based binning for histograms and histogram-based probability density
models. _Digital Signal Processing_ 95. doi: [10.1016/j.dsp.2019.102581](https://doi.org/10.1016/j.dsp.2019.102581)