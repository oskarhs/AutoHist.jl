# Algorithms for constructing irregular histograms

Constructing data-adaptive irregular histograms is in general a difficult problem from a computational perspective, and as a result computing the exact optimal partition is often impractical for larger sample sizes. This package solves the problem of irregular histogram construction via heuristics that combine a greedy search procedure with dynamic programming techniques to quickly compute a nearly optimal partition. Examples showcasing the use of the provided algorithms in toy problems can be found [here](examples/algorithm_choice.md). Note that the default option for the `alg` keyword offers a reasonable tradeoff between accuracy and computational efficiency, and in cases where one simply wants to draw an irregular histogram quickly for a given dataset, these default choices will typically yield a reasonable density estimate within a reasonable amount of time. However, if better performance or additional accuracy is desired, then a more fine-tuned approach to selecting the algorithm used and its hyperparameters is needed.

## Problem description
All the irregular histogram methods supported by this library are the product of solving an optimization problem of the form
```math
    \max_{\boldsymbol{\tau}\in \mathcal{T}}\big\{\sum_{j=1}^k \Phi(\tau_{j-1}, \tau_{j}] + \Psi(k)\big\},
```
where ``\mathcal{T}`` is a set of candidate cutpoints. While it is generally desirable to use a moderate number of possible cutpoints relative to the number of samples, this comes at a heavy computational cost if an exact solution of the above optimization is desired as the runtime complexity of the algorithm is cubic in the number of candidates. We note that for the special case where ``\Psi(k) = 0`` for all ``k``, a more efficient algorithm is available, allowing for a quadratic-time solution. In both cases, computing the exact solution becomes unfeasible for larger sample sizes.

To ease the computational burden we adopt a greedy search heuristic to construct a subset of possible cutpoints when the number of candidates is large, and then run the optimization algorithm of choice on this smaller set. To find the optimal partition for a given set of candidate cutpoints, this package includes two solvers; an [exact dynamic programming algorithm](#Dynamic-programming) and a heuristic [greedy pruned dynamic programming algorithm](#Greedy-pruned-dynamic-programming). The algorithm used can be controlled via the `alg` keyword passed to [`fit`](@ref) or [`histogram_irregular`](@ref).

## Dynamic programming
The choice `alg = DP()` results in the use of the exact dynamic programming algorithm of [Kanazawa (1988)](https://doi.org/10.1080/03610928808829688) or, if applicable, the exact optimal partitioning algorithm of [Jackson et al. (2005)](https://doi.org/10.1109/LSP.2001.838216) to find the optimal histogram partition. Runtime complexity is cubic in the number of candidate cutpoints for the former and quadratic for the latter.
```@docs
DP
```

## Greedy pruned dynamic programming
The choice `alg = GPDP()` results in the use of the heuristic greedy dynamic programming algorithm of [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034) for finding an approximately optimal histogram partition. Runtime complexity is quadratic in the number of candidate cutpoints.
```@docs
GPDP
```