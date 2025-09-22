# Algorithms for constructing irregular histograms

Constructing data-adaptive irregular histograms is in general a difficult problem from a computational perspective, and as a result computing the exact optimal partition is often impractical for larger sample sizes. This package solves the problem of irregular histogram construction via heuristics that combine a greedy search procedure with dynamic programming techniques to quickly compute a nearly optimal partition. Examples showcasing the use of the provided algorithms in toy problems can be found [here](examples/algorithm_choice.md). Note that the default option for the `alg` keyword offers a reasonable tradeoff between accuracy and computational efficiency, and in cases where one simply wants to draw an irregular histogram quickly for a given dataset, these default choices will typically yield a reasonable density estimate within a reasonable amount of time. However, if better performance or additional accuracy is desired, then a more fine-tuned approach to selecting the algorithm used and its hyperparameters is needed.

## Problem description
All the irregular histogram methods supported by this library are the product of solving an optimization problem of the form
```math
    \max_{1\leq k\leq k_n}\max_{\boldsymbol{t}_{0:k}} \Big\{\sum_{j=1}^k \Phi(\tau_{n, t_{j-1}}, \tau_{n, t_{j}}] + \Psi(k)\Big\},
```
where the inner maximum is over integer vectors ``\boldsymbol{t}_{0:k}`` satisfying ``0 = t_0 < t_1 < \cdots < t_{k-1} < t_k = k_n`` and ``\{\tau_{n,j}\colon 0\leq j \leq k_n\}`` are the candidate cutpoints between the partition intervals. While it is generally desirable to use a moderate number of possible cutpoints relative to the number of samples, this comes at a heavy computational cost if an exact solution of the above optimization is desired as the runtime complexity of the algorithm is cubic in the number of candidates. We note that for the special case where ``\Psi(k) = \beta k`` for some scalar ``\beta\in \mathbb{R}`` and all ``k\in \mathbb{N}``, a more efficient dynamic programming algorithm is available, allowing for a quadratic-time solution. In both cases, computing the exact solution becomes unfeasible for larger sample sizes.

To ease the computational burden we adopt a greedy search heuristic to construct a subset of possible cutpoints when the number of candidates is large, and then run the optimization algorithm of choice on this smaller set. To find the optimal partition for a given set of candidate cutpoints, this package includes two solvers based on dynamic programming. The algorithm used can be controlled via the `alg` keyword passed to the rule argument of [`fit`](@ref), see the [methods page](methods.md#irregular-histograms).


## Segment neighbourhood
The choice `alg = SegNeig()` results in the use of the exact dynamic programming algorithm of [Kanazawa (1988)](https://doi.org/10.1080/03610928808829688) to find the optimal partition.
Runtime complexity is cubic in the number of candidate cutpoints.
```@docs
SegNeig
```

## Optimal partitioning
The choice `alg = OptPart()` results in the use of the exact optimal partitioning algorithm of [Jackson et al. (2005)](https://doi.org/10.1109/LSP.2001.838216) to find the optimal histogram partition. This algorithm solves a less general optimization problem than the `SegNeig` algorithm. Runtime complexity is quadratic in the number of candidate cutpoints.
```@docs
OptPart
```