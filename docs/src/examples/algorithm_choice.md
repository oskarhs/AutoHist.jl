# Choice of algorithm
In this section, we empirically assess the efficiency of the different algorithms provided for irregular histograms, and show how heuristics can be used to speed up the computations.[^1]

[^1]: **Note:** The benchmarks presented here were performed on a Windows machine with a 12th Gen Intel® Core™ i7-1255U CPU. Results may vary on systems with different hardware configurations.

## Dynamic programming
As a toy problem, we consider standard normal random samples of using a data-based grid. In this case, the number of candidate cutpoints are ``k_n = n+1``, where ``n`` is the sample size. For smaller samples, we can just compute the exact solution using the default dynamic programming algorithm, available as [`DP`](@ref). The code snippet below illustrates how this algorithm can be explicitly specified when calling `fit`:
```@example algorithm; continued = true
using AutoHist, Distributions
n = 500
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    grid = :data,
    alg = DP(greedy=false)
)
```
Benchmarking the above code snippet yields a mean runtime of around ``60\ \text{ms}`` on my machine. Since dynamic programming is quick for samples of this size, the greedy algorithm will only be used if the number of candidate cutpoints exceeds ``500``. Thus, changing `greedy` to `true` in the above code snippet would produce the same histogram, as there are ``499`` possible cutpoints.

#### Speeding up computations via heuristics

The ``\mathcal{O}(k_n^3)`` runtime of dynamic programming means that computing the optimal solution quickly becomes computationally prohibitive, even for moderate samples. As an example, when doubling the number of samples in the above code snippet to ``n = 1000``, the mean runtime increases to ``450\ \text{ms}``, a rougly ``8``-fold increase. To ensure that the code retains good performance even for larger samples, we have implemented a greedy search heuristic which selects a subset of the candidate cutpoints, and the dynamic programming algorithm is subsequently run on this smaller set. Adopting the heuristic approach can improve performance considerably, but comes at the cost of no longer being guaranteed to find the optimal solution. To showcase the computational advantages of the heuristic approach, we run a benchmark on a normal sample of size ``n = 10^6``.
```@example algorithm; continued = true
n = 10^6
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    grid = :data,
    alg = DP(greedy=true) # NB! greedy=true is the default option
)
```
Benchmarking the above, I find a mean runtime of approximately ``720\ \text{ms}`` which is not much slower than the time it took to compute the exact solution for random samples of size ``n = 10^3``.

The number candidate cutpoints constructed by the greedy search heuristic can be controlled through the `gr_maxbins` keyword argument, which equals the number of selected gridpoints plus one. By default, the greedy algorithm will produce a subset consisting of ``\max\{500, n^{1/3}\}+1`` cutpoints by default (including the edges). For `gr_maxbins1 < gr_maxbins2`, the cutpoint subset formed by the greedy algorithm for `gr_maxbins1` is a subset of that selected with `gr_maxbins2` bins. Thus, increasing the number of candidate cutpoints added this grid will never lead to a worse solution of the original optimization problem. If additional precision is desired in the above example, we can increase `gr_maxbins` to ``1000``:
```@example algorithm; continued = true
n = 10^6
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    grid = :data,
    alg = DP(greedy=true, gr_maxbins=10^3)
)
```
The above code snippet averages a runtime of about ``1.1\ \text{s}`` on my machine.

#### A note on cross-validation based criteria
For the [L2CV](../methods.md#l2cv-(irregular)) and [KLCV](../methods.md#klcv-(irregular)) criteria, it becomes possible to compute the exact solution via a quadratic-time dynamic programming algorithm instead of the cubic-time algorithm used for the other problems. In practice, this means that computing the exact solution is often feasible even as the number of candidate cutpoints becomes quite large. For instance, benchmarking the following code snippet results in a mean runtime of ``4.2\ \text{ms}``, roughly ``2`` orders of magnitude below the runtime from computing the exact solution with the cubic time algorithm.
```@example algorithm
n = 10^3
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    rule = :l2cv,
    grid = :data,
    alg = DP(greedy=false)
)
nothing # hide
```
Due to the superior runtime complexity of the exact algorithm for cross-validation criteria, the exact solution is by default computed when the number of total candidate cutpoints is less than ``3001``, and the greedy search heuristic is used thereafter to build a smaller candidate set as previously. By default `gr_maxbins` is set to ``\max\{3000, \sqrt{n}\}``, but this can be replaced with a user-defined value if better performance or additional accuracy is desired.

## Greedy pruned dynamic programming
The other algorithm supported by this library is a heuristic variant of the cubic-time dynamic programming algorithm, available as [`GPDP`](@ref). This algorithm introduces a greedy step within the dynamic programming recursion which reduces the amount of computation and yields an overall ``\mathcal{O}(k_n^2)`` runtime. The faster asymptotic runtime comes at the cost of no longer being guaranteed to compute the optimal solution. As previously, the `GPDP` algorithm can be combined with the same greedy search heuristic that was used for the [exact dynamic programming algorithm](#dynamic-programming) to reduce the computational burden for larger sample sizes.

The main advantage of using the greedy variant of dynamic programming is that the faster asymptotic runtime makes it feasible to consider a larger number of candidate cutpoints in the dynamic programming algorithm. This means that when encountering larger sample sizes, we can allow for the greedy preprocessing algorithm to select a larger number of candidate cutpoints.

Since this algorithm includes some additional operations at each iteration of the dynamic programming recursion, most of the computational benefits are mostly relevant for larger sample sizes, and running the exact dynamic programming algorithm will therefore often be faster unless the number of candidate changepoints is moderately large.

To showcase the use of the `GPDP` algorithm, consider the following code snippet, which fits an `AutmaticHistogram` to a normal random sample using the `GPDP` algorithm without a greedy search.
```@example algorithm
n = 10^3
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    grid = :data,
    alg = GPDP(greedy=false)
)
nothing # hide
```
The above code block clocks in with an average runtime of about ``104\ \text{ms}``. By default, the greedy search heuristic will only be used if the `gr_maxbins` keyword takes a value greater than ``1000``. By default, `gr_maxbins` is set to ``\max\{1000, \sqrt{n}\}``, so the greedy algorithm is only used for sample sizes greater than ``1000``.

The `GPDP` algorithm comes with the additional keyword `max_cand`, which controls the level of greediness of the algorithm. When set to a larger value, the greedy algorithm will consider a larger set of candidate cutpoints in each iteration of the dynamic programming recursion, potentially leading to solutions of higher quality at the cost of increased computation. By default, `max_cand` is set to a value of ``15``. The below code snippet shows how to run the `GPDP` algorithms with `max_cand` set to ``30``.
```@example algorithm
n = 10^3
fit(
    AutomaticHistogram, 
    rand(Normal(), n);
    grid = :data,
    alg = GPDP(greedy=false, max_cand = 30)
)
nothing # hide
```
This code snippet clocks in at about ``150\ \text{ms}`` on average.

<!-- Could do a final example here that compares statistical efficiency vs. runtime for the two approaches -->