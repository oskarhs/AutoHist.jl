# ------------------------------
# Random irregular histogram 
struct RIH <: AbstractIrregularRule
    a::Real
    logprior::Function
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
end

"""
    RIH(;
        a::Real,
        logprior::Function:=k-> 0.0,
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=SegNeig()
    )

The random irregular histogram criterion.

Consists of maximizing the marginal log-posterior of the partition ``\\mathcal{I} = (\\mathcal{I}_1, \\ldots, \\mathcal{I}_k)``,
```math
    \\sum_{j=1}^k \\big\\{\\log \\Gamma(a_j + N_j) - \\log \\Gamma(a_j) - N_j\\log|\\mathcal{I}_j|\\big\\} + \\log p_n(k) - \\log \\binom{k_n-1}{k-1}
```
Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins, which can be controlled by
supplying a function to the `logprior` keyword argument. 
The default value is ``p_n(k) \\propto 1``. Here, ``a_j = a/k``, for a scalar ``a > 0``, not depending on ``k``.

# Keyword arguments
- `a`: Specifies Dirichlet concentration parameter in the Bayesian histogram model. Must be a fixed positive number, and defaults to `a=5.0`.
- `logprior`: Unnormalized logprior distribution on the number ``k`` of bins. Defaults to a uniform prior, e.g. `logprior(k) = 0` for all `k`.
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only `SegNeig()` for this rule. See [`DP`](@ref) for further details.

# References
This approach to irregular histograms first appeared in [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).
"""
RIH(; a::Real=5.0, logprior::Function=k->0.0, grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=SegNeig()) = RIH(a, logprior, grid, maxbins, alg)

function fit_autohist(x::AbstractVector{T}, rule::RIH, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    if rule.a ≤ 0.0
        throw(DomainError("Supplied value of a must be positive."))
    end
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = (N .+ rule.a*p0) ./ ((n + rule.a)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, rule.a)
    return h
end


# ------------------------------
# RMG penA
struct RMG_penA <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
end

"""
    RMG_penA(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k N_j \\log (N_j/|\\mathcal{I}_j|) - \\log \\binom{k_n-1}{k-1} - k - 2\\log(k) - \\sqrt{2(k-1)\\Big[\\log \\binom{k_n-1}{k-1}+ 2\\log(k)\\Big]}.
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only `SegNeig()` is supported for this rule. See [`DP`](@ref) for further details.

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
RMG_penA(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=SegNeig()) = RMG_penA(grid, maxbins, alg)

function fit_autohist(x::AbstractVector{T}, rule::RMG_penA, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ ((n)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# RMG penB
struct RMG_penB <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
end

"""
    RMG_penB(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k N_j \\log (N_j/|\\mathcal{I}_j|) - \\log \\binom{k_n-1}{k-1} - k - \\log^{2.5}(k).
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only `SegNeig()` is supported for this rule. See [`DP`](@ref) for further details.

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
RMG_penB(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=SegNeig()) = RMG_penB(grid, maxbins, alg)

function fit_autohist(x::AbstractVector{T}, rule::RMG_penB, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ ((n)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# RMG penB
struct RMG_penR <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
end

"""
    RMG_penB(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k N_j \\log (N_j/|\\mathcal{I}_j|) - \\log \\binom{k_n-1}{k-1} - k - \\log^{2.5}(k).
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only `SegNeig()` is supported for this rule. See [`DP`](@ref) for further details.

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
RMG_penR(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=SegNeig()) = RMG_penR(grid, maxbins, alg)

function fit_autohist(x::AbstractVector{T}, rule::RMG_penR, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ ((n)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# NML_I
struct NML_I <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
end

"""
    NML_I(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=DP()
    )

A quick-to-evalutate version of the normalized maximum likelihood criterion.

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
\\begin{aligned}
    &\\sum_{j=1}^k N_j\\log \\frac{N_j}{|\\mathcal{I}_j|} - \\frac{k-1}{2}\\log(n/2) - \\log\\frac{\\sqrt{\\pi}}{\\Gamma(k/2)} - n^{-1/2}\\frac{\\sqrt{2}k\\Gamma(k/2)}{3\\Gamma(k/2-1/2)} \\\\
    &- n^{-1}\\left(\\frac{3+k(k-2)(2k+1)}{36} - \\frac{\\Gamma(k/2)^2 k^2}{9\\Gamma(k/2-1/2)^2} \\right)  - \\log \\binom{k_n-1}{k-1}
\\end{aligned}
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only `SegNeig()` is supported for this rule. See [`DP`](@ref) for further details.

# References
This a variant of this criterion first suggested by [Kontkanen and Myllymäki (2007)](https://proceedings.mlr.press/v2/kontkanen07a.html).
"""
NML_I(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=SegNeig()) = NML_I(grid, maxbins, alg)

function fit_autohist(x::AbstractVector{T}, rule::NML_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ ((n)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end

# ------------------------------
# L2CV_I
struct L2CV_I <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
    use_min_length::Bool
end

"""
    L2CV_I(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=OptPart(),
        use_min_length::Bool=false
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a L2 leave-one-out cross-validation criterion,
```math
    \\frac{n+1}{n}\\sum_{j=1}^k \\frac{N_j^2}{|\\mathcal{I}_j|} - 2\\sum_{j=1}^k \\frac{N_j}{|\\mathcal{I}_j|}.
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, `OptPart()` and `SegNeig()` are supported, with the formed algorithm being the default. See [`DP`](@ref) for further details.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.

# References
This approach dates back to [Rudemo (1982)](https://www.jstor.org/stable/4615859).
"""
L2CV_I(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=OptPart(), use_min_length::Bool=false) = L2CV_I(grid, maxbins, alg, use_min_length)

function fit_autohist(x::AbstractVector{T}, rule::L2CV_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# KLCV_I
struct KLCV_I <: AbstractIrregularRule
    grid::Symbol
    maxbins::Union{Int, Symbol}
    alg::AbstractAlgorithm
    use_min_length::Bool
end

"""
    KLCV_I(;
        grid::Symbol=:regular,
        maxbins::Union{Int, Symbol}=:default,
        alg::AbstractAlgorithm=OptPart(),
        use_min_length::Bool=false
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a Kullback-Leibler leave-one-out cross-validation criterion,
```math
    \\sum_{j=1}^k N_j\\log(N_j-1) - \\sum_{j=1}^k N_j\\log |I_j|,
```
where the maximmization is over all partitions with ``N_j \\geq 2`` for all ``j``.

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, `OptPart()` and `SegNeig()` are supported, with the formed algorithm being the default. See [`DP`](@ref) for further details.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.

# References
This approach to irregular histograms was, to the best of our knowledge, first considered in [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).
"""
KLCV_I(; grid::Symbol=:regular, maxbins::Union{Int, Symbol}=:default, alg::AbstractAlgorithm=OptPart(), use_min_length::Bool=false) = KLCV_I(grid, maxbins, alg, use_min_length)

function fit_autohist(x::AbstractVector{T}, rule::KLCV_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end