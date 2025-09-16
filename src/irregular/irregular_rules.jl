# ------------------------------
# Random irregular histogram 
struct RIH{T<:Real, A<:AbstractAlgorithm} <: AbstractIrregularRule
    a::T
    logprior::Function
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
end

"""
    RIH(;
        a::Real                    = 5.0,
        logprior::Function         = k -> 0.0,
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = SegNeig()
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
- `alg`: Algorithm used to fit the model. Currently, only [`SegNeig`](@ref) is supported for this rule.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = RIH(a = 5.0, logprior = k-> -log(k), grid = :data);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.18220071105959446, 0.3587941358096334, 0.8722292888743843, 1.0]
density: [0.11858322346056327, 0.6490600487586273, 1.6066011289666577, 0.30436411114439915]
counts: [10, 57, 414, 19]
type: irregular
closed: right
a: 5.0
```

# References
This approach to irregular histograms first appeared in [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).
"""
function RIH(; a::Real=5.0, logprior::Function=k->0.0, grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=SegNeig())
    if typeof(alg) != SegNeig
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule RIH. Only the SegNeig algorithm is supported for this rule."))
    end
    if a ≤ 0.0
        throw(DomainError("Supplied value of a must be positive."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option: $(grid). Please select one of :data, :regular or :quantile."))
    end
    if maxbins == :default
        return RIH(a, logprior, grid, 0, true, alg)
    else
        return RIH(a, logprior, grid, maxbins, false, alg)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::RIH, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
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
struct RMG_penA{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
end

"""
    RMG_penA(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k N_j \\log (N_j/|\\mathcal{I}_j|) - \\log \\binom{k_n-1}{k-1} - k - 2\\log(k) - \\sqrt{2(k-1)\\Big[\\log \\binom{k_n-1}{k-1}+ 2\\log(k)\\Big]}.
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only [`SegNeig`](@ref) is supported for this rule.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = RMG_penA(grid = :data);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.18875598171056715, 0.3644223879547405, 0.8696799410193466, 1.0]
density: [0.116552597701164, 0.6717277510418354, 1.6229346697072502, 0.30693663211078037]
counts: [11, 59, 410, 20]
type: irregular
closed: right
a: NaN
```

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
function RMG_penA(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=SegNeig())
    if typeof(alg) != SegNeig
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule RMG_penA. Only the SegNeig algorithm is supported for this rule."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if maxbins == :default
        return RMG_penA(grid, 0, true, alg)
    else
        return RMG_penA(grid, maxbins, false, alg)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::RMG_penA, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# RMG penB
struct RMG_penB{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
end

"""
    RMG_penB(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k N_j \\log (N_j/|\\mathcal{I}_j|) - \\log \\binom{k_n-1}{k-1} - k - \\log^{2.5}(k).
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only [`SegNeig`](@ref) is supported for this rule.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = RMG_penB(grid = :data);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.1948931612779725, 0.375258352661302, 0.8268306249022703, 0.9222490305512866, 1.0]
density: [0.12314439276691318, 0.7096713008662634, 1.6962954704872724, 0.7545714006671028, 0.12861575966067232]
counts: [12, 64, 383, 36, 5]
type: irregular
closed: right
a: NaN
```

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
function RMG_penB(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=SegNeig())
    if typeof(alg) != SegNeig
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule RMG_penB. Only the SegNeig algorithm is supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return RMG_penB(grid, 0, true, alg)
    else
        return RMG_penB(grid, maxbins, false, alg)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::RMG_penB, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# RMG penR
struct RMG_penR{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
end

"""
    RMG_penR(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = SegNeig()
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized log-likelihood,
```math
    \\sum_{j=1}^k \\big\\{N_j \\log (N_j/|\\mathcal{I}_j|) - \\frac{N_j}{2n}\\big\\} - \\log \\binom{k_n-1}{k-1} - \\log^{2.5}(k).
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, only [`SegNeig`](@ref) is supported for this rule.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = RMG_penR(grid = :data);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.18875598171056715, 0.3699070396003733, 0.8285645195146814, 0.9222490305512866, 1.0]
density: [0.116552597701164, 0.6845115973621804, 1.6875338000474953, 0.7471886144834441, 0.12861575966067232]
counts: [11, 62, 387, 35, 5]
type: irregular
closed: right
a: NaN
```

# References
This approach was suggested by [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021).
"""
function RMG_penR(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=SegNeig())
    if typeof(alg) != SegNeig
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule RMG_penR. Only the SegNeig algorithm is supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return RMG_penR(grid, 0, true, alg)
    else
        return RMG_penR(grid, maxbins, false, alg)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::RMG_penR, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# NML_I
struct NML_I{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
end

"""
    NML_I(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = DP()
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
- `alg`: Algorithm used to fit the model. Currently, only [`SegNeig`](@ref) is supported for this rule.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = NML_I(grid = :data);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.18875598171056715, 0.3644223879547405, 0.8696799410193466, 1.0]
density: [0.116552597701164, 0.6717277510418354, 1.6229346697072502, 0.30693663211078037]
counts: [11, 59, 410, 20]
type: irregular
closed: right
a: NaN
```

# References
This a variant of this criterion first suggested by [Kontkanen and Myllymäki (2007)](https://proceedings.mlr.press/v2/kontkanen07a.html).
"""
function NML_I(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=SegNeig())
    if typeof(alg) != SegNeig
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule NML_I. Only the SegNeig algorithm is supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return NML_I(grid, 0, true, alg)
    else
        return NML_I(grid, maxbins, false, alg)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::NML_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end

# ------------------------------
# L2CV_I
struct L2CV_I{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
    use_min_length::Bool
end

"""
    L2CV_I(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = OptPart(),
        use_min_length::Bool       = false
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a L2 leave-one-out cross-validation criterion,
```math
    \\frac{n+1}{n}\\sum_{j=1}^k \\frac{N_j^2}{|\\mathcal{I}_j|} - 2\\sum_{j=1}^k \\frac{N_j}{|\\mathcal{I}_j|}.
```

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, [`OptPart`](@ref) and [`SegNeig`](@ref) are supported for this rule, with the former algorithm being the default.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = L2CV_I(grid = :data, use_min_length=true);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.149647045210915, 0.2499005080461325, 0.3490626376697454, 0.4600140220788484, 0.7765683248449301, 0.8535131937737716, 0.9121099383916996, 0.9560732934980348, 1.0]
density: [0.08018868653963065, 0.3590898407087615, 0.7664216197097746, 1.2798398213438569, 1.8448651460332468, 1.2476465466304287, 0.6826317786220794, 0.2729545998246794, 0.045530388214070294]
counts: [6, 18, 38, 71, 292, 48, 20, 6, 1]
type: irregular
closed: right
a: NaN
```

# References
This approach dates back to [Rudemo (1982)](https://www.jstor.org/stable/4615859).
"""
function L2CV_I(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=OptPart(), use_min_length::Bool=false)
    if !(typeof(alg) in [SegNeig, OptPart])
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule L2CV_I. Only the SegNeig and OptPart algorithms are supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return L2CV_I(grid, 0, true, alg, use_min_length)
    else
        return L2CV_I(grid, maxbins, false, alg, use_min_length)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::L2CV_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# KLCV_I
struct KLCV_I{A<:AbstractAlgorithm} <: AbstractIrregularRule
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
    use_min_length::Bool
end

"""
    KLCV_I(;
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = OptPart(),
        use_min_length::Bool       = false
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a Kullback-Leibler leave-one-out cross-validation criterion,
```math
    \\sum_{j=1}^k N_j\\log(N_j-1) - \\sum_{j=1}^k N_j\\log |\\mathcal{I}_j|,
```
where the maximmization is over all partitions with ``N_j \\geq 2`` for all ``j``.

# Keyword arguments
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, [`OptPart`](@ref) and [`SegNeig`](@ref) are supported for this rule, with the former algorithm being the default.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = KLCV_I(grid = :data, use_min_length=true);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.13888886265725095, 0.23836051747480758, 0.33883651300547, 0.45084951551151237, 0.7900230337213711, 0.8722292888743843, 0.9352920770792058, 1.0]
density: [0.07200001359848368, 0.321699684786505, 0.7165890680628054, 1.2319998295961743, 1.8220762141507212, 1.1191362485586687, 0.5074307830485909, 0.09272434856770642]
counts: [5, 16, 36, 69, 309, 46, 16, 3]
type: irregular
closed: right
a: NaN
```

# References
This approach to irregular histograms was, to the best of our knowledge, first considered in [Simensen et al. (2025)](https://doi.org/10.48550/ARXIV.2505.22034).
"""
function KLCV_I(; grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=OptPart(), use_min_length::Bool=false)
    if !(typeof(alg) in [SegNeig, OptPart])
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule KLCV_I. Only the SegNeig and OptPart algorithms are supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return KLCV_I(grid, 0, true, alg, use_min_length)
    else
        return KLCV_I(grid, maxbins, false, alg, use_min_length)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::KLCV_I, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end


# ------------------------------
# Bayesian Blocks
struct BayesBlocks{T<:Real, A<:AbstractAlgorithm} <: AbstractIrregularRule
    p0::T
    grid::Symbol
    maxbins::Int
    use_default_maxbins::Bool
    alg::A
    use_min_length::Bool
end

"""
    BayesBlocks(;
        p0::Real                   = 0.05,
        grid::Symbol               = :regular,
        maxbins::Union{Int,Symbol} = :default,
        alg::AbstractAlgorithm     = OptPart(),
        use_min_length::Bool       = false
    )

Consists of finding the partition ``\\mathcal{I}`` that maximizes a penalized Poisson likelihood criterion,
```math
    \\sum_{j=1}^k N_j\\log\\big(N_j/|\\mathcal{I}_j|\\big) + k\\big(4 - \\log(73.53 p_0 n^{-0.478})\\big).
```

# Keyword arguments
- `p0`: Hyperparameter controlling the penalty for the number of bins. Must be a number in the open interval ``(0,1)``, with `p0 = 0.05` serving as the default.
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)` if `grid` is `regular` or `quantile`. Ignored if `grid=:data`.
- `alg`: Algorithm used to fit the model. Currently, [`OptPart`](@ref) and [`SegNeig`](@ref) are supported for this rule, with the former algorithm being the default.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = BayesBlocks(grid = :data, use_min_length=true);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram{Vector{Float64}, Vector{Float64}, Vector{Int64}}
breaks: [0.0, 0.1948931612779725, 0.375258352661302, 0.8268306249022703, 0.9222490305512866, 1.0]
density: [0.12314439276691318, 0.7096713008662634, 1.6962954704872724, 0.7545714006671028, 0.12861575966067232]
counts: [12, 64, 383, 36, 5]
type: irregular
closed: right
a: NaN
```

# References
The Bayesian Blocks method was first introduced in [Scargle et al. (2013)](https://doi.org/10.1088/0004-637X/764/2/167).
"""
function BayesBlocks(; p0::Real=0.05, grid::Symbol=:regular, maxbins::Union{Int,Symbol}=:default, alg::AbstractAlgorithm=OptPart(), use_min_length::Bool=false)
    if !(0.0 < p0 < 1.0)
        throw(DomainError(p0, ""))
    end
    if !(typeof(alg) in [SegNeig, OptPart])
        throw(ArgumentError("Algorithm $(typeof(alg)) not supported for rule BayesBlocks. Only the SegNeig and OptPart algorithms are supported for this rule."))
    end
    if !(grid in (:data, :regular, :quantile))
        throw(ArgumentError("Invalid grid option :$(grid). Please select one of :data, :regular or :quantile."))
    end
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1
        throw(DomainError(maxbins, "maxbins bins must be positive."))
    end
    if maxbins == :default
        return BayesBlocks(p0, grid, 0, true, alg, use_min_length)
    else
        return BayesBlocks(p0, grid, maxbins, false, alg, use_min_length)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::BayesBlocks, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    bin_edges_norm = maximize_additive_crit(x, rule, xmin, xmax, closed)
    
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm
    N = bin_irregular_int(x, bin_edges, closed == :right)
    dens = N ./ (n*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, NaN)
    return h
end