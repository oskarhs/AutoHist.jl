# ------------------------------
# Random regular histogram 
struct RRH{T} <: AbstractRegularRule
    a::T
    logprior::Function
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    RRH(;
        a::Union{Real, Function}    = 5.0,
        logprior::Function          = k->0.0,
        maxbins::Union{Int, Symbol} = :default
    )
    Knuth(; maxbins::Union{Int, Symbol} = :default)

The random regular histogram criterion.

The number ``k`` of bins is chosen as the maximizer of the marginal log-posterior,
```math
   n\\log (k) + \\sum_{j=1}^k \\big\\{\\log \\Gamma(a_j + N_j) - \\log \\Gamma(a_j)\\big\\} + \\log p_n(k).
```
Here ``p_n(k)`` is the prior distribution on the number ``k`` of bins, which can be controlled by
supplying a function to the `logprior` keyword argument. 
The default value is ``p_n(k) \\propto 1``. Here, ``a_j = a/k``, for a scalar ``a > 0``,
possibly depending on ``k``. The value of ``a`` can be set by supplying a fixed, positive scalar or
a function ``a(k)`` to the keyword argument `a`. The default value is `a=5.0` for `RRH()`.

The rule `Knuth()` is a special case of the `RRH` criterion, which corresponds to the particular
choices ``a_j = 0.5`` and ``p_n(k)\\propto 1``.

# Keyword arguments
- `a`: Specifies Dirichlet concentration parameter in the Bayesian histogram model. Can either be a fixed positive number or a function computing aₖ for different values of k. Defaults to `5.0` if not supplied.
- `logprior`: Unnormalized logprior distribution on the number ``k`` of bins. Defaults to a uniform prior, e.g. `logprior(k) = 0` for all `k`.
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = RRH(a = k->0.5*k, logprior = k->0.0);

julia> h = fit(AutomaticHistogram, x, rule)
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04950495049504951, 0.2079207920792079, 0.5643564356435643, 1.0, 1.495049504950495, 1.8712871287128714, 1.9702970297029703, 1.6732673267326732, 0.9603960396039604, 0.2079207920792079]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: 5.0

julia> h == fit(AutomaticHistogram, x, Knuth())
true
```

# References
The `Knuth` criterion for histograms was proposed by [Knuth (2019)](https://doi.org/10.1016/j.dsp.2019.102581). The random regular histogram criterion is a generalization.
"""
RRH, Knuth

function RRH(; a::Union{Real, Function}=5.0, logprior::Function=k->0.0, maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return RRH(a, logprior, 0, true)
    else
        return RRH(a, logprior, maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::RRH, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    a_vec = Vector{Float64}(undef, k_max)
    if !isa(rule.a, Function) # create constant function if typeof(a) <: Real
        if rule.a ≤ 0.0
            throw(DomainError("Supplied value of a must be positive."))
        else 
            a_vec[1:end] .= rule.a
        end
    else 
        for k = 1:k_max
            a_curr = rule.a(k)
            a_vec[k] = a_curr
        end
        if minimum(a_vec) ≤ 0.0
            throw(DomainError("Supplied function a(k) must return strictly positive values."))
        end
    end
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = logposterior_k(N, k, ones(k)/k, a_vec[k], n, rule.logprior)
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    a_opt = a_vec[k_opt]
    dens = k_opt/(xmax-xmin) * (N .+ a_opt/k_opt) / (a_opt + n) # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, a_opt)
    return h
end


# ------------------------------
# Knuth's rule
Knuth(; maxbins::Union{Int, Symbol}=:default) = RRH(; a=k->0.5*k, logprior=k->0.0, maxbins=maxbins)


# ------------------------------
# AIC histogram
struct AIC <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    AIC(; maxbins::Union{Int, Symbol} = :default)

AIC criterion for regular histograms.

The number ``k`` of bins is chosen as the maximizer of the penalized log-likelihood,
```math
    n\\log (k) + \\sum_{j=1}^k N_j \\log (N_j/n) - k,
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, AIC())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04, 0.2, 0.56, 1.0, 1.5, 1.88, 1.98, 1.68, 0.96, 0.2]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: NaN
```

# References
The aic criterion was proposed by [Taylor (1987)](https://doi.org/10.1093/biomet/74.3.636) for histograms.
"""
function AIC(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return AIC(0, true)
    else
        return AIC(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::AIC, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_AIC(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# BIC histogram
struct BIC <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    BIC(; maxbins::Union{Int, Symbol} = :default)

BIC criterion for regular histograms.

The number ``k`` of bins is chosen as the maximizer of the penalized log-likelihood,
```math
    n\\log (k) + \\sum_{j=1}^k N_j \\log (N_j/n) - \\frac{k}{2}\\log(n).
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, BIC())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 9)
density: [0.048, 0.336, 0.816, 1.44, 1.904, 1.904, 1.264, 0.288]
counts: [3, 21, 51, 90, 119, 119, 79, 18]
type: regular
closed: right
a: NaN
```
"""
function BIC(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return BIC(0, true)
    else
        return BIC(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::BIC, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_BIC(N, k, n) # Note: negative of BIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end

# ------------------------------
# Birgé-Rozenholc histogram
struct BR <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    BR(; maxbins::Union{Int, Symbol} = :default)

Birgé-Rozenholc criterion for regular histograms.

The number ``k`` of bins is chosen as the maximizer of the penalized log-likelihood,
```math
    n\\log (k) + \\sum_{j=1}^k N_j \\log (N_j/n) - k - \\log^{2.5} (k),
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, BR())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04, 0.2, 0.56, 1.0, 1.5, 1.88, 1.98, 1.68, 0.96, 0.2]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: NaN
```

# References
This criterion was proposed by [Birgé and Rozenholc (2006)](https://doi.org/10.1051/ps:2006001).
"""
function BR(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return BR(0, true)
    else
        return BR(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::BR, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_BR(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# MDL regular histogram
struct MDL <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    MDL(; maxbins::Union{Int, Symbol} = :default)

MDL criterion for regular histograms.

The number ``k`` of bins is chosen as the minimizer of an encoding length of the data, and is equivalent to the maximizer of
```math
    n\\log(k) + \\sum_{j=1}^k \\big(N_j-\\frac{1}{2}\\big)\\log\\big(N_j-\\frac{1}{2}\\big) - \\big(n-\\frac{k}{2}\\big)\\log\\big(n-\\frac{k}{2}\\big) - \\frac{k}{2}\\log(n),
```
where ``n`` is the sample size and the maximmization is over all regular partitions with ``N_j \\geq 1`` for all ``j``.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, MDL())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04, 0.2, 0.56, 1.0, 1.5, 1.88, 1.98, 1.68, 0.96, 0.2]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: NaN
```

# References
The minimum description length principle was first applied to histogram estimation by [Hall and Hannan (1988)](https://doi.org/10.1093/biomet/75.4.705).
"""
function MDL(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return MDL(0, true)
    else
        return MDL(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::MDL, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_MDL(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# NML regular histogram
struct NML_R <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    NML_R(; maxbins::Union{Int, Symbol} = :default)

NML_R criterion for regular histograms.

The number ``k`` of bins is chosen by maximizing a penalized likelihood,
```math
\\begin{aligned}
    &\\sum_{j=1}^k N_j\\log \\frac{N_j}{|\\mathcal{I}_j|} - \\frac{k-1}{2}\\log(n/2) - \\log\\frac{\\sqrt{\\pi}}{\\Gamma(k/2)} - n^{-1/2}\\frac{\\sqrt{2}k\\Gamma(k/2)}{3\\Gamma(k/2-1/2)} \\\\
    &- n^{-1}\\left(\\frac{3+k(k-2)(2k+1)}{36} - \\frac{\\Gamma(k/2)^2 k^2}{9\\Gamma(k/2-1/2)^2} \\right).
\\end{aligned}
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, NML_R())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 24)
density: [0.046, 0.0, 0.138, 0.184, 0.368, 0.506, 0.69, 0.874, 1.104, 1.334  …  1.978, 1.978, 1.978, 1.84, 1.61, 1.334, 0.966, 0.644, 0.276, 0.046]
counts: [1, 0, 3, 4, 8, 11, 15, 19, 24, 29  …  43, 43, 43, 40, 35, 29, 21, 14, 6, 1]
type: regular
closed: right
a: NaN
```

# References
This is a regular variant of the normalized maximum likelihood criterion considered by [Kontkanen and Myllymäki (2007)](https://proceedings.mlr.press/v2/kontkanen07a.html).
"""
function NML_R(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return NML_R(0, true)
    else
        return NML_R(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::NML_R, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_NML(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# L2CV for regular histograms
struct L2CV_R <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    L2CV_R(; maxbins::Union{Int, Symbol} = :default)

L2 cross-validation criterion for regular histograms.

The number ``k`` of bins is chosen by maximizing a leave-one-out L2 cross-validation criterion,
```math
    -2k + k\\frac{n+1}{n^2}\\sum_{j=1}^k N_j^2.
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, L2CV_R())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04, 0.2, 0.56, 1.0, 1.5, 1.88, 1.98, 1.68, 0.96, 0.2]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: NaN
```

# References
This approach to histogram density estimation was first considered by [Rudemo (1982)](https://www.jstor.org/stable/4615859).
"""
function L2CV_R(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return L2CV_R(0, true)
    else
        return L2CV_R(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::L2CV_R, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_L2CV(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end

# ------------------------------
# KLCV for regular histograms
struct KLCV_R <: AbstractRegularRule
    maxbins::Int
    use_default_maxbins::Bool
end

"""
    KLCV_R(; maxbins::Union{Int, Symbol} = :default)

Kullback-Leibler cross-validation criterion for regular histograms.

The number ``k`` of bins is chosen by maximizing a leave-one-out Kullback-Leibler cross-validation criterion,
```math
    n\\log(k) + \\sum_{j=1}^k N_j\\log (N_j-1),
```
where ``n`` is the sample size and the maximmization is over all regular partitions with ``N_j \\geq 2`` for all ``j``.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, KLCV_R())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 11)
density: [0.04, 0.2, 0.56, 1.0, 1.5, 1.88, 1.98, 1.68, 0.96, 0.2]
counts: [2, 10, 28, 50, 75, 94, 99, 84, 48, 10]
type: regular
closed: right
a: NaN
```

# References
This approach was first studied by [Hall (1990)](https://doi.org/10.1007/BF01203164).
"""
function KLCV_R(; maxbins::Union{Int, Symbol}=:default)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    if maxbins == :default
        return KLCV_R(0, true)
    else
        return KLCV_R(maxbins, false)
    end
end

function fit_autohist(x::AbstractVector{T}, rule::KLCV_R, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k_max = ifelse(rule.use_default_maxbins, min(ceil(Int, 4.0*n / log(n)^2), 1000), rule.maxbins)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_KLCV(N, k, n) # Note: negative of AIC is computed
    end
    k_opt = argmax(criterion)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# Sturges' rule
"""
    Sturges()

Sturges' rule for regular histograms.

The number ``k`` of bins is chosen as
```math
    k = \\lceil \\log_2(n) \\rceil + 1,
```
where ``n`` is the sample size.

This is the default procedure used by the `hist()` function in base R.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, Sturges())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 10)
density: [0.054, 0.252, 0.666, 1.206, 1.71, 1.98, 1.782, 1.098, 0.252]
counts: [3, 14, 37, 67, 95, 110, 99, 61, 14]
type: regular
closed: right
a: NaN
```

# References
This classical rule is due to [Sturges (1926)](https://doi.org/10.1080/01621459.1926.10502161).
"""
struct Sturges <: AbstractRegularRule end

function fit_autohist(x::AbstractVector{T}, rule::Sturges, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k = ceil(Int64, log2(n))
    N = bin_irregular_int(x, LinRange(xmin, xmax, k+1), closed == :right)
    dens = k/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# Scott's rule
"""
    Scott()

Scott's rule for regular histograms.

The number ``k`` of bins is computed according to the formula
```math
    k = \\big\\lceil \\hat{\\sigma}^{-1}(24\\sqrt{\\pi})^{-1/3}n^{1/3}\\big\\rceil,
```
where ``\\hat{\\sigma}`` is the sample standard deviation and ``n`` is the sample size.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, Scott())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 14)
density: [0.026, 0.13, 0.338, 0.624, 0.988, 1.378, 1.716, 1.924, 2.002, 1.768, 1.3, 0.676, 0.13]
counts: [1, 5, 13, 24, 38, 53, 66, 74, 77, 68, 50, 26, 5]
type: regular
closed: right
a: NaN
```

# References
This classical rule is due to [Scott (1979)](https://doi.org/10.1093/biomet/66.3.605).
"""
struct Scott <: AbstractRegularRule end

function fit_autohist(x::AbstractVector{T}, rule::Scott, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    h_scott = std(x)*((24.0*sqrt(π))/n)^(1.0/3.0)
    k = ceil(Int64, (xmax-xmin)/h_scott)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k+1), closed == :right)
    dens = k/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# Freedman and Diaconis' rule
"""
    FD()

Freedman and Diaconis' rule for regular histograms.

The number ``k`` of bins is computed according to the formula
```math
    k = \\big\\lceil\\frac{n^{1/3}}{2\\text{IQR}(\\boldsymbol{x})}\\big\\rceil,
```
where ``\\text{IQR}(\\boldsymbol{x})`` is the sample interquartile range and ``n`` is the sample size.

This is the default procedure used by the `histogram()` function in `Plots.jl`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x, FD())
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 16)
density: [0.03, 0.09, 0.24, 0.48, 0.78, 1.08, 1.44, 1.71, 1.92, 2.01, 1.89, 1.59, 1.08, 0.54, 0.12]
counts: [1, 3, 8, 16, 26, 36, 48, 57, 64, 67, 63, 53, 36, 18, 4]
type: regular
closed: right
a: NaN
```

# References
This rule dates back to [Freedman and Diaconis (1982)](https://doi.org/10.1007/BF01025868).
"""
struct FD <: AbstractRegularRule end

function fit_autohist(x::AbstractVector{T}, rule::FD, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    h_fd = 2.0*iqr(x)/n^(1.0/3.0)
    k = ceil(Int64, (xmax - xmin)/h_fd)
    N = bin_irregular_int(x, LinRange(xmin, xmax, k+1), closed == :right)
    dens = k/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k+1), dens, N, :regular, closed, NaN)
    return h
end


# ------------------------------
# Wand's rule (see wand_num_bins.jl for the implementation)
struct Wand <: AbstractRegularRule
    level::Int
    scalest::Symbol
end

"""
    Wand(;
        level::Int      = 2,
        scalest::Symbol = :minim
    )

Wand's rule for regular histograms.

A more sophisticated version of Scott's rule, Wand's rule proceeds by determining the bin width ``h`` as
```math
    h = \\Big(\\frac{6}{\\hat{C}(f_0) n}\\Big)^{1/3},
```
where ``\\hat{C}(f_0)`` is an estimate of the functional ``C(f_0) = \\int \\{f_0'(x)\\}^2\\, \\text{d}x``. The corresponding number of bins is ``k = \\lceil h^{-1}\\rceil``.

# Keyword arguments
`level`: The `level` keyword controls the number of stages of functional estimation used to compute ``\\hat{C}``, and can take values `0, 1, 2, 3, 4, 5`, with the default value being `level=2`. The choice `level=0` corresponds to a varation on Scott's rule, with a custom scale estimate.
`scalest`: Estimate of scale parameter. Possible choices are `:minim` `:stdev` and `:iqr`. The latter two use sample standard deviation or the sample interquartile range, respectively, to estimate the scale. The default choice `:minim` uses the minimum of the above estimates.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 500)) .^(1/3)).^(1/3);

julia> rule = Wand(scalest=:stdev, level=5);

julia> fit(AutomaticHistogram, x, rule)
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 13)
density: [0.024, 0.144, 0.408, 0.72, 1.128, 1.536, 1.872, 1.992, 1.848, 1.416, 0.744, 0.168]
counts: [1, 6, 17, 30, 47, 64, 78, 83, 77, 59, 31, 7]
type: regular
closed: right
a: NaN
```

# References
The full details on this method are given in [Wand (1997)](https://doi.org/10.2307/2684697).
"""
function Wand(; level::Int=2, scalest::Symbol=:minim)
    if !(scalest in [:minim, :stdev, :iqr])     # check that supplied scale-estimate is a valid option
        throw(ArgumentError("Supplied scalest value, :$(scalest), is not supported. Use one of :minim, :stdev or :iqr."))
    end
    if !(level in [0, 1, 2, 3, 4, 5])           # check that supplied level is a valid option
        throw(ArgumentError("Supplied level, $(level), is not supported. Use one of 0, 1, 2, 3, 4 or 5."))
    end
    return Wand(level, scalest)
end

function fit_autohist(x::AbstractVector{T}, rule::Wand, xmin::T, xmax::T, closed::Symbol) where {T <: Real}
    n = length(x)
    k = wand_num_bins(x, rule.level, rule.scalest, 401, (xmin, xmax))
    N = bin_irregular_int(x, LinRange(xmin, xmax, k+1), closed == :right)
    dens = k/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k+1), dens, N, :regular, closed, NaN)
    return h
end