struct RRH <: AbstractRegularRule
    a::Union{Real, Function}
    logprior::Function
    maxbins::Union{Int, Symbol}
end

"""
    RRH(; a::Union{Real, Function}, logprior::Function, maxbins::Union{Int, Symbol}=:default)
    Knuth(; maxbins::Union{Int, Symbol}=:default)

The random regular histogram criterion for regular histograms.

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

# References
The `Knuth` criterion for histograms was proposed by [Knuth (2019)](https://doi.org/10.1016/j.dsp.2019.102581). The random regular histogram criterion is a generalization.
"""
RRH, Knuth

RRH(; a::Union{Real, Function}=5.0, logprior::Function=k->0.0, maxbins::Union{Int, Symbol}=:default) = RRH(a, logprior, maxbins)

function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}, rule::RRH; support::Tuple{Real,Real}=(-Inf,Inf), closed::Symbol=:right)
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end
    if maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("Maximal number of bins must be positive."))
    end
    xmin, xmax = check_support(x, support[1], support[2])
    n = length(x)
    aₖ = Vector{Float64}(undef, k_max)
    if !isa(a, Function) # create constant function if typeof(a) <: Real
        if a ≤ 0.0
            throw(DomainError("Supplied value of a must be positive."))
        else 
            aₖ[1:end] .= a
        end
    else 
        for k = 1:k_max
            a_curr = a(k)
            aₖ[k] = a_curr
        end
        if min(a_curr) ≤ 0.0
            throw(DomainError("Supplied function a(k) must return strictly positive values."))
        end
    end
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = logposterior_k(N, k, ones(k)/k, aₖ[k], n, logprior)
    end
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    a_opt = aₖ[k_opt]
    dens = k_opt/(xmax-xmin) * (N .+ a_opt/k_opt) / (a_opt + n) # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, ifelse(a_opt > 0.0, a_opt, NaN))
    return h
end

struct Knuth <: AbstractRegularRule
    maxbins::Union{Int, Symbol}
end
Knuth(; maxbins::Union{Int, Symbol}=:default) = RRH(a=k->0.5*k, logprior=k->0.0, maxbins)


struct AIC <: AbstractRegularRule
    maxbins::Union{Int, Symbol}
end

"""
    AIC(; maxbins::Union{Int, Symbol}=:default)

AIC criterion for regular histograms.

The number ``k`` of bins is chosen as the maximizer of the penalized log-likelihood,
```math
    n\\log (k) + \\sum_{j=1}^k N_j \\log (N_j/n) - k,
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.

# References
The aic criterion was proposed by [Taylor (1987)](https://doi.org/10.1093/biomet/74.3.636) for histograms.
"""
AIC(; maxbins::Union{Int, Symbol}=:default) = AIC(maxbins)

function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}, rule::AIC; support::Tuple{Real,Real}=(-Inf,Inf), closed::Symbol=:right)
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end
    if maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("Maximal number of bins must be positive."))
    end
    xmin, xmax = check_support(x, support[1], support[2])
    n = length(x)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_AIC(N, k, n) # Note: negative of AIC is computed
    end
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, ifelse(a_opt > 0.0, a_opt, NaN))
    return h
end


struct BIC <: AbstractRegularRule
    maxbins::Union{Int, Symbol}
end

"""
    BIC(; maxbins::Union{Int, Symbol}=:default)

BIC criterion for regular histograms.

The number ``k`` of bins is chosen as the maximizer of the penalized log-likelihood,
```math
    n\\log (k) + \\sum_{j=1}^k N_j \\log (N_j/n) - \\frac{k}{2}\\log(n).
```
where ``n`` is the sample size.

# Keyword arguments
- `maxbins`: Maximal number of bins for which the above criterion is evaluated. Defaults to `maxbins=:default`, which sets maxbins to the ceil of `min(1000, 4n/log(n)^2)`.
"""
BIC(; maxbins::Union{Int, Symbol}=:default) = BIC(maxbins)

function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}, rule::BIC; support::Tuple{Real,Real}=(-Inf,Inf), closed::Symbol=:right)
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end
    if maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("Maximal number of bins must be positive."))
    end
    xmin, xmax = check_support(x, support[1], support[2])
    n = length(x)
    criterion = Vector{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    z = @. (x - xmin) / (xmax - xmin)         # Scale data to the interval [0,1]
    Threads.@threads for k = 1:k_max
        N = bin_regular(z, 0.0, 1.0, k, closed == :right)
        @inbounds criterion[k] = compute_BIC(N, k, n) # Note: negative of BIC is computed
    end
    N = bin_irregular_int(x, LinRange(xmin, xmax, k_opt+1), closed == :right)
    dens = k_opt/(xmax-xmin) * N / n # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, ifelse(a_opt > 0.0, a_opt, NaN))
    return h
end
