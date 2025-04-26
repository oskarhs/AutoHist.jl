"""
    histogram_regular(x::AbstractVector{<:Real}; rule::Str="bayes", right::Bool=true, maxbins::Int=1000, support::Tuple{Real,Real}=(-Inf,Inf), logprior::Function=k->0.0, a::Union{Real,Function}=1.0)

Create a regular histogram based on optimization criterion from Bayesian probability, penalized likelihood or LOOCV.
Returns a StatsBase.Histogram object with regular bins, with the optimal bin number corresponding to the supplied criterion.

...
# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Defaults to the method Bayesian method of Simensen et al. (2025)
- `right`: Boolean indicating whether the drawn intervals should be right-inclusive or not. Defaults to `true`.
- `maxbins`: The maximal number of bins to be considered by the optimization criterion. Ignored if the specified argument is not a positive integer. Defaults to `maxbins=1000`
- `support`: Tuple specifying the the support of the histogram estimate. If the first element is -Inf, then `minimum(x)` is taken as the leftmost cutpoint. Likewise, if the second element is `Inf`, then the rightmost cutpoint is `maximum(x)`. Default value is `(-Inf, Inf)`, which estimates the support of the data.
- `logprior`: Unnormalized logprior distribution of the number k of bins. Only used in the case where the supplied rule is `"bayes"`. Defaults to a uniform prior.
- `a`: Specifies Dirichlet concentration parameter in the Bayesian histogram model. Can either be a fixed positive number or a function computing aₖ for different values of k. Defaults to `1.0` if not supplied. Uses default if suppled value is negative.

# Returns
- `H`: StatsBase.Histogram object with weights corresponding to densities, e.g. `:isdensity` is set to true.

# Examples
```
julia> x = [0.037, 0.208, 0.189, 0.656, 0.45, 0.846, 0.986, 0.751, 0.249, 0.447]
julia> H1 = histogram_regular(x)
julia> H2 = histogram_regular(x; logprior=k->-log(k), a=k->0.5*k)
```
...
"""
function histogram_regular(x::AbstractVector{<:Real}; rule::String="bayes", right::Bool=true, maxbins::Int=1000, support::Tuple{Real,Real}=(-Inf,Inf), logprior::Function=k->0.0, a::Union{Real,Function}=1.0)
    rule = lowercase(rule)
    if !(rule in ["aic", "bic", "br", "bayes", "mdl", "sc", "klcv", "nml", "l2cv"])
        rule = "bayes"
    end
    if rule == "bayes"
        if !isa(a, Function) # create constant function if typeof(a) <: Real
            if a ≤ 0.0
                a_func = k -> 1.0
            else 
                a_func = k -> a
            end
        else 
            a_func = k -> ifelse(a(k) > 0.0, a(k), 1.0)
        end
    end

    n = length(x)
    if maxbins < 1
        maxbins = 10^3 # Default maximal number of bins
    end
    k_max = min(ceil(Int, 4.0*n / log(n)^2), maxbins)

    criterion = Array{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule

    # Scale data to the interval [0,1]:
    if support[1] == -Inf
        xmin = minimum(x) 
    else
        xmin = support[1]
    end
    if support[2] == Inf
        xmax = maximum(x)
    else 
        xmax = support[2]
    end
    z = @. (x - xmin) / (xmax - xmin)
    N = Vector{Float64}(undef, k_max)
    if rule == "aic"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_AIC(view(N, 1:k), k, n) # Note: negative of AIC is computed
        end
    elseif rule == "bic"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_BIC(view(N, 1:k), k, n) # Note: negative of BIC is computed
        end
    elseif rule == "br"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_BR(view(N, 1:k), k, n)
        end
    elseif rule == "bayes"
        for k = 1:k_max
            aₖ = a_func(k)
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = logposterior_k(view(N, 1:k), k, ones(k)/k, aₖ, n, logprior)
        end
    elseif rule == "mdl"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_MDL(view(N, 1:k), k, n)
        end
    elseif rule == "sc"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_SC(view(N, 1:k), k, n)
        end
    elseif rule == "klcv"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_KLCV(view(N, 1:k), k, n)
        end
    elseif rule == "l2cv"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_L2CV(view(N, 1:k), k, n)
        end
    elseif rule == "nml"
        for k = 1:k_max
            N[1:k] = bin_regular(z, 0.0, 1.0, k, right)
            criterion[k] = compute_NML(view(N, 1:k), k, n)
        end
    end

    # Create a StatsBase.Histogram object with the chosen number of bins
    k_opt = argmax(criterion)
    N = bin_regular(z, 0.0, 1.0, k_opt, right)
    if right
        H_opt = Histogram(LinRange(xmin, xmax, k_opt+1), N, :right, true)
    else
        H_opt = Histogram(LinRange(xmin, xmax, k_opt+1), N, :left, true)
    end
    if rule == "bayes"
        aₖ = a_func(k_opt)
        H_opt.weights = k_opt/(xmax-xmin) * (N .+ aₖ/k_opt) / (aₖ + n) # Estimated density
    else
        H_opt.weights = k_opt/(xmax-xmin) * N / n # Estimated density
    end
    return H_opt
end