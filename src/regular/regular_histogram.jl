"""
    histogram_regular(x::AbstractVector{<:Real}; rule::Symbol=:bayes, closed::Symbol=:right, maxbins::Int=1000, support::Tuple{Real,Real}=(-Inf,Inf), logprior::Function=k->0.0, a::Union{Real,Function}=1.0, level::Int=2, scalest::Symbol=:minim)

Create a regular histogram based on an asymptotic risk estimate, or optimization criteria from Bayesian probability, penalized likelihood or LOOCV.
Returns a AutomaticHistogram object with regular bins, with the optimal bin number corresponding to the supplied criterion.

...
# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Defaults to the method Bayesian method of Simensen et al. (2025)
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `maxbins`: The maximal number of bins to be considered by the optimization criterion. Ignored if the specified argument is not a positive integer. Defaults to `maxbins=1000`
- `support`: Tuple specifying the the support of the histogram estimate. If the first element is `-Inf`, then `minimum(x)` is taken as the leftmost cutpoint. Likewise, if the second element is `Inf`, then the rightmost cutpoint is `maximum(x)`. Default value is `(-Inf, Inf)`, which estimates the support of the data.
- `logprior`: Unnormalized logprior distribution of the number k of bins. Only used in the case where the supplied rule is `:bayes`. Defaults to a uniform prior.
- `a`: Specifies Dirichlet concentration parameter in the Bayesian histogram model. Can either be a fixed positive number or a function computing aₖ for different values of k. Defaults to `1.0` if not supplied. Uses default if suppled value is negative.
- `scalest`: Estimate of scale parameter used in computing Wands' rule. Only used if `rule` is set to `:wand`. Possible values are `:minim` `:stdved` and `:iqr`. Default value is `scalest=:minim`.
- `level`: Specifies the level used for the Kernel functional estimate in Wands' rule. Only used if `rule==:wand`. Possible values are 0,1,2,3,4 and 5. Default value is `level=2`.

# Returns
- `h`: AutomaticHistogram object with weights corresponding to densities, e.g. `isdensity` is set to true.

# Examples
```
julia> x = [0.037, 0.208, 0.189, 0.656, 0.45, 0.846, 0.986, 0.751, 0.249, 0.447]
julia> h1 = histogram_regular(x)
julia> h2 = histogram_regular(x; logprior=k->-log(k), a=k->0.5*k)
```
...
"""
function histogram_regular( x::AbstractVector{<:Real}; rule::Symbol=:bayes, closed::Symbol=:right, maxbins::Int=1000,
                            support::Tuple{Real,Real}=(-Inf,Inf), logprior::Function=k->0.0, a::Union{Real,Function}=1.0,
                            scalest::Symbol=:minim, level::Int=2)
    if !(rule in [:aic, :bic, :br, :bayes, :mdl, :klcv, :nml, :l2cv, :sturges, :fd, :scott, :wand])
        throw(ArgumentError("The supplied rule, :$rule, is not supported. The rule kwarg must be one of :aic, :bic, :br, :bayes, :mdl, :klcv, :nml, :l2cv, :sturges, :fd, :scott or :wand"))
    end
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end
    n = length(x)

    if maxbins < 1          # maximal number of bins must be positive
        throw(DomainError("Maximal number of bins must be positive."))
    end

    xmin, xmax = extrema(x)

    if support[1] > -Inf   # estimate lower bound of support if unknown,
        if xmin > support[1]
            xmin = support[1]   # use known lower bound
        else 
            throw(DomainError("The supplied lower bound is greater than the smallest value of the sample."))
        end
    end
    if support[2] < Inf
        if xmax < support[2]
            xmax = support[2]   # use known upper bound
        else 
            throw(DomainError("The supplied upper bound is smaller than the smallest value of the sample."))
        end
    end

    k_opt = if rule == :sturges # Sturges' rule
        k = min(ceil(Int64, log2(length(n))) + 1, maxbins)
        k
    elseif rule == :fd # Freedman and Diaconis' rule
        h_fd = 2.0*iqr(x)/n^(1.0/3.0)
        k = min(ceil(Int64, (xmax - xmin)/h_fd), maxbins)
        k
    elseif rule == :scott # Scott's normal reference rule
        h_scott = std(x)*((24.0*sqrt(π))/n)^(1.0/3.0)
        k = min(ceil(Int64, (xmax-xmin)/h_scott), maxbins)
        k
    elseif rule == :wand # Wand's rule (using 2 steps)
        if !(scalest in [:minim, :stdev, :iqr])     # check that supplied scale-estimate is a valid option
            throw(ArgumentError("Supplied scalest value, :$scalest, is not supported. Use one of :minim, :stdev or :iqr."))
        end
        if !(level in [0, 1, 2, 3, 4, 5])           # check that supplied level is a valid option
            throw(ArgumentError("Supplied level, $level, is not supported. Use one of 0, 1, 2, 3, 4 or 5."))
        end
        k = min(wand_num_bins(x, level, scalest, 401, (xmin, xmax)), maxbins) # add functionality for specifying level later
        k
    else 
        k_max = min(ceil(Int, 4.0*n / log(n)^2), maxbins)

        if rule == :bayes
            aₖ = Array{Float64}(undef, k_max)
            if !isa(a, Function) # create constant function if typeof(a) <: Real
                if a ≤ 0.0
                    throw(DomainError("Supplied value of a must be positive."))
                else 
                    aₖ[1:end] .= a
                end
            else 
                for k = 1:k_max
                    a_curr = a(k)
                    if a_curr ≤ 0.0
                        throw(DomainError("Supplied function a(k) must return strictly positive values."))
                    else
                        aₖ[k] = a_curr
                    end
                end
            end
        end

        criterion = Array{Float64}(undef, k_max) # Criterion to be maximized depending on the specified rule
    
        # Scale data to the interval [0,1]:
        z = @. (x - xmin) / (xmax - xmin)
    
        if rule == :aic
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_AIC(N, k, n) # Note: negative of AIC is computed
            end
        elseif rule == :bic
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_BIC(N, k, n) # Note: negative of BIC is computed
            end
        elseif rule == :br
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_BR(N, k, n)
            end
        elseif rule == :bayes
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = logposterior_k(N, k, ones(k)/k, aₖ[k], n, logprior)
            end
        elseif rule == :mdl
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_MDL(N, k, n)
            end
        elseif rule == :klcv
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_KLCV(N, k, n)
            end
        elseif rule == :l2cv
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_L2CV(N, k, n)
            end
        elseif rule == :nml
            Threads.@threads for k = 1:k_max
                N = bin_regular(z, 0.0, 1.0, k, closed == :right)
                @inbounds criterion[k] = compute_NML(N, k, n)
            end
        end
        k_opt = argmax(criterion)
        k_opt
    end

    # Create a StatsBase.Histogram object with the chosen number of bins
    #N = convert.(Int64, bin_regular(x, xmin, xmax, k_opt, closed == :right))
    N = bin_regular_int(x, xmin, xmax, k_opt, closed == :right)
    a_opt = 0.0
    if rule == :bayes
        a_opt = aₖ[k_opt]
    end
    dens = k_opt/(xmax-xmin) * (N .+ a_opt/k_opt) / (a_opt + n) # Estimated density
    h = AutomaticHistogram(LinRange(xmin, xmax, k_opt+1), dens, N, :regular, closed, ifelse(a_opt > 0.0, a_opt, NaN))
    return h

    #H_opt = Histogram(LinRange(xmin, xmax, k_opt+1), dens, closed, true)
    #return H_opt
end