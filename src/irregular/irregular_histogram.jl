"""
    histogram_irregular(x::AbstractVector{<:Real}; rule::Str="bayes", grid::String="regular", right::Bool=true, greedy::Bool=true, maxbins::Int=-1, support::Tuple{Real,Real}=(-Inf,Inf), use_min_length::Bool=false, logprior::Function=k->0.0, a::Real=1.0)

Create an irregular histogram based on optimization of a criterion based on Bayesian probability, penalized likelihood or LOOCV.
Returns a tuple where the first argument is a StatsBase.Histogram object, the second the value of the maxinized criterion.

# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Defaults to the Bayesian method of Simensen et al. (2025).
- `grid`: String indicating how the finest possible mesh should be constructed. Options are `"data"`, which uses each unique data point as a grid point, `"regular"` (default) which constructs a fine regular grid, and `"quantile"` which constructs the grid based on the sample quantiles.
- `right`: Boolean indicating whether the drawn intervals should be right-inclusive or not. Defaults to `true`.
- `greedy`: Boolean indicating whether or not the greedy binning strategy of Rozenholc et al. (2006) should be used prior to running the dynamical programming algorithm. Defaults to `true`. The algorithm can be quite slow for large datasets when this keyword is set to `false`.
- `maxbins`: The maximal number of bins to be considered by the optimization criterion, only used if grid is set to "regular" or "quantile". Defaults to `maxbins=min(4*n/log(n)^2, 1000)`. If the specified argument is not a positive integer, the default value is used.
- `support`: Tuple specifying the the support of the histogram estimate. If the first element is -Inf, then `minimum(x)` is taken as the leftmost cutpoint. Likewise, if the second elemen is `Inf`, then the rightmost cutpoint is `maximum(x)`. Default value is `(-Inf, Inf)`, which estimates the support of the data.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.
- `logprior`: Unnormalized logprior distribution for the number k of bins. Defaults to a uniform prior. Only used in when `rule="bayes"`.
- `a`: Dirichlet concentration parameter in the Bayesian irregular histogram model. Set to the default value (1.0) if the supplied value is not a positive real number. Only used when `rule="bayes"`.

# Examples
```
julia> x = [0.037, 0.208, 0.189, 0.656, 0.45, 0.846, 0.986, 0.751, 0.249, 0.447]
julia> H1, criterion1 = histogram_irregular(x)
julia> H2, criterion2 = histogram_irregular(x; grid="quantile", support=(0.0, 1.0), logprior=k->-log(k), a=sqrt(10))
```
...
"""
function histogram_irregular(x::AbstractVector{<:Real}; rule::String="bayes", grid::String="regular", 
                            right::Bool=true, greedy::Bool=true, maxbins::Int=-1, support::Tuple{Real,Real}=(-Inf,Inf),
                            use_min_length::Bool=false, logprior::Function=k->0.0, a::Real=1.0)
    rule = lowercase(rule)
    if !(rule in ["pena", "penb", "penr", "bayes", "klcv", "l2cv", "nml"])
        rule = "bayes" # Set penalty to default
    end
    if rule == "bayes"
        if a ≤ 0.0
            a = 1.0
        end
    end

    grid = lowercase(grid)
    if !(grid in ["data", "regular", "quantile"])
        grid = "regular"
    end

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
    y = @. (x - xmin) / (xmax - xmin)
    n = length(x)

    if grid == "data"
        maxbins = n
    elseif typeof(maxbins) != Int || maxbins ≤ 0
        maxbins = min(n, 1000, ceil(Int, 4.0*n/log(n)^2))
    end

    # Calculate gridpoints (left-open grid, breaks at data points)
    # right == true means to include observation in the right endpoint, i.e. right-closed
    finestgrid = Array{Float64}(undef, maxbins+1)
    N_cum = Array{Float64}(undef, length(finestgrid)) # cumulative cell counts
    N_cum[1] = 0.0
    if grid == "data"
        sort!(y)
        if right
            finestgrid[1] = - eps()
            finestgrid[2] = 0.5*(y[1]+y[2])
            finestgrid[3:n] = y[2:n-1] .+ eps()
            finestgrid[n+1] = 1.0 + eps()
        else
            finestgrid[1] = -eps()
            finestgrid[2:n-1] = y[2:n-1] .- eps()
            finestgrid[n] = 0.5 * (y[n] - y[n-1])
            finestgrid[n+1] = 1.0 + eps()
        end
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, right))
    elseif grid == "regular"
        N_cum[2:end] = cumsum(bin_regular(y, 0.0, 1.0, maxbins, right))
        finestgrid[1:end] = LinRange(0.0, 1.0, maxbins+1)
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0+eps()
    elseif grid == "quantile"
        sort!(y)
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0 + eps()
        if right
            finestgrid[2:end-1] = quantile(y, LinRange(1.0/maxbins, 1.0-1.0/maxbins, maxbins-1); sorted=true) .+ eps()
        else 
            finestgrid[2:end-1] = quantile(y, LinRange(1.0/maxbins, 1.0-1.0/maxbins, maxbins-1); sorted=true) .- eps()
        end
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, right))
    end

    if greedy
        gr_maxbins = min(maxbins, max(floor(Int, (log(n)*n)^(1.0/3.0)), 100))
        grid_ind = greedy_grid(N_cum, finestgrid, maxbins, gr_maxbins)
        mesh = finestgrid[grid_ind]
        k_max = length(mesh) - 1
        # convert grid_ind to array of integers equal to true
        chosen_ind = findall(grid_ind)

        # Update bin counts to the newly constructed grid
        N_cum = N_cum[grid_ind]
    else
        k_max = maxbins
        mesh = finestgrid
    end

    if rule in ["pena", "penb", "nml", "penr"]
        control = Dict{String, Float64}()
    elseif rule == "bayes"
        control = Dict{String, Float64}("a" => a)
    elseif rule in ["klcv", "l2cv"]
        minlength = 0.0
        if use_min_length
            minlength = log(n)^(1.5)/n
        end
        control = Dict{String, Float64}("minlength" => minlength)
    end

    optimal, ancestor = dynamic_algorithm(rule, N_cum, mesh, n, k_max, control)
    psi = Array{Float64}(undef, k_max)
    if rule == "penb"    
        for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
        end
    elseif rule == "bayes"
        for k = 1:k_max
            @inbounds psi[k] = logprior(k) - logabsbinomial(maxbins-1, k-1)[1] + loggamma(a) - loggamma(a + n)
        end
    elseif rule == "penr"
        for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
        end
    elseif rule == "pena"
        for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - 2.0*log(k) -
                    2.0 * sqrt(1.0*0.5*(k-1)*(logabsbinomial(maxbins-1, k-1)[1] + 1.0*log(k)))
        end
    elseif rule == "nml"
        for k = 1:k_max
            @inbounds psi[k] =  -( 0.5*k*log(0.5*n) - loggamma(0.5*k) +
            1.0/sqrt(n) * sqrt(2.0)*k/3.0 * exp(loggamma(0.5*k) - loggamma(0.5*k-0.5)) +
            1.0/n * ((3.0 + k*(k-2.0)*(2.0*k+1.0))/36.0 - k^2/9.0*exp(2.0*loggamma(0.5*k) - 2.0*loggamma(0.5*k-0.5)))
            )
            psi[k] = psi[k] - logabsbinomial(maxbins-1, k-1)[1]
        end
    end # NB! no penalties for klcv and l2cv
    k_opt = argmax(optimal + psi)
    criterion_opt = optimal[k_opt] + psi[k_opt]

    bin_edges_norm = compute_bounds(ancestor, mesh, k_opt)
    bin_edges =  xmin .+ (xmax - xmin) * bin_edges_norm
    N = bin_irregular(x, bin_edges, right)
    if right
        H = Histogram(bin_edges, N, :right, true)
    else
        H = Histogram(bin_edges, N, :left, true)
    end
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    if rule == "bayes"
        H.weights = (H.weights .+ a*p0) ./ ((n + a)*(bin_edges[2:end] - bin_edges[1:end-1]))
    else
        H.weights = H.weights ./ (n * (bin_edges[2:end] - bin_edges[1:end-1]) )
    end
    H.isdensity = true
    return H, criterion_opt
end