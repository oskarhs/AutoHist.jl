"""
    histogram_irregular(x::AbstractVector{<:Real}; rule::Symbol=:bayes, grid::Symbol=:regular, closed::Symbol=:right, greedy::Bool=true, maxbins::Int=-1, support::Tuple{Real,Real}=(-Inf,Inf), use_min_length::Bool=false, logprior::Function=k->0.0, a::Real=1.0)

Create an irregular histogram based on optimization of a criterion based on Bayesian probability, penalized likelihood or LOOCV.
Returns an AutomaticHistogram object with the optimal partition corresponding to the supplied rule.

# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Defaults to the Bayesian method of Simensen et al. (2025).
- `grid`: Symbol indicating how the finest possible mesh should be constructed. Options are `:data`, which uses each unique data point as a grid point, `:regular` (default) which constructs a fine regular grid, and `:quantile` which constructs the grid based on the sample quantiles.
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `greedy`: Boolean indicating whether or not the greedy binning strategy of Rozenholc et al. (2006) should be used prior to running the dynamical programming algorithm. Defaults to `true`. The algorithm can be quite slow for large datasets when this keyword is set to `false`.
- `maxbins`: The maximal number of bins to be considered by the optimization criterion, only used if grid is set to "regular" or "quantile". Defaults to `maxbins=min(4*n/log(n)^2, 1000)`. If the specified argument is not a positive integer, the default value is used.
- `support`: Tuple specifying the the support of the histogram estimate. If the first element is -Inf, then `minimum(x)` is taken as the leftmost cutpoint. Likewise, if the second elemen is `Inf`, then the rightmost cutpoint is `maximum(x)`. Default value is `(-Inf, Inf)`, which estimates the support of the data.
- `use_min_length`: Boolean indicating whether or not to impose a restriction on the minimum bin length of the histogram. If set to true, the smallest allowed bin length is set to `(maximum(x)-minimum(x))/n*log(n)^(1.5)`.
- `logprior`: Unnormalized logprior distribution for the number k of bins. Defaults to a uniform prior. Only used in when `rule` keyword is set to `:bayes`.
- `a`: Dirichlet concentration parameter in the Bayesian irregular histogram model. Set to the default value (5.0) if the supplied value is not a positive real number. Only used when `rule` is set to `:bayes`.

# Returns
- `h`: AutomaticHistogram object with weights corresponding to densities, e.g. `:isdensity` is set to true.

# Examples
```
julia> x = [0.037, 0.208, 0.189, 0.656, 0.45, 0.846, 0.986, 0.751, 0.249, 0.447]
julia> h1 = histogram_irregular(x)
julia> h2 = histogram_irregular(x; grid=:quantile, support=(0.0, 1.0), logprior=k->-log(k), a=sqrt(10))
```
...
"""
function histogram_irregular(x::AbstractVector{<:Real}; rule::Symbol=:bayes, grid::Symbol=:regular, 
                            closed::Symbol=:right, greedy::Bool=true, maxbins::Int=-1, support::Tuple{Real,Real}=(-Inf,Inf),
                            use_min_length::Bool=false, logprior::Function=k->0.0, a::Real=5.0)
    if !(rule in [:pena, :penb, :penr, :bayes, :klcv, :l2cv, :nml])
        throw(ArgumentError("The supplied rule, :$rule, is not supported. The rule kwarg must be one of :pena, :penb, :penr, :bayes, :klcv, :l2cv or :nml"))
    end
    if rule == :bayes && a ≤ 0.0
        throw(DomainError("Supplied value of a must be positive."))
    end
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end

    if !(grid in [:data, :regular, :quantile])
        throw(ArgumentError("The supplied grid, :$grid, is not supported. The grid kwarg must be one of :data, :regular or :quantile"))
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
    y = @. (x - xmin) / (xmax - xmin)
    n = length(x)

    if grid == :data
        maxbins = n-1
    elseif maxbins ≤ 0
        maxbins = min(n, ceil(Int, 4.0*n/log(n)^2))
    end

    # Calculate gridpoints (left-open grid, breaks at data points)
    # closed == :right means to include observation in the right endpoint, i.e. right-closed
    finestgrid = Array{Float64}(undef, maxbins+1)
    N_cum = Array{Float64}(undef, length(finestgrid)) # cumulative cell counts
    N_cum[1] = 0.0
    if grid == :data
        sort!(y)
        finestgrid[1] = -eps()
        finestgrid[2:end-1] = y[2:n-1]
        finestgrid[end] = 1.0 + eps()
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    elseif grid == :regular
        N_cum[2:end] = cumsum(bin_regular(y, 0.0, 1.0, maxbins, closed == :right))
        finestgrid[1:end] = LinRange(0.0, 1.0, maxbins+1)
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0+eps()
    elseif grid == :quantile
        sort!(y)
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0 + eps()
        finestgrid[2:end-1] = quantile(y, LinRange(1.0/maxbins, 1.0-1.0/maxbins, maxbins-1); sorted=true)
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
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

    # Set objective for the dynamic programming part according to the suppplied rule
    if rule in [:pena, :penb, :nml]
        phi = let N_cum = N_cum, mesh = mesh
            f(i,j) = phi_penB(i, j, N_cum, mesh)
        end
    elseif rule == :bayes
        phi = let N_cum = N_cum, mesh = mesh, a = a
            f(i,j) = phi_bayes(i, j, N_cum, mesh, a)
        end
    elseif rule == :penr
        phi = let N_cum = N_cum, mesh = mesh, n = n
            f(i,j) = phi_penR(i, j, N_cum, mesh, n)
        end
    elseif rule == :klcv
        minlength = 0.0
        if use_min_length
            minlength = log(n)^(1.5)/n
        end
        phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
            f(i,j) = phi_KLCV(i, j, N_cum, mesh, n; minlength=minlength)
        end
    elseif rule == :l2cv
        minlength = 0.0
        if use_min_length
            minlength = log(n)^(1.5)/n
        end
        phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
            f(i,j) = phi_L2CV(i, j, N_cum, mesh, n; minlength=minlength)
        end
    end

    optimal, ancestor = dynamic_algorithm(phi, k_max)
    psi = Array{Float64}(undef, k_max)
    if rule == :penb   
        @simd for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
        end
    elseif rule == :bayes
        @simd for k = 1:k_max
            @inbounds psi[k] = logprior(k) - logabsbinomial(maxbins-1, k-1)[1] + loggamma(a) - loggamma(a + n)
        end
    elseif rule == :penr
        @simd for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
        end
    elseif rule == :pena
        @simd for k = 1:k_max
            @inbounds psi[k] = -logabsbinomial(maxbins-1, k-1)[1] - k - 2.0*log(k) -
                    2.0 * sqrt(1.0*0.5*(k-1)*(logabsbinomial(maxbins-1, k-1)[1] + 1.0*log(k)))
        end
    elseif rule == :nml
        @simd for k = 1:k_max
            @inbounds psi[k] =  -( 0.5*k*log(0.5*n) - loggamma(0.5*k) +
            1.0/sqrt(n) * sqrt(2.0)*k/3.0 * exp(loggamma(0.5*k) - loggamma(0.5*k-0.5)) +
            1.0/n * ((3.0 + k*(k-2.0)*(2.0*k+1.0))/36.0 - k^2/9.0*exp(2.0*loggamma(0.5*k) - 2.0*loggamma(0.5*k-0.5)))
            )
            psi[k] = psi[k] - logabsbinomial(maxbins-1, k-1)[1]
        end
    end # NB! no penalties for klcv and l2cv
    k_opt = argmax(optimal + psi)

    bin_edges_norm = compute_bounds(ancestor, mesh, k_opt)
    bin_edges = @. xmin + (xmax - xmin) * bin_edges_norm

    N = convert.(Int64, bin_irregular(x, bin_edges, closed == :right))
    a_opt = 0.0
    if rule == :bayes
        a_opt = a
    end
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    dens = (N .+ a_opt*p0) ./ ((n + a_opt)*(bin_edges[2:end] - bin_edges[1:end-1]))
    h = AutomaticHistogram(bin_edges, dens, N, :irregular, closed, ifelse(a_opt > 0.0, a_opt, NaN))
    return h
    
    #H = Histogram(bin_edges, dens, closed, true)
    
#=     N = bin_irregular(x, bin_edges, closed == :right)
    H = Histogram(bin_edges, dens, closed, true)
    p0 = bin_edges_norm[2:end] - bin_edges_norm[1:end-1]
    if rule == :bayes
        H.weights = (H.weights .+ a*p0) ./ ((n + a)*(bin_edges[2:end] - bin_edges[1:end-1]))
    else
        H.weights = H.weights ./ (n * (bin_edges[2:end] - bin_edges[1:end-1]) )
    end
    H.isdensity = true =#
    #return H
end