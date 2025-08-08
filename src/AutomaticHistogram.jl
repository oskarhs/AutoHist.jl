"""
    AutomaticHistogram

A type for representing a histogram where the histogram partition has been chosen automatically based on the sample. Can be fitted to data using the [`fit`](@ref) method.

# Fields
- `breaks`: AbstractVector consisting of the cut points in the chosen partition.
- `density`: Estimated density in each bin.
- `counts`: The bin counts for the partition corresponding to `breaks`.
- `type`: Symbol indicating whether the histogram was fit using an irregular procedure (`type==:irregular`) or a regular one (`type==:regular`).
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `a`: Value of the Dirichlet concentration parameter corresponding to the chosen partition. Only of relevance if a Bayesian method was used to fit the histogram, and is otherwise set to `NaN`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x)
AutomaticHistogram
breaks: [0.0001220703125, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9202995853147645, 1.0]
density: [0.006626835974128547, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.733509595364622, 3.545742066060377]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1090, 1414]
type: irregular
closed: right
a: 5.0
```
"""
struct AutomaticHistogram 
    breaks::AbstractVector{Float64}
    density::AbstractVector{Float64}
    counts::AbstractVector{Int}
    type::Symbol
    closed::Symbol
    a::Float64
end

AutomaticHistogram(breaks::AbstractVector{Float64}, density::AbstractVector{Float64}, counts::AbstractVector{Int}, type::Symbol, closed::Symbol) = AutomaticHistogram(breaks, density, counts, type, closed, NaN)

"""
    length(h::AutomaticHistogram)

Returns the number of bins of `h`.
"""
Base.length(h::AutomaticHistogram) = length(h.density)

#= """
    ==(h1::AutomaticHistogram, h2::AutomaticHistogram)

Return `true` if all fields of `h1` and `h2` are equal. Otherwise, return `false`.
""" =#
Base.:(==)(h1::AutomaticHistogram, h2::AutomaticHistogram) = all(isequal(getfield(h1, f), getfield(h2, f)) for f in fieldnames(AutomaticHistogram))

#= """
    ≈(h1::AutomaticHistogram, h2::AutomaticHistogram)

Return `true` if all numeric fields of `h1` and `h2` are approximately equal, and all symbol fields are exactly equal. Otherwise, return `false`.
""" =#
function Base.:(≈)(h1::AutomaticHistogram, h2::AutomaticHistogram)
    if ( !isnan(getfield(h1, :a)) ) && isnan(getfield(h2, :a))
        return false
    elseif ( !isnan(getfield(h2, :a)) ) && isnan(getfield(h1, :a))
        return false
    elseif length(h1) != length(h2)
        return false
    else 
        return all(isapprox(getfield(h1, f), getfield(h2, f)) for f in (:breaks, :density, :counts)) && all(isequal(getfield(h1, f), getfield(h2, f)) for f in (:type, :closed))
    end
end

function Base.show(io::IO, h::AutomaticHistogram)
    println(io, typeof(h))
    println(io, "breaks: ", h.breaks)
    println(io, "density: ", h.density)
    println(io, "counts: ", h.counts)
    println(io, "type: ", h.type)
    println(io, "closed: ", h.closed)
    println(io, "a: ", h.a)
end


"""
    fit(AutomaticHistogram, x::AbstractVector{<:Real}, rule::AbstractRule=RIH(); support::Tuple{Real,Real}=(-Inf,Inf), closed::Symbol=:right)

Fit a histogram to a one-dimensional vector `x` with an automatic and data-based selection of the histogram partition.

# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Default value is `rule=RIH()`, the random irregular histogram.
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `support`: Tuple specifying the the support of the histogram estimate. If the first element is `-Inf`, then `minimum(x)` is taken as the leftmost cutpoint. Likewise, if the second element is `Inf`, then the rightmost cutpoint is `maximum(x)`. Default value is `(-Inf, Inf)`, which estimates the support of the data.

# Returns
- `h`: An object of type [`AutomaticHistogram`](@ref), corresponding to the fitted histogram.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = (1.0 .- (1.0 .- LinRange(0.0, 1.0, 5000)) .^(1/3)).^(1/3);

julia> fit(AutomaticHistogram, x) == fit(AutomaticHistogram, x, RIH())
true

julia> h = fit(AutomaticHistogram, x, Wand(scalest=:stdev, level=4))
AutomaticHistogram
breaks: LinRange{Float64}(0.0, 1.0, 27)
density: [0.0052, 0.0312, 0.0884, 0.1612, 0.2652, 0.4004, 0.5408, 0.7176, 0.8944, 1.0868  …  2.0072, 1.9656, 1.8616, 1.69, 1.4508, 1.1596, 0.8372, 0.5044, 0.2184, 0.0364]
counts: [1, 6, 17, 31, 51, 77, 104, 138, 172, 209  …  386, 378, 358, 325, 279, 223, 161, 97, 42, 7]
type: regular
closed: right
a: NaN
```
"""
function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}, rule::AbstractRule=RIH(); support::Tuple{Real,Real}=(-Inf,Inf), closed::Symbol=:right)
    if !(closed in [:right, :left]) # if supplied symbol is nonsense, just use default
        throw(ArgumentError("The supplied value of the closed keyword, :$closed, is invalid. Valid values are :left or :right."))
    end
    xmin, xmax = extrema(x)

    if support[1] > -Inf       # estimate lower bound of support if unknown,
        if xmin > support[1]
            xmin = support[1]  # use known lower bound
        else 
            throw(DomainError("The supplied lower bound is greater than the smallest value of the sample."))
        end
    end
    if support[2] < Inf
        if xmax < support[2]
            xmax = support[2]  # use known upper bound
        else 
            throw(DomainError("The supplied upper bound is smaller than the smallest value of the sample."))
        end
    end
    return fit_autohist(x, rule, xmin, xmax, closed)
end

"""
    convert(Histogram, h::AutomaticHistogram)

Convert an `h` to a StatsBase.Histogram, normalized to be a probability density.
"""
Base.convert(::Type{Histogram}, h::AutomaticHistogram) = Histogram(h.breaks, h.density, h.closed, true)

"""
    loglikelihood(h::AutomaticHistogram)

Compute the log-likelihood (up to proportionality) of an `h`.

The value of the log-likelihood is
    ``\\sum_j N_j \\log (d_j)``
where ``N_j``, ``d_j`` are the bin counts and estimated densities for bin j.
"""
function loglikelihood(h::AutomaticHistogram)
    llik = 0.0
    @inbounds for j in eachindex(h.density)
        if h.counts[j] > 0
            llik += h.counts[j]*log(h.density[j])
        end
    end
    return llik
end

"""
    logmarginallikelihood(h::AutomaticHistogram, a::Real)
    logmarginallikelihood(h::AutomaticHistogram)

Compute the log-marginal likelihood (up to proportionality) of `h` when the value of the Dirichlet concentration parameter equals `a`. This can be automatically inferred if the histogram was fitted with the `rule` argument set to [`RIH`](@ref) or [`RRH`](@ref), and does not have to be explicitly passed as an argument in this case.

Assumes that the Dirichlet prior is centered on the uniform distribution, so that ``a_j = a/k`` for a scalar ``a>0`` and all ``j``.
The value of the log-marginal likelihood is
``\\sum_j \\{ \\log \\Gamma (a_j + N_j) - \\log \\Gamma (a_j) - N_j\\log |\\mathcal{I}_j| \\} - \\log \\Gamma (a+n) + \\log \\Gamma (a)``
, where ``N_j`` is the bin count for bin ``j`` .
"""
function logmarginallikelihood(h::AutomaticHistogram, a::Real)
    k = length(h.counts)
    logmarglik = loggamma(a) - loggamma(a+sum(h.counts))
    @inbounds for j in eachindex(h.counts)
        if h.counts[j] > 0
            len_bin = h.breaks[j+1] - h.breaks[j]
            logmarglik += loggamma(h.counts[j] + a/k) - loggamma(a/k) - h.counts[j]*len_bin
        end
    end
    return logmarglik
end

function logmarginallikelihood(h::AutomaticHistogram)
    if isnan(h.a)
        throw(ArgumentError("If h was not fit with rule=:bayes, the value of a must be explicitly specified to compute the log-marginal likelihood."))
    end
    return logmarginallikelihood(h, h.a)
end

"""
    minimum(h::AutomaticHistogram)

Return the minimum of the support of `h`.
"""
Base.minimum(h::AutomaticHistogram) = h.breaks[1]

"""
    maximum(h::AutomaticHistogram)

Return the maximum of the support of `h`.
"""
Base.maximum(h::AutomaticHistogram) = h.breaks[end]

"""
    extrema(h::AutomaticHistogram)

Return the minimum and the maximum of the support of `h` as a 2-tuple.
"""
Base.extrema(h::AutomaticHistogram) = (h.breaks[1], h.breaks[end])

"""
    insupport(h::AutomaticHistogram, x::Real)

Return `true` if `x` is in the support of `h`, and `false` otherwise.
"""
insupport(h::AutomaticHistogram, x::Real) = (h.breaks[1] <= x <= h.breaks[end])

"""
    pdf(h::AutomaticHistogram, x::Real)

Evaluate the probability density function of `h` at `x`.
"""
function pdf(h::AutomaticHistogram, x::Real)
    val = 0.0
    if insupport(h, x)
        k = length(h)
        if h.type == :irregular
            if h.closed == :right
                idx = max(1, searchsortedfirst(h.breaks, x) - 1)
                @inbounds val = h.density[idx]
            else 
                idx = min(k, searchsortedlast(h.breaks, x))
                @inbounds val = h.density[idx]
            end
        else
            xmin, xmax = extrema(h)
            edges_inc = k/(xmax-xmin)
            if h.closed == :right
                idx = max(0, floor(Int, (x-xmin)*edges_inc-10.0*edges_inc*eps())) + 1
                @inbounds val = h.density[idx]
            else 
                idx = min(k-1, floor(Int, (x-xmin)*edges_inc+10.0*edges_inc*eps())) + 1
                @inbounds val = h.density[idx]
            end
        end
    end
    return val
end

# Allows for pdf.(h, x) where x is a vector.
function Base.Broadcast.broadcasted(::typeof(pdf), h::AutomaticHistogram, x::AbstractVector)
    vals = Vector{Float64}(undef, length(x))
    @inbounds for i in eachindex(x)
        vals[i] = pdf(h, x[i])
    end
    return vals
end

"""
    cdf(h::AutomaticHistogram, x::Real)

Evaluate the cumulative distribution function of `h` at `x`.
"""
function cdf(h::AutomaticHistogram, x::Real)
    val = 0.0
    xmin, xmax = extrema(h)
    if xmin < x < xmax
        k = length(h)
        if h.type == :irregular
            if h.closed == :right
                idx = max(1, searchsortedfirst(h.breaks, x) - 1)
                @inbounds val = h.density[idx] * (x - h.breaks[idx]) + sum(diff(h.breaks[1:idx]) .* h.density[1:idx-1])
            else 
                idx = min(k, searchsortedlast(h.breaks, x))
                @inbounds val = h.density[idx] * (x - h.breaks[idx]) + sum(diff(h.breaks[1:idx]) .* h.density[1:idx-1])
            end
        else
            edges_inc = k/(xmax-xmin)
            if h.closed == :right
                idx = max(0, floor(Int, (x-xmin)*edges_inc-10.0*edges_inc*eps())) + 1
                @inbounds val = h.density[idx] * (x - h.breaks[idx]) + sum(diff(h.breaks[1:idx]) .* h.density[1:idx-1])
            else 
                idx = min(k-1, floor(Int, (x-xmin)*edges_inc+10.0*edges_inc*eps())) + 1
                @inbounds val = h.density[idx] * (x - h.breaks[idx]) + sum(diff(h.breaks[1:idx]) .* h.density[1:idx-1])
            end
        end
    elseif x ≥ xmax
        val = 1.0
    end
    return val
end

# Allows for cdf.(h, x) where x is a vector.
function Base.Broadcast.broadcasted(::typeof(cdf), h::AutomaticHistogram, x::AbstractVector)
    vals = Vector{Float64}(undef, length(x))
    @inbounds for i in eachindex(x)
        vals[i] = cdf(h, x[i])
    end
    return vals
end

"""
    peaks(h::AutomaticHistogram)

Return the location of the modes/peaks of `h` as a Vector, sorted in increasing order.

Formally, the modes/peaks of the histogram `h` are defined as the midpoints of an interval ``\\mathcal{J}``, where the density of `h` is constant on ``\\mathcal{J}``, and the density of `h` is strictly smaller than this value in the histogram bins adjacent to ``\\mathcal{J}``. Note that according this definition, ``\\mathcal{J}`` is in general a nonempty union of intervals in the histogram partition.
"""
function peaks(h::AutomaticHistogram)
    breaks = h.breaks
    dens = h.density
    k = length(dens)

    # Make new breaks and density vectors where equal-density bins are 'concatenated'
    dens1 = [dens[1]]
    breaks1 = [breaks[1]]
    for j in 1:(k-1)
        if !(isapprox(dens[j], dens[j+1]; atol=1e-10))
            push!(dens1, dens[j+1])
            push!(breaks1, breaks[j+1])
        end
    end
    push!(breaks1, breaks[end])

    k1 = length(dens1)
    # Now identify modes of the histogram density
    hist_modes = Float64[]
    if k1 == 1
        push!(hist_modes, 0.5*(breaks1[1]+breaks1[2]))
    else
        if dens1[1] > dens1[2]
            push!(hist_modes, 0.5*(breaks1[1]+breaks1[2]))
        end
        for j in 2:(k1-1)
            if dens1[j] > dens1[j-1] && dens1[j] > dens1[j+1]
                push!(hist_modes, 0.5*(breaks1[j]+breaks1[j+1]))
            end
        end
        if dens1[end] > dens1[end-1]
            push!(hist_modes, 0.5*(breaks1[end-1]+breaks1[end]))
        end
    end
    return hist_modes
end


"""
    distance(h1::AutomaticHistogram, h2::AutomaticHistogram, dist::Symbol=:iae; p::Real=1.0)

Compute a statistical distance between two histogram probability densities.

# Arguments
- `h1`, `h2`: The two histograms for which the distance should be computed
- `dist`: The name of the distance to compute. Valid options are `:iae` (default), `:ise`, `:hellinger`, `:sup`, `:kl`, `:lp`. For the ``L_p``-metric, a given power `p` can be specified as a keyword argument. 

# Keyword arguments
- `p`: Power of the ``L_p``-metric, which should be a number in the interval ``[1, \\infty]``. Ignored if `dist != :lp`. Defaults to `p=1.0`.
"""
function distance(h1::AutomaticHistogram, h2::AutomaticHistogram, dist::Symbol=:iae; p::Real=1.0)
    if !(dist in [:iae, :ise, :lp, :hell, :sup, :kl])
        throw(ArgumentError("The supplied distance, :$dist, is not supported."))
    end
    if dist == :iae
        return lp_distance(h1, h2, 1.0)
    elseif dist == :ise
        return lp_distance(h1, h2, 2.0)^2
    elseif dist == :sup
        return supremum_distance(h1, h2)
    elseif dist == :lp
        if p ≤ 1.0
            throw(DomainError("Power parameter p must be a number in the interval [1, p]."))
        elseif p < Inf
            return lp_distance(h1, h2, p)
        else
            return supremum_distance(h1, h2)
        end
    elseif dist == :kl
        return kl_divergence(h1, h2)
    else
        return hellinger_distance(h1, h2)
    end
end