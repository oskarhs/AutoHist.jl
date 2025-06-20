"""
    AutomaticHistogram

A type for representing a histogram where the histogram partition has been chosen automatically based on the sample. Can be fitted to data using the [`fit`](@ref), [`histogram_irregular`](@ref) or [`histogram_regular`](@ref) methods.

# Fields
- `breaks`: AbstractVector consisting of the cut points in the chosen partition.
- `density`: Estimated density in each bin.
- `counts`: The bin counts for the partition corresponding to `breaks`.
- `type`: Symbol indicating whether the histogram was fit using an irregular procedure (`type==:irregular`) or a regular one (`type==:regular`).
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `a`: Value of the Dirichlet concentration parameter corresponding to the chosen partition. Only of relevance if a Bayesian method was used to fit the histogram, and is otherwise set to `NaN`.

# Examples
```jldoctest 
julia> using AutoHist

julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x)
AutomaticHistogram
breaks: [0.00012207031249977798, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9239223314368207, 1.0000000000000002]
density: [0.00662683597412854, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.7513848134529346, 3.564842182949155]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1147, 1357]
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
    fit(AutomaticHistogram, x::AbstractVector{x<:Real}; rule=:default, type=:irregular, kwargs...)

Fit a histogram to a one-dimensional vector `x` with an automatic and data-based selection of the histogram partition.

...
# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Default value is `rule=:default` which uses the default rule for the regular or irregular histogram procedure depending on the value of `type`.
- `type`: Symbol indicating whether the fitted method is a regular and irregular one. The rules `:bayes`, `:l2cv`, `:klcv` and `:nml` are implemented for both regular and irregular histograms, and this keyword specifies whether the regular or irregular version should be used. For other rules, this function infers the type automatically from the `rule` keyword, and misspecifying the rule in this case has not effect. Possible values are `:irregular` (default) and `:regular`.
- `kwargs`: Additional keyword arguments passed to [`histogram_regular`](@ref) or [`histogram_irregular`](@ref) depending on the specified or inferred type.

# Returns
- `h`: An object of type [`AutomaticHistogram`](@ref), corresponding to the fitted histogram.

# Examples
```julia
julia> x = randn(10^3)
julia> h1 = fit(AutomaticHistogram, x)                                      # fits an irregular histogram
julia> h2 = fit(AutomaticHistogram, x; rule=:wand, scalest=:stdev, level=4) # fits a regular histogram
```
"""
function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}; rule::Symbol=:default, type::Symbol=:irregular, kwargs...)
    if rule == :default
        if type == :irregular
            rule = :bayes
        elseif type == :regular
            rule = :knuth
        else 
            throw(ArgumentError("The supplied type, :$type, is not supported and the supplied rule=:default could as such not be inferred. Valid types are :irregular and :regular."))
        end
    elseif rule in [:aic, :bic, :br, :mdl, :sturges, :fd, :scott, :wand]
        type = :regular
    elseif rule in [:pena, :penb, :penr]
        type = :irregular
    elseif rule in [:bayes, :l2cv, :klcv, :nml]
        if !(type in [:regular, :irregular]) # if type is not specified correctly for any of these methods, throw an ArgumentError as type is ambiguous in this case.
            throw(ArgumentError("Unable to infer type automatically from the supplied rule and the supplied type is not supported. Valid types are :irregular and :regular."))
        end
    else
        throw(ArgumentError("The supplied rule, :$rule, is not supported. The rule kwarg must be one of :aic, :bic, :br, :bayes, :mdl, :klcv, :nml, :l2cv, :sturges, :fd, :scott, :wand, :pena, :penb or :penr."))
    end

    # Fit the histogram here
    if type == :irregular
        h = histogram_irregular(x; rule=rule, kwargs...)
    elseif type == :regular
        h = histogram_regular(x; rule=rule, kwargs...)
    end
    return h
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
    ∑ⱼ Nⱼ log (dⱼ),
where Nⱼ, dⱼ are the bin counts and estimated densities for bin j.
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

Compute the log-marginal likelihood (up to proportionality) of `h` when the value of the Dirichlet concentration parameter equals `a`. This can be automatically inferred if the histogram was fitted with `rule=:bayes`, and does not have to be explicitly passed as an argument in this case.

Assumes that the Dirichlet prior is centered on the uniform distribution, so that aⱼ = a/k for a scalar a>0 and all j.
The value of the log-marginal likelihood is
    ∑ⱼ {log Γ(Nⱼ+aⱼ) - log Γ(aⱼ) - Nⱼlog(dⱼ)} - log Γ(a+n) + log Γ(a),
where where Nⱼ is the bin count for bin j.
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
    vals = Array{Float64}(undef, length(x))
    @inbounds for i in eachindex(x)
        vals[i] = pdf(h, x[i])
    end
    return vals
end

"""
    peaks(h::AutomaticHistogram)

Return the location of the modes/peaks of `h` as a Vector, sorted in increasing order.

Formally, the modes/peaks of the histogram `h` are defined as the midpoints of an interval ``\\mathcal{J}``, where the density of `h` is constant on ``\\mathcal{J}``, and the density of `h` is strictly smaller than this value in the histogram bins adjacent to ``\\mathcal{J}``. Note that according this definition, ``\\mathcal{J}`` is in general a nonempty union of intervals in the histogram partition.
"""
function peaks(h::AutomaticHistogram)
    # implementation here (see peak_id_loss from loss_functions.jl)
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