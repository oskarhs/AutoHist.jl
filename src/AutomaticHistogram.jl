"""
    AutomaticHistogram

A type for representing a histogram where the histogram partition has been chosen automatically based on the sample.

...
# Fields
- `breaks`: AbstractVector consisting of the cut points in the chosen partition.
- `density`: Estimated density in each bin.
- `counts`: The bin counts for the partition corresponding to `breaks`.
- `type`: Symbol indicating whether the histogram was fit using an irregular procedure (`type==:irregular`) or a regular one (`type==:regular`).
- `closed`: Symbol indicating whether the drawn intervals should be right-inclusive or not. Possible values are `:right` (default) and `:left`.
- `a`: Value of the Dirichlet concentration parameter corresponding to the chosen partition. Only of relevance if a Bayesian method was used to fit the histogram, and is otherwise set to `NaN`.

# Examples
"""
struct AutomaticHistogram 
    breaks::AbstractVector{Float64}
    density::AbstractVector{Float64}
    counts::AbstractVector{Int}
    type::Symbol
    closed::Symbol
    a::Float64
end

AutomaticHistogram(breaks::AbstractVector{Float64}, density::AbstractVector{Float64}, counts::AbstractVector{Int}, type::Symbol, closed::Symbol) = AutomaticHistogram(breaks, density, counts, type, NaN)

Base.:(==)(h1::AutomaticHistogram, h2::AutomaticHistogram) = all(getfield(h1, f) == getfield(h2, f) for f in fieldnames(AutomaticHistogram))

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
    fit(AutomaticHistogram, x::AbstractVector{x<:Real})

Fit a histogram to a one-dimensional vector x with an automatic and data-based selection of the histogram partition.

...
# Arguments
- `x`: 1D vector of data for which a histogram is to be constructed.

# Keyword arguments
- `rule`: The criterion used to determine the optimal number of bins. Defaults to the method Bayesian method of Simensen et al. (2025).
- `type`: Symbol indicating whether the fitted method is a regular and irregular one. The rules `:bayes`, `:l2cv`, `:klcv` and `:nml` are implemented for both regular and irregular histograms, and this keyword specifies whether the regular or irregular version should be used. For other rules, this function infers the type automatically from the `rule` keyword, and misspecifying the rule in this case has not effect. Possible values are `:irregular` (default) and `regular`.
- `kwargs`: Additional keyword arguments passed to `histogram_regular` or `histogram_irregular` depending on the specified or inferred type.

# Returns
- `h`: An object of type `AutomaticHistogram`, corresponding to the fitted histogram.
"""
function fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}; rule=:bayes, type=:irregular, kwargs...)
    if rule in [:aic, :bic, :br, :mdl, :sturges, :fd, :scott, :wand]
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

Convert an `AutomaticHistogram` to a StatsBase.Histogram, normalized to be a density.
"""
Base.convert(::Type{Histogram}, h::AutomaticHistogram) = Histogram(h.breaks, h.density, h.closed, true)

"""
    loglikelihood(h::AutomaticHistogram)

Compute the log-likelihood (up to proportionality) of an `AutomaticHistogram`.

...
The value of the log-likelihood is
    ∑ⱼ Nⱼ log (dⱼ),
where Nⱼ, dⱼ are the bin counts and estimated densities for bin j.
"""
function loglikelihood(h::AutomaticHistogram)
    llik = -n*log(n)
    for j in eachindex(h.density)
        if h.counts[j] > 0
            llik += h.counts[j]*log(h.density[j])
        end
    end
    return llik
end

"""
    logmarginallikelihood(h::AutomaticHistogram, a::Real)
    logmarginallikelihood(h::AutomaticHistogram)

Compute the log-marginal likelihood (up to proportionality) of an `AutomaticHistogram` when the value of the Dirichlet concentration parameter equals a. This can be automatically inferred if the histogram was fitted with `rule=:bayes`, and does not have to be explicitly passed as an argument in this case.

...
Assumes that the Dirichlet prior is centered on the uniform distribution, so that aⱼ = a/k for a scalar a>0 and all j.
The value of the log-marginal likelihood is
    ∑ⱼ {log Γ(Nⱼ+aⱼ) - log Γ(aⱼ) - Nⱼlog(dⱼ)} - log Γ(a+n) + log Γ(a),
where where Nⱼ is the bin count for bin j.
"""
function logmarginallikelihood(h::AutomaticHistogram, a::Real)
    logmarglik = loggamma(a) - loggamma(a+sum(h.counts))
    for j in eachindex(h.counts)
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
    return loglikelihood(h, h.a)
end