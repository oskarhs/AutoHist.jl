using Integrals, Distributions, AutoHist
import MergeSorted: mergesorted

#include("test_distributions.jl")

"""
    discontinuities(d::ContinuousUnivariateDistribution)

Return the points at which the probability density function of `d` is discontinuous as a vector in increasing order. This method is intended to be overridden. The fallback value is an empty vector, corresponding to no discontinuities.
"""
function discontinuities(d::ContinuousUnivariateDistribution)
    return Float64[]
end

"""
    distance(h::AutomaticHistogram, d::ContinuousUnivariateDistribution, name::Symbol, alg::SciMLBase.AbstractIntegralAlgorithm)

Numerically approximate a given statistical distance between a histogram density estimate and a distribution.

This method uses Integrals.jl to approximate common loss functions used in density estimation, with special care taken to avoid issues encountered at points of discontinuity of the estimator or the true density.
This is accomplished by performing the numeric quadrature piecewise over intervals where both the histogram estimate and the true density are continuous.

# Arguments
- `h`: The histogram density estimate.
- `d`: The distribution we are trying to estimate.
- `name`: The distance we would like to compute. Valid options are `:lp`, `:hellinger`, `:ise` and `iae`. For `name=:lp`, the lp-metric between the two densities is computed for some `p >= 1` controlled via a keyword argument (`p = 2` by default).
  `name=:hellinger` corresponds to the Hellinger metric. The options `:ise` (integrated squared error) and `:iae` (integrated absolute error) are aliases for `:lp` with `p = 2` and `p = 1`, respectively.
- `alg`: The numeric quadrature algorithm to use. Supports any algorithm that can be passed to `Integrals.solve()`, see [the Integrals.jl documentation.](https://docs.sciml.ai/Integrals/stable/solvers/IntegralSolvers/), 

# Keyword arguments
- `p`: The value of `p` in the lp metric. Ignored if `name != :lp`. Defaults to `p = 2.0`.

# Return
`dist`: The approximate value of the distance.

!!! note
    To get more accurate approximations from this method, the `discontinuities` method should be implemented for `d` prior to usage, as
    otherwise, the method will treat `d` as continuous over its support.

*** Write a jldoctest here later ***

"""
function distance(h::AutomaticHistogram, d::ContinuousUnivariateDistribution, name::Symbol, alg::SciMLBase.AbstractIntegralAlgorithm; p = 2.0)
    if !(name in [:lp, :hellinger, :ise, :iae]) # write a test for this
        throw(ArgumentError("The supplied distance, :$name, is not supported. Valid options are :lp, :hellinger, :ise or :iae."))
    end
    if name == :ise
        name = :lp
        p = 2.0
    elseif name == :iae
        name = :lp
        p = 1.0
    end

    if name == :lp
        if p < 1.0 # write a test for this
            throw(DomainError("p must be >= 1.0."))
        end
        dist = _lp_loss(d, h.breaks, h.density, alg, p)
    else
        dist = _hellinger_loss(d, h.breaks, h.density, alg)
    end
    return dist
end

function _hellinger_loss(d::ContinuousUnivariateDistribution, breaks_hist::AbstractVector{<:Real}, dens_hist::AbstractVector{<:Real}, alg::SciMLBase.AbstractIntegralAlgorithm)
    disc = mergesorted(discontinuities(d), breaks_hist) # union of discontinuities of f₀ and \hat{f} in ascending order
    # Now perform numerical quadrature piecewise over each interval where √f₀-√\hat{f} is continuous
    m = length(disc) - 1
    hell = cdf(d, disc[1]-10.0*eps())
    for j = 1:m
        bin_ind = searchsortedfirst(breaks_hist, disc[j]+10.0*eps()) - 1
        if bin_ind == 0 || bin_ind == (length(dens_hist)+1)
            dens = 0.0
        else
            dens = dens_hist[bin_ind]
        end
        ip = IntegralProblem((t,p) -> (sqrt(pdf(d, t)) - sqrt(dens))^2, disc[j]+10.0*eps(), disc[j+1]-10.0*eps())
        hell += solve(ip, alg).u
    end
    hell += 1.0 - cdf(d, disc[end]+10*eps())
    return sqrt(hell)
end

function _lp_loss(d::ContinuousUnivariateDistribution, breaks_hist::AbstractVector{<:Real}, dens_hist::AbstractVector{<:Real}, alg::SciMLBase.AbstractIntegralAlgorithm, p::Real)
    disc = mergesorted(discontinuities(d), breaks_hist)
    supp = support(d)
    # Now perform numerical quadrature piecewise over each interval where f₀-\hat{f} is continuous
    m = length(disc) - 1
    ip = IntegralProblem((t,p) -> pdf(d, t)^p, supp.lb, disc[1]-10.0*eps())
    lp = solve(ip, alg).u
    for j = 1:m
        bin_ind = searchsortedfirst(breaks_hist, disc[j]+10.0*eps()) - 1
        if bin_ind == 0 || bin_ind == (length(dens_hist)+1)
            dens = 0.0
        else
            dens = dens_hist[bin_ind]
        end
        ip = IntegralProblem((t,p) -> (pdf(d, t) - dens)^p, disc[j]+10.0*eps(), disc[j+1]-10.0*eps())
        lp += solve(ip, alg).u
    end
    ip = IntegralProblem((t,p) -> pdf(d, t)^p, disc[end]+10.0*eps(), supp.ub)
    lp += solve(ip, alg).u
    return lp^(1.0/p)
end


function peak_id_loss(d::ContinuousUnivariateDistribution, breaks_hist::AbstractVector{<:Real}, dens_hist::AbstractVector{<:Real};
                      verbose::Bool=false)
    # Identify locations of modes in the histogram density
    # First combine bins that have the same density (only needed for regular histograms)
    true_modes = modes(d)
    delta = pid_tolerance(d)
    xmin = breaks_hist[1]
    xmax = breaks_hist[end]

    k = length(dens_hist)
    dens_hist1 = [dens_hist[1]]
    breaks_hist1 = [breaks_hist[1]]
    for j in 1:(k-1)
        if !(isapprox(dens_hist[j], dens_hist[j+1], atol=100*eps(), rtol=0.0))
            push!(dens_hist1, dens_hist[j+1])
            push!(breaks_hist1, breaks_hist[j+1])
            #println("$j, $(dens_hist[j])")
        end
    end
    push!(breaks_hist1, breaks_hist[end])

    k = length(dens_hist1)

    #println(breaks_hist1)
    #println(dens_hist1)
    mids = zeros(k)
    for j = 1:k
        mids[j] = 0.5*(breaks_hist1[j]+breaks_hist1[j+1])
    end

    # Now identify modes of the histogram density
    hist_modes = Float64[]
    mode_val = Float64[]
    if k == 1
        push!(hist_modes, 0.5*(breaks_hist1[1]+breaks_hist1[2]))
        push!(mode_val, dens_hist1[1])
    else
        if dens_hist1[1] > dens_hist1[2]
            push!(hist_modes, 0.5*(breaks_hist1[1]+breaks_hist1[2]))
            push!(mode_val, dens_hist1[1])
        end
        for j in 2:(k-1)
            if dens_hist1[j] > dens_hist1[j-1] && dens_hist1[j] > dens_hist1[j+1]
                push!(hist_modes, 0.5*(breaks_hist1[j]+breaks_hist1[j+1]))
                push!(mode_val, dens_hist1[j])
            end
        end
        if dens_hist1[end] > dens_hist1[end-1]
            push!(hist_modes, 0.5*(breaks_hist1[end-1]+breaks_hist1[end]))
            push!(mode_val, dens_hist1[end])
        end
    end

    # Finally, see if the modes of the histogram match any of those from the true density
    C = 0 # number of correctly identified modes
    for j in eachindex(true_modes)
        if minimum(abs.(true_modes[j] .- hist_modes)) < delta[j]
            C = C + 1
        end
    end
    if verbose
        println("Not identified: $(length(true_modes) - C), Spurious: $(length(hist_modes) - C)")
    end
    return (length(true_modes) - C) + (length(hist_modes) - C), hist_modes, mode_val, (length(hist_modes) - C)
end