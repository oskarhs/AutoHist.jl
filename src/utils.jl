# Specialized function to perform regular binning quickly for known min, max.
# Supports both left- and right-inclusive binning
function bin_regular(x::AbstractVector{<:Real}, xmin::Real, xmax::Real, k::Int, right::Bool)
    R = xmax - xmin
    bincounts = zeros(Float64, k)
    edges_inc = k/R
    if right
        for val in x
            idval = min(max(1, ceil(Int, (val-xmin)*edges_inc)), k)
            bincounts[idval] += 1.0
        end
    else
        for val in x
            idval = max(min(k-1, floor(Int, (val-xmin)*edges_inc)) + 1, 1)
            bincounts[idval] += 1.0
        end
    end
    return bincounts
end

function bin_irregular(x::AbstractVector{<:Real}, edges::AbstractVector{<:Real}, right::Bool)
    k = length(edges)-1
    bincounts = zeros(Float64, length(edges)-1)
    if right
        for val in x
            idval = min(max(1, searchsortedfirst(edges, val) - 1), k)
            bincounts[idval] += 1.0
        end
    else
        for val in x
            idval = max(min(k, searchsortedlast(edges, val)), 1)
            bincounts[idval] += 1.0
        end
    end
    return bincounts
end

function bin_irregular_int(x::AbstractVector{<:Real}, edges::AbstractVector{<:Real}, right::Bool)
    k = length(edges)-1
    bincounts = zeros(Int64, length(edges)-1)
    if right
        for val in x
            idval = min(max(1, searchsortedfirst(edges, val) - 1), k)
            bincounts[idval] += 1
        end
    else
        for val in x
            idval = max(min(k, searchsortedlast(edges, val)), 1)
            bincounts[idval] += 1
        end
    end
    return bincounts
end

# Merge two sorted vectors in a way such that the resulting vector is also sorted.
function merge_sorted_vec(x1::AbstractVector{T}, x2::AbstractVector{T}) where T
    k1 = length(x1)
    k2 = length(x2)
    new = Vector{T}(undef, k1+k2)
    j1 = 1
    j2 = 1
    i = 1
    @inbounds while j1 ≤ k1 && j2 ≤ k2
        new[i], j1, j2 = ifelse(
            x1[j1] < x2[j2],
            (x1[j1], j1+1, j2),
            (x2[j2], j1, j2+1)
        )
        i = i+1
    end
    @inbounds while j1 ≤ k1
        new[i] = x1[j1]
        j1 = j1+1
        i = i+1
    end
    @inbounds while j2 ≤ k2
        new[i] = x2[j2]
        j2 = j2+1
        i = i+1
    end
    return new
end

# Compute the Hellinger distance between two histogram probability densities
function hellinger_distance(h1::AutomaticHistogram, h2::AutomaticHistogram)
    breaks1 = h1.breaks
    dens1 = h1.density
    breaks2 = h2.breaks
    dens2 = h2.density
    disc = merge_sorted_vec(breaks1, breaks2) # union of points of discontinuity for the two histograms
    m = length(disc)-1
    hell = 0.0
    for j = 1:m
        bin_ind1 = searchsortedfirst(breaks1, disc[j]+10*eps()) - 1
        bin_ind2 = searchsortedfirst(breaks2, disc[j]+10*eps()) - 1
        if bin_ind1 == 0 || bin_ind1 == (length(dens1)+1)
            d1 = 0.0
        else
            d1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            d2 = 0.0
        else
            d2 = dens2[bin_ind2]
        end
        hell += (disc[j+1]-disc[j])*(sqrt(d1)-sqrt(d2))^2
    end
    return sqrt(hell)
end

# Compute the lp distance between two histogram probability densities
function lp_distance(h1::AutomaticHistogram, h2::AutomaticHistogram, p::Real)
    breaks1 = h1.breaks
    dens1 = h1.density
    breaks2 = h2.breaks
    dens2 = h2.density
    disc = merge_sorted_vec(breaks1, breaks2) # union of points of discontinuity for the two histograms
    m = length(disc)-1
    lp = 0.0
    for j = 1:m
        bin_ind1 = searchsortedfirst(breaks1, disc[j]+10*eps()) - 1
        bin_ind2 = searchsortedfirst(breaks2, disc[j]+10*eps()) - 1
        if bin_ind1 == 0 || bin_ind1 == (length(dens1)+1)
            d1 = 0.0
        else
            d1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            d2 = 0.0
        else
            d2 = dens2[bin_ind2]
        end
        lp += (disc[j+1]-disc[j])*abs(d1-d2)^p
    end
    return lp^(1.0/p)
end

function supremum_distance(h1::AutomaticHistogram, h2::AutomaticHistogram)
    breaks1 = h1.breaks
    dens1 = h1.density
    breaks2 = h2.breaks
    dens2 = h2.density
    disc = merge_sorted_vec(breaks1, breaks2) # union of points of discontinuity for the two histograms
    m = length(disc)-1
    sup = 0.0
    for j = 1:m
        bin_ind1 = searchsortedfirst(breaks1, disc[j]+10*eps()) - 1
        bin_ind2 = searchsortedfirst(breaks2, disc[j]+10*eps()) - 1
        if bin_ind1 == 0 || bin_ind1 == (length(dens1)+1)
            d1 = 0.0
        else
            d1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            d2 = 0.0
        else
            d2 = dens2[bin_ind2]
        end
        sup = max(abs(d1-d2), sup)
    end
    return sup
end

# Compute K(f₁, f₂) (asymmetric)
function kl_divergence(h1::AutomaticHistogram, h2::AutomaticHistogram)
    breaks1 = h1.breaks
    dens1 = h1.density
    breaks2 = h2.breaks
    dens2 = h2.density
    disc = merge_sorted_vec(breaks1, breaks2) # union of points of discontinuity for the two histograms
    m = length(disc)-1
    kl = 0.0
    for j = 1:m
        bin_ind1 = searchsortedfirst(breaks1, disc[j]+10*eps()) - 1
        bin_ind2 = searchsortedfirst(breaks2, disc[j]+10*eps()) - 1
        if bin_ind1 == 0 || bin_ind1 == (length(dens1)+1)
            d1 = 0.0
        else
            d1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            d2 = 0.0
        else
            d2 = dens2[bin_ind2]
        end
        kl += (disc[j+1]-disc[j])*ifelse(d1 == 0.0, 0.0, d1*log(d1) - d1*log(d2))
    end
    return kl
end

# Check that maxbins is positive and compute default if applicable
#= function get_maxbins_regular(maxbins::Union{Symbol, Int}, n::Int)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    k_max = ifelse(maxbins == :default, min(ceil(Int, 4.0*n / log(n)^2), 1000), maxbins)
    return k_max
end =#