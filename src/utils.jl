# Specialized function to perform regular binning quickly for known min, max.
# Supports both left- and right-inclusive binning
function bin_regular(x::AbstractVector{<:Real}, xmin::Real, xmax::Real, k::Int, right::Bool)
    R = xmax - xmin
    bincounts = zeros(Float64, k)
    edges_inc = k/R
    if right
        for val in x
            idval = max(0, floor(Int, (val-xmin)*edges_inc-10.0*edges_inc*eps())) + 1
            @inbounds bincounts[idval] += 1.0
        end
    else
        for val in x
            idval = min(k-1, floor(Int, (val-xmin)*edges_inc+10.0*edges_inc*eps())) + 1
            @inbounds bincounts[idval] += 1.0
        end
    end
    return bincounts
end

function bin_regular_int(x::AbstractVector{<:Real}, xmin::Real, xmax::Real, k::Int, right::Bool)
    R = xmax - xmin
    bincounts = zeros(Int64, k)
    edges_inc = k/R
    if right
        for val in x
            idval = max(0, floor(Int, (val-xmin)*edges_inc-10.0*edges_inc*eps())) + 1
            @inbounds bincounts[idval] += 1
        end
    else
        for val in x
            idval = min(k-1, floor(Int, (val-xmin)*edges_inc+10.0*edges_inc*eps())) + 1
            @inbounds bincounts[idval] += 1
        end
    end
    return bincounts
end

function bin_irregular(x::AbstractVector{<:Real}, edges::AbstractVector{<:Real}, right::Bool)
    bincounts = zeros(Float64, length(edges)-1)
    if right
        for val in x
            idval = max(1, searchsortedfirst(edges, val) - 1)
            @inbounds bincounts[idval] += 1.0
        end
    else
        k = length(edges)-1
        for val in x
            idval = min(k, searchsortedlast(edges, val))
            @inbounds bincounts[idval] += 1.0
        end
    end
    return bincounts
end

function bin_irregular_int(x::AbstractVector{<:Real}, edges::AbstractVector{<:Real}, right::Bool)
    bincounts = zeros(Int64, length(edges)-1)
    if right
        for val in x
            idval = max(1, searchsortedfirst(edges, val) - 1)
            @inbounds bincounts[idval] += 1
        end
    else
        k = length(edges)-1
        for val in x
            idval = min(k, searchsortedlast(edges, val))
            @inbounds bincounts[idval] += 1
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
    while j1 ≤ k1 && j2 ≤ k2
        new[i], j1, j2 = ifelse(
            x1[j1] < x2[j2],
            (x1[j1], j1+1, j2),
            (x2[j2], j1, j2+1)
        )
        i = i+1
    end
    while j1 ≤ k1
        new[i] = x1[j1]
        j1 = j1+1
        i = i+1
    end
    while j2 ≤ k2
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
            h1 = 0.0
        else
            h1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            h2 = 0.0
        else
            h2 = dens2[bin_ind2]
        end
        hell += (disc[j+1]-disc[j])*(sqrt(h1)-sqrt(h2))^2
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
            h1 = 0.0
        else
            h1 = dens1[bin_ind1]
        end
        if bin_ind2 == 0 || bin_ind2 == (length(dens2)+1)
            h2 = 0.0
        else
            h2 = dens2[bin_ind2]
        end
        lp += (disc[j+1]-disc[j])*abs(h1-h2)^p
    end
    return lp^(1.0/p)
end