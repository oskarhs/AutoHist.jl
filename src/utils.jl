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

#= function bin_regular_int(x::AbstractVector{<:Real}, xmin::Real, xmax::Real, k::Int, right::Bool)
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
end =#

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
function get_maxbins_regular(maxbins::Union{Symbol, Int}, n::Int)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins bins must be positive."))
    end
    k_max = min(ifelse(maxbins == :default, ceil(Int, 4.0*n / log(n)^2), maxbins), 1000)
    return k_max
end

# Check that maxbins is positive and compute default if applicable
function get_maxbins_irregular!(y::AbstractVector{<:Real}, maxbins::Union{Symbol, Int}, n::Int, grid::Symbol)
    if typeof(maxbins) <: Symbol && maxbins != :default
        throw(ArgumentError("maxbins must either be a positive integer or :default."))
    elseif typeof(maxbins) <: Int && maxbins < 1             # maximal number of bins must be positive
        throw(DomainError("maxbins must be positive."))
    end
    if !(grid in [:data, :regular, :quantile])
        throw(ArgumentError("The supplied grid, :$grid, is not supported. The grid kwarg must be one of :data, :regular or :quantile"))
    end

    k_max = if grid == :data
        sort!(y)
        z = unique(y)
        length(z)-1
    else
        ifelse(maxbins == :default, min(n, ceil(Int, 4.0*n/log(n)^2)), maxbins)     
    end
    return k_max
end


function maximize_additive_crit(x::AbstractVector{T}, rule::A, xmin::T, xmax::T, closed::Symbol) where {T <: Real, A<:AbstractIrregularRule}
    y = @. (x - xmin) / (xmax - xmin)
    n = length(x)

    if typeof(rule.grid) <: Symbol && rule.grid in (:data, :regular, :quantile)
        grid = rule.grid
    else
        throw(ArgumentError("Invalid grid option: $(rule.grid)"))
    end

    if grid == :data
        sort!(y)
        z = unique(y)
        maxbins = length(z)-1
    elseif typeof(rule.maxbins) <: Symbol && rule.maxbins == :default
        maxbins = min(n, ceil(Int, 4.0*n/log(n)^2))
    elseif isa(rule.maxbins, Int) && rule.maxbins > 0
        maxbins = rule.maxbins
    else
        if typeof(rule.maxbins) <: Symbol && rule.maxbins != :default
            throw(ArgumentError("maxbins must either be a positive integer or :default."))
        else
            throw(DomainError("maxbins must be positive."))
        end
    end

    if hasfield(typeof(rule), :use_min_length)
        use_min_length = rule.use_min_length
    else
        use_min_length = false
    end

    finestgrid = Vector{Float64}(undef, maxbins+1)
    N_cum = Vector{Float64}(undef, length(finestgrid)) # cumulative cell counts
    N_cum[1] = 0.0
    if grid == :data
        finestgrid[1] = -eps()
        finestgrid[2:end-1] = z[2:end-1]
        finestgrid[end] = 1.0 + eps()
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    elseif grid == :regular
        N_cum[2:end] = cumsum(bin_regular(y, 0.0, 1.0, maxbins, closed == :right))
        finestgrid[1:end] = LinRange(0.0, 1.0, maxbins+1)
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0+eps()
    elseif grid == :quantile
        finestgrid[1] = -eps()
        finestgrid[end] = 1.0 + eps()
        finestgrid[2:end-1] = quantile(y, LinRange(1.0/maxbins, 1.0-1.0/maxbins, maxbins-1); sorted=false)
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    end

    if rule.alg.greedy
        if rule.alg.gr_maxbins == :default
            gr_maxbins = ifelse(
                typeof(rule.alg) == OptPart,
                min(maxbins, max(ceil(Int, n^(1.0/2.0)), 3000)),
                min(maxbins, max(ceil(Int, n^(1.0/3.0)), 500))
            )
        else
            gr_maxbins = min(maxbins, rule.alg.gr_maxbins)
        end
        if gr_maxbins == maxbins # no reason to use the greedy algorithm in this case
            k_max = maxbins
            mesh = finestgrid
        else
            grid_ind = greedy_grid(
                get_phi(ifelse(typeof(rule) != KLCV_I, rule, RMG_penB()), N_cum, finestgrid, n), 
                maxbins,
                gr_maxbins
            )
            mesh = finestgrid[grid_ind]
            k_max = length(mesh) - 1

            # Update bin counts to the newly constructed grid
            N_cum = N_cum[grid_ind]
        end
    else
        k_max = maxbins
        mesh = finestgrid
    end
    phi = get_phi(rule, N_cum, mesh, n)

    return dynprog(rule.alg, rule, phi, mesh, k_max, maxbins, n)
end