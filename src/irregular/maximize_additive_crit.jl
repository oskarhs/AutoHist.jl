# Maximize an additive criterion for irregular histograms
function maximize_additive_crit(x::AbstractVector{T}, rule::A, xmin::T, xmax::T, closed::Symbol) where {T <: Real, A<:AbstractIrregularRule}
    y = @. (x - xmin) / (xmax - xmin)
    n = length(x)
    grid = rule.grid

    if grid == :data
        sort!(y)
        z = unique(y)
        maxbins = length(z)-1
    elseif rule.use_default_maxbins
        maxbins = min(n, ceil(Int, 4.0*n/log(n)^2))
    else
        maxbins = rule.maxbins
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
        finestgrid[1] = 0.0
        finestgrid[2:end-1] = z[2:end-1]
        finestgrid[end] = 1.0
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    elseif grid == :regular
        finestgrid[1:end] = LinRange(0.0, 1.0, maxbins+1)
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    elseif grid == :quantile
        finestgrid[1] = 0.0
        finestgrid[end] = 1.0
        finestgrid[2:end-1] = quantile(y, LinRange(1.0/maxbins, 1.0-1.0/maxbins, maxbins-1); sorted=false)
        N_cum[2:end] = cumsum(bin_irregular(y, finestgrid, closed == :right))
    end

    if rule.alg.greedy
        if rule.alg.use_default_gr_maxbins
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