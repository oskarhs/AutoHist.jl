# Uses the greedy algorithm of Rozenholc et al. to construct a coarser grid to make the algorithm better suited for large datasets.
# To be used prior to running the dynamic programming algorithm.
function greedy_grid(N_cum::AbstractVector{<:Real}, finestgrid::AbstractVector{<:Real}, maxbins::Int, gr_maxbins::Int)
    # Update increments between the values i and j
    n = N_cum[end]
    compute_loglik_increments! = let N_cum = N_cum, finestgrid=finestgrid, n = n
        function (incr, i, j)
            @inbounds if finestgrid[i] < finestgrid[j]
                # Log-likelihood contribution 
                @inbounds loglik_old = (N_cum[j] - N_cum[i]) * log((N_cum[j]-N_cum[i])/(n*(finestgrid[j]-finestgrid[i])))
                @turbo for l = (i+1):(j-1)
                    loglik_new = (N_cum[l] - N_cum[i]) * log((N_cum[l]-N_cum[i])/(n*(finestgrid[l]-finestgrid[i]))) +
                            (N_cum[j] - N_cum[l]) * log((N_cum[j]-N_cum[l])/(n*(finestgrid[j]-finestgrid[l])))
                    incr[l] = loglik_new - loglik_old
                end
            end
        end
    end
    grid_ind = fill(false, maxbins+1) # Array of booleans storing which indices to use
    grid_ind[1] = true
    grid_ind[end] = true
    incr = Array{Float64}(undef, maxbins+1) # Array of increments from splitting at index d at each step
    incr[1] = -Inf # Would create a bin of lebesgue measure 0
    incr[end] = -Inf

    # First iteration
    i = 1
    j = maxbins+1
    compute_loglik_increments!(incr, i, j)
    num_bins = 1

    # Terminate when num_bins has reached the limit or 
    # when increases to the log-likelihood are no longer possible
    while num_bins < gr_maxbins && maximum(incr) > 0
        # Update increments for indices in (i, d) and (d, j)
        d = argmax(incr)
        @inbounds grid_ind[d] = true # Include finestgrid[d] in the grid
        @inbounds incr[d] = -Inf # Included in grid
        num_bins = num_bins + 1

        # Set i to maximal index < than d s.t. grid_ind[i] == true
        i = findlast(@views grid_ind[1:d-1])
        # Set j to minimal index > than d s.t. grid_ind[j] == true
        j = findfirst(@views grid_ind[d+1:end]) + d

        compute_loglik_increments!(incr, i, d)
        compute_loglik_increments!(incr, d, j)
    end
    return grid_ind
end