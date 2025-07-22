# Uses the greedy algorithm of Rozenholc et al. (2010) to construct a coarser grid to make the algorithm better suited for large datasets.
# To be used prior to running the dynamic programming algorithm.
function greedy_grid(phi::Function, maxbins::Int, gr_maxbins::Int)
    # Update increments between the values i and j
    function compute_objective_increments!(incr, i, j)
        # objective contribution
        @inbounds obj_old = phi(i, j)
        @inbounds for l = (i+1):(j-1)
            obj_new = phi(i, l)
            obj_new += phi(l, j)        
            incr[l] = obj_new - obj_old
        end
    end
    grid_ind = Int64[1, maxbins+1] # Indices included in partiiton
    sizehint!(grid_ind, gr_maxbins)
    incr = Vector{Float64}(undef, maxbins+1) # Array of increments from splitting at index d at each step
    incr[1] = -Inf # Would create a bin of lebesgue measure 0
    incr[end] = -Inf

    # First iteration
    i = 1
    j = maxbins+1
    compute_objective_increments!(incr, i, j)
    num_bins = 1
    dctn = Dict{Int64, Tuple{Int64, Int64, Float64}}() # store next, argmax, max over interval to next
    max_new, amax_new = findmax(incr)
    dctn[1] = (maxbins+1, amax_new, max_new)

    # Terminate when num_bins has reached the limit
    while num_bins < gr_maxbins
        # Update increments for indices in (i, d) and (d, j)
        curr_val, next_val, argmax_val = findmax_dctn(dctn)
        push!(grid_ind, argmax_val)
        @inbounds incr[argmax_val] = -Inf # Included in grid
        num_bins = num_bins + 1

        # Update increments for  curr_val < i < argmax val
        compute_objective_increments!(incr, curr_val, argmax_val)
        # Update increments for  argmax_val < i < next_val
        compute_objective_increments!(incr, argmax_val, next_val)

        # Update dictionary entries so that curr_val -> argmax_val
        max_new, amax_new = if curr_val+1 < argmax_val findmax(@views incr[curr_val+1:argmax_val-1]) else (-Inf, 0) end
        dctn[curr_val] = (argmax_val, amax_new+curr_val, max_new)
        max_new, amax_new = if argmax_val+1 < next_val findmax(@views incr[argmax_val+1:next_val-1]) else (-Inf, 0) end
        dctn[argmax_val] = (next_val, amax_new+argmax_val, max_new)
    end
    println(grid_ind)
    return sort!(grid_ind)
end

# Find dictionary entry yielding the maximal increment
function findmax_dctn(dctn::Dict{Int64, Tuple{Int64, Int64, Float64}})
    max_val = -Inf
    argmax_val = 0
    curr_val = 0
    next_val = 0
    @inbounds for (key, value) in dctn
        if value[3] > max_val
            max_val = value[3]
            argmax_val = value[2]
            next_val = value[1]
            curr_val = key
        end
    end
    return curr_val, next_val, argmax_val
end