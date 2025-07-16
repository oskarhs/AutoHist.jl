# Greedy, linear-time pruned dynamic programming algorithm
#= function optimal_partitioning_greedy(phi::Function, k_max::Int, max_cand::Int=15)
    cum_weight = Vector{Float64}(undef, k_max+1)
    cum_weight[1] = 0.0
    ancestor = Vector{Int64}(undef, k_max)
    ancestor[1] = 0

    function optimal_pruned_path!(ancestor, cum_weight, ctps_cand, k)
        # optimize via the DP recursion for the non-pruned candidate changepoints
        obj = Vector{Float64}(undef, length(ctps_cand))
        @inbounds for (j,i) in enumerate(ctps_cand)
            obj[j] = cum_weight[i+1] + phi(i+1, k+1)
        end
        amax = argmax(obj)
        @inbounds ancestor[k] = ctps_cand[amax]
        @inbounds cum_weight[k+1] = obj[amax]

        # prune the remaining candidate changepoints
        noprune = partialsortperm(obj, 1:min(max_cand-1, length(obj)), rev=true)
        R = ctps_cand[noprune]
        push!(R, k)
        return R
    end

    # vector of candidate cutpoints at each time stp
    ctps_cand = Int64[0]

    for k in 1:k_max
        ctps_cand = optimal_pruned_path!(ancestor, cum_weight, ctps_cand, k)
    end
    optimal = cum_weight[k_max+1] # Get optimum

    return optimal, ancestor
end =#

# Greedy, quadratic-time pruned dynamic programming algorithm
function segment_neighborhood_greedy(phi::Function, k_max::Int, max_cand::Int=15)
    cum_weight = Matrix{Float64}(undef, k_max, k_max)
    ancestor = Matrix{Int64}(undef, k_max, k_max)
    ancestor[:, 1] .= 0
    ancestor[1,:] .= 0
    weight = Matrix{Float64}(undef, k_max+1, k_max+1)

    function optimal_pruned_path!(ancestor, cum_weight, k)
        ancestor0 = Vector{Int64}(undef, k_max-k+1) # these don't have to be reallocated
        cum_weight0 = Vector{Float64}(undef, k_max-k+1)

        cpts_cand = Int64[k]
        sizehint!(cpts_cand, max_cand)
        @inbounds for i in k:k_max
            # Compute optimum among the non-pruned candidates
            obj = Vector{Float64}(undef, length(cpts_cand))
            for (l,j) in enumerate(cpts_cand)
                #obj[l] = cum_weight[j, k-1] + weight[j+1, i+1]
                obj[l] = cum_weight[j-1, k-1] + weight[j, i+1]
            end
            amax = argmax(obj)
            ancestor0[i-k+1] = cpts_cand[amax]
            cum_weight0[i-k+1] = obj[amax]

            # Greedy pruning step
            num_cand = min(length(obj), max_cand-1)
            if num_cand < max_cand - 1
                push!(cpts_cand, i)        # add candidates until max_cand is reached
            else
                cpts_cand[argmin(obj)] = i # replace worst candidate
            end
        end
        @inbounds ancestor[k:k_max, k] = ancestor0
        @inbounds cum_weight[k:k_max, k] = cum_weight0
    end

    # Compute weights
    for i in 1:k_max
        @simd for j in (i+1):(k_max+1)
            @inbounds weight[i, j] = phi(i, j)
        end
    end
    @inbounds cum_weight[:,1] = @views weight[1,2:k_max+1]
    for k in 2:k_max
        optimal_pruned_path!(ancestor, cum_weight, k)
    end
    optimal = @views cum_weight[k_max,:] # Get weight function for each partition

    return optimal, ancestor
end