# The dynamical programming algorithm of Kanazawa (1988)
function dynamic_algorithm(phi::Function, k_max::Int)
    cum_weight = Matrix{Float64}(undef, k_max, k_max)
    ancestor = Array{Int64}(undef, k_max, k_max)
    ancestor[:, 1] .= 0
    ancestor[1,:] .= 0
    weight = Matrix{Float64}(undef, k_max+1, k_max+1)

    function optimal_path!(ancestor, cum_weight, k)
        ancestor0 = Array{Int64}(undef, k_max-k+1) # these don't have to be reallocated
        cum_weight0 = Array{Float64}(undef, k_max-k+1)

        @inbounds for i = k:k_max
            obj = @views cum_weight[(k-1):(i-1), k-1] .+ weight[k:i, i+1]
            ancestor0[i-k+1] = argmax(obj)
            cum_weight0[i-k+1] = obj[ancestor0[i-k+1]]
        end
        @inbounds ancestor[k:k_max, k] = ancestor0 .+ (k-2)
        @inbounds cum_weight[k:k_max, k] = cum_weight0
    end

    # Compute weights for each possible interval
    for i in 1:k_max
        @simd for j in (i+1):(k_max+1)
            @inbounds weight[i, j] = phi(i, j)
        end
    end

    # Compute cumulative weights
    @inbounds cum_weight[:,1] = @views weight[1,2:k_max+1]
    for k in 2:k_max
        optimal_path!(ancestor, cum_weight, k)
    end
    optimal = @views cum_weight[k_max,:] # Get weight function for each partition

    return optimal, ancestor
end

# Compute optimal partition based on the output of the DP algorithm
function compute_bounds(ancestor::Matrix{Int}, grid::AbstractVector{<:Real}, k::Int)
    L = Array{Int64}(undef, k+1)
    L[k+1] = size(ancestor, 1)
    @inbounds for i = k:-1:1
        L[i] = ancestor[L[i+1], i]
    end
    bounds = grid[L .+ 1]
    return bounds
end