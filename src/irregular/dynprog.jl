# The dynamic programming algorithm of Jackson et al. (2005).
function dynprog_quadratic(phi::Function, k_max::Int)
    cum_weight = Vector{Float64}(undef, k_max+1)
    cum_weight[1] = 0.0
    ancestor = Vector{Int64}(undef, k_max)
    ancestor[1] = 0
    weight = Matrix{Float64}(undef, k_max+1, k_max+1)

    function optimal_path!(ancestor, cum_weight, k)
        obj = Vector{Float64}(undef, k)
        @inbounds for i in 1:k
            obj[i] = cum_weight[i] + weight[i, k+1]
        end
        amax = argmax(obj)
        ancestor[k] = amax - 1
        cum_weight[k+1] = obj[amax]
    end

    # Compute weights for each possible interval
    for i in 1:k_max
        @simd for j in (i+1):(k_max+1)
            @inbounds weight[i, j] = phi(i, j)
        end
    end

    for k in 1:k_max
        optimal_path!(ancestor, cum_weight, k)
    end
    optimal = cum_weight[k_max+1] # Get optimum

    return optimal, ancestor
end

function compute_bounds_quadratic(ancestor::Vector{Int}, grid::AbstractVector{<:Real}, k_max::Int)
    # Start recursion at k_max (last cutpoint), then reverse the result to get the correct order
    L = Int64[k_max]
    j = k_max
    @inbounds while j > 0
        j = ancestor[j]
        L = push!(L, j)
    end
    bounds = grid[reverse(L) .+ 1]
    return bounds
end


# The dynamical programming algorithm of Kanazawa (1988).
function dynprog_cubic(phi::Function, k_max::Int)
    cum_weight = Matrix{Float64}(undef, k_max, k_max)
    ancestor = Matrix{Int64}(undef, k_max, k_max)
    ancestor[:, 1] .= 0
    ancestor[1,:] .= 0
    weight = Matrix{Float64}(undef, k_max+1, k_max+1)

    function optimal_path!(ancestor, cum_weight, k)
        ancestor0 = Vector{Int64}(undef, k_max-k+1) # these don't have to be reallocated
        cum_weight0 = Vector{Float64}(undef, k_max-k+1)

        @inbounds for i in k:k_max
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

# Compute optimal partition based on the output of the cubic DP algorithm
function compute_bounds_cubic(ancestor::Matrix{Int}, grid::AbstractVector{<:Real}, k::Int)
    L = Vector{Int64}(undef, k+1)
    L[k+1] = size(ancestor, 1)
    @inbounds for i in k:-1:1
        L[i] = ancestor[L[i+1], i]
    end
    bounds = grid[L .+ 1]
    return bounds
end