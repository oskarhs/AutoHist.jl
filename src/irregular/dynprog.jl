# The optimal partitioning (dynamic programming) algorithm of Jackson et al. (2005).
function optimal_partitioning(phi::Function, k_max::Int)
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

# The segment neighborhood (dynamical programming) algorithm of Kanazawa (1988).
function segment_neighborhood(phi::Function, k_max::Int)
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


function dynprog(alg::OptPart, rule::AbstractIrregularRule, phi::Function, mesh::AbstractVector{<:Real}, k_max::Int, maxbins::Int, n::Int)
    optimal, ancestor = optimal_partitioning(phi, k_max)
    bin_edges_norm = compute_bounds_op(ancestor, mesh, k_max)
    return bin_edges_norm
end


function dynprog(alg::SegNeig, rule::AbstractIrregularRule, phi::Function, mesh::AbstractVector{<:Real}, k_max::Int, maxbins::Int, n::Int)
    optimal, ancestor = segment_neighborhood(phi, k_max)
    psi = get_psi(rule, maxbins, n)
    dep_k = Vector{Float64}(undef, k_max)
    @simd for k in 1:k_max
        @inbounds dep_k[k] = psi(k)
    end
    k_opt = argmax(optimal + dep_k)
    bin_edges_norm = compute_bounds_sn(ancestor, mesh, k_opt)
    return bin_edges_norm
end