# The dynamic programming algorithm of Jackson et al. (2005).
function dynprog_quadratic(phi::Function, k_max::Int)
    cum_weight = Vector{Float64}(undef, k_max+1)
    cum_weight[1] = 0.0
    ancestor = Vector{Int64}(undef, k_max)
    ancestor[1] = 0
    weight = Matrix{Float64}(undef, k_max+1, k_max+1)

    function optimal_path!(ancestor, cum_weight, k)
        #ancestor0 = Array{Int64}(undef, k_max-k+1) # these don't have to be reallocated
        obj = Vector{Float64}(undef, k)
        for i = 1:k
            obj[i] = cum_weight[i] + weight[i, k+1]
        end
        temp = argmax(obj)
        ancestor[k] = temp - 1
        cum_weight[k+1] = obj[temp]
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

# Compute optimal partition based on the output of the cubic DP algorithm
function compute_bounds_cubic(ancestor::Matrix{Int}, grid::AbstractVector{<:Real}, k::Int)
    L = Array{Int64}(undef, k+1)
    L[k+1] = size(ancestor, 1)
    @inbounds for i = k:-1:1
        L[i] = ancestor[L[i+1], i]
    end
    bounds = grid[L .+ 1]
    return bounds
end

#= function main()
    n = 10^5
    x = rand(Xoshiro(1), Beta(2.0, 2.0), n)
    k_max = 500
    h = fit(Histogram, x, LinRange(0, 1, k_max+1))
    println(length(h.weights))
    N = h.weights
    grid = h.edges[1]
    N_cum = Vector{Int64}(undef, k_max+1)
    N_cum[1] = 0; N_cum[2:end] = cumsum(N)

    #= function phi(i::Int, j::Int)
        @inbounds N_bin = N_cum[j] - N_cum[i]
        @inbounds len_bin = grid[j] - grid[i]
        contrib = ((n+1)/n * N_bin^2 - 2.0*N_bin) / len_bin
        return contrib
    end =#
    minlength = 0.02
    function phi1(i::Int, j::Int)
        @inbounds N_bin = N_cum[j] - N_cum[i]
        @inbounds len_bin = grid[j] - grid[i]
        if len_bin > minlength
            contrib = ((n+1)/n * N_bin^2 - 2.0*N_bin) / len_bin
        else
            contrib = -Inf64
        end
        return contrib
    end 
    
    opt1, anc1 = dynprog_quadratic(phi1, k_max)
    opt2, anc2 = dynprog_cubic(phi1, k_max)

    k_opt = argmax(opt2)

    #println(opt1, " ", maximum(opt2))
    @assert opt1 ≈ maximum(opt2)
    println(compute_bounds_quadratic(anc1, grid, k_max))
    #println(anc1, " ", anc2)

    # Sanity check 
#=     part = compute_bounds_quadratic(anc1, grid, k_max)
    tot = 0.0
    for l = 1:length(part)-1
        tot += phi1(part[l], part[l+1])
    end =#
    #println(tot)

    @assert compute_bounds_quadratic(anc1, grid, k_max) == compute_bounds_cubic(anc2, grid, k_opt)

    @btime dynprog_quadratic($phi1, $k_max)
    @btime dynprog_cubic($phi1, $k_max)
    return nothing
end

main() =#