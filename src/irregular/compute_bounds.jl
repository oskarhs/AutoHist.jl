# Compute optimal partition based on the output of the optimal partitioning algorithm
function compute_bounds_op(ancestor::Vector{Int}, grid::AbstractVector{<:Real}, k_max::Int)
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

# Compute optimal partition based on the output of the segment neighborhood algorithm
function compute_bounds_sn(ancestor::Matrix{Int}, grid::AbstractVector{<:Real}, k::Int)
    # Start recursion at k = k_opt, build the index vector "in reverse"
    L = Vector{Int64}(undef, k+1)
    L[k+1] = size(ancestor, 1)
    @inbounds for i in k:-1:1
        L[i] = ancestor[L[i+1], i]
    end
    bounds = grid[L .+ 1]
    return bounds
end