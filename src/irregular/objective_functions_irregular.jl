# Φ corresponding to penB of Rozenholc et al. (2010)
function phi_penB(i::Int, j::Int, N_cum::AbstractArray{<:Real}, grid::AbstractArray{<:Real})
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(N_bin > 0.0, N_bin * log(N_bin / len_bin), 0.0)
    return contrib
end

# Φ corresponding to Bayesian histogram with fixed concentration parameter a (not dep. on k)
function phi_bayes(i::Int, j::Int, N_cum::AbstractArray{<:Real}, grid::AbstractArray{<:Real}, a::Real)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i] # Note: p0 = len_bin on the interval 0-1
    contrib = loggamma(a*len_bin + N_bin) - loggamma(a*len_bin) - N_bin * log(len_bin)
    return contrib
end

# Φ corresponding to penR of Rozenholc et al. (2010)
function phi_penR(i::Int, j::Int, N_cum::AbstractArray{<:Real}, grid::AbstractArray{<:Real}, n::Real)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(N_bin > 0, N_bin * log(N_bin / len_bin) - 0.5 * N_bin / (n*len_bin), 0.0)
    return contrib
end

# Φ corresponding to Kullback-Leibler LOOCV
function phi_KLCV(i::Int, j::Int, N_cum::AbstractArray{<:Real}, grid::AbstractArray{<:Real}, n::Real; minlength::Real=0.0)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    # NB! log throws error for negative args
    contrib = ifelse(N_bin ≥ 2 && len_bin ≥ minlength, N_bin * log(max((N_bin-1.0) / len_bin, 0.0)), -Inf64)
    return contrib
end

# Φ corresponding to L2 LOOCV
function phi_L2CV(i::Int, j::Int, N_cum::AbstractArray{<:Real}, grid::AbstractArray{<:Real}, n::Real; minlength::Real=0.0)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(len_bin ≥ minlength, ((n+1)/n * N_bin^2 - 2.0*N_bin) / len_bin, -Inf64)
    return contrib
end