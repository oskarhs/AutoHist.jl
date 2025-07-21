# Φ corresponding to penB of Rozenholc et al. (2010)
@inline function phi_penB(i::Int, j::Int, N_cum::AbstractVector{<:Real}, grid::AbstractVector{<:Real})
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(N_bin > 0.0, N_bin * log(N_bin / len_bin), 0.0)
    return contrib
end

# Φ corresponding to Bayesian histogram with fixed concentration parameter a (not dep. on k)
@inline function phi_bayes(i::Int, j::Int, N_cum::AbstractVector{<:Real}, grid::AbstractVector{<:Real}, a::Real)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i] # Note: p0 = len_bin on the interval 0-1
    contrib = loggamma(a*len_bin + N_bin) - loggamma(a*len_bin) - N_bin * log(len_bin)
    return contrib
end

# Φ corresponding to penR of Rozenholc et al. (2010)
@inline function phi_penR(i::Int, j::Int, N_cum::AbstractVector{<:Real}, grid::AbstractVector{<:Real}, n::Real)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(N_bin > 0, N_bin * log(N_bin / len_bin) - 0.5 * N_bin / (n*len_bin), 0.0)
    return contrib
end

# Φ corresponding to Kullback-Leibler LOOCV
@inline function phi_KLCV(i::Int, j::Int, N_cum::AbstractVector{<:Real}, grid::AbstractVector{<:Real}, n::Real; minlength::Real=0.0)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    # NB! log throws error for negative args
    contrib = ifelse(N_bin ≥ 2 && len_bin ≥ minlength, N_bin * log(max((N_bin-1.0) / len_bin, 0.0)), -Inf64)
    return contrib
end

# Φ corresponding to L2 LOOCV
@inline function phi_L2CV(i::Int, j::Int, N_cum::AbstractVector{<:Real}, grid::AbstractVector{<:Real}, n::Real; minlength::Real=0.0)
    @inbounds N_bin = N_cum[j] - N_cum[i]
    @inbounds len_bin = grid[j] - grid[i]
    contrib = ifelse(len_bin ≥ minlength, ((n+1)/n * N_bin^2 - 2.0*N_bin) / len_bin, -Inf64)
    return contrib
end

function get_objective(rule::Symbol, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real, a::Real=5.0, use_min_length::Bool=false)
    if rule in [:pena, :penb, :nml]
        phi = let N_cum = N_cum, mesh = mesh
            f(i,j) = phi_penB(i, j, N_cum, mesh)
        end
    elseif rule == :bayes
        phi = let N_cum = N_cum, mesh = mesh, a = a
            f(i,j) = phi_bayes(i, j, N_cum, mesh, a)
        end
    elseif rule == :penr
        phi = let N_cum = N_cum, mesh = mesh, n = n
            f(i,j) = phi_penR(i, j, N_cum, mesh, n)
        end
    elseif rule == :klcv
        minlength = ifelse(use_min_length, log(n)^(1.5)/n, 0.0)
        phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
            f(i,j) = phi_KLCV(i, j, N_cum, mesh, n; minlength=minlength)
        end
    elseif rule == :l2cv
        minlength = ifelse(use_min_length, log(n)^(1.5)/n, 0.0)
        phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
            f(i,j) = phi_L2CV(i, j, N_cum, mesh, n; minlength=minlength)
        end
    end
    return phi
end

function get_phi(rule::RIH, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    phi = let N_cum = N_cum, mesh = mesh
        f(i,j) = phi_bayes(i, j, N_cum, mesh, rule.a)
    end
    return phi
end

function get_phi(rule::RMG_penB, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    phi = let N_cum = N_cum, mesh = mesh
        f(i,j) = phi_penB(i, j, N_cum, mesh)
    end
    return phi
end

function get_phi(rule::RMG_penA, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    phi = let N_cum = N_cum, mesh = mesh
        f(i,j) = phi_penB(i, j, N_cum, mesh)
    end
    return phi
end

function get_phi(rule::NML_I, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    phi = let N_cum = N_cum, mesh = mesh
        f(i,j) = phi_penB(i, j, N_cum, mesh)
    end
    return phi
end

function get_phi(rule::RMG_penR, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    phi = let N_cum = N_cum, mesh = mesh
        f(i,j) = phi_penR(i, j, N_cum, mesh, n)
    end
    return phi
end

function get_phi(rule::L2CV_I, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    minlength = ifelse(rule.use_min_length, log(n)^(1.5)/n, 0.0)
    phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
        f(i,j) = phi_L2CV(i, j, N_cum, mesh, n; minlength=minlength)
    end
    return phi
end

function get_phi(rule::KLCV_I, N_cum::AbstractVector{<:Real}, mesh::AbstractVector{<:Real}, n::Real)
    minlength = ifelse(rule.use_min_length, log(n)^(1.5)/n, 0.0)
    phi = let N_cum = N_cum, mesh = mesh, n = n, minlength=minlength
        f(i,j) = phi_KLCV(i, j, N_cum, mesh, n; minlength=minlength)
    end
    return phi
end

get_psi(rule::RIH, maxbins::Int, n::Int) = k-> rule.logprior(k) - logabsbinomial(maxbins-1, k-1)[1] + loggamma(rule.a) - loggamma(rule.a + n)
get_psi(rule::RMG_penB, maxbins::Int, n::Int) = k-> -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
get_psi(rule::RMG_penA, maxbins::Int, n::Int) = k-> -logabsbinomial(maxbins-1, k-1)[1] - k - 2.0*log(k) -
                        2.0 * sqrt(1.0*0.5*(k-1)*(logabsbinomial(maxbins-1, k-1)[1] + 2.0*log(k)))
get_psi(rule::RMG_penR, maxbins::Int, n::Int) = k-> -logabsbinomial(maxbins-1, k-1)[1] - k - log(k)^(2.5)
get_psi(rule::NML_I, maxbins::Int, n::Int) = k-> -( 0.5*k*log(0.5*n) - loggamma(0.5*k) +
                1.0/sqrt(n) * sqrt(2.0)*k/3.0 * exp(loggamma(0.5*k) - loggamma(0.5*k-0.5)) +
                1.0/n * ((3.0 + k*(k-2.0)*(2.0*k+1.0))/36.0 - k^2/9.0*exp(2.0*loggamma(0.5*k) - 2.0*loggamma(0.5*k-0.5)))
                ) - logabsbinomial(maxbins-1, k-1)[1]
get_psi(rule::L2CV_I, maxbins::Int, n::Int) = k-> 0.0
get_psi(rule::KLCV_I, maxbins::Int, n::Int) = k-> 0.0