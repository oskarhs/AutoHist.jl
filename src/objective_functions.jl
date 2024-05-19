using Distributions, SpecialFunctions


# Unnormalized log-posterior with Geometric(τ) prior on k, Dirichlet(ap_0) prior on p
function logposterior_k(N, k, a, p0, n, logprior)
    logpost = loggamma(a) - loggamma(a + n) + n*log(k) + logprior(k)
    #logpost = loggamma(a) - loggamma(a + n) + n*log(k) - log(1+k)
    @inbounds for j = 1:k
        logpost = logpost + loggamma(a*p0[j] + N[j]) - loggamma(a*p0[j])
    end
    return logpost
end

function compute_BIC(N, k, n)
    BIC = -2.0*n*log(k) + k * log(n)
    @inbounds for i in eachindex(N)
        if N[i] > 0
            BIC = BIC - 2.0*N[i] * log(N[i]/n)
        end
    end
    return BIC
end

function compute_AIC(N, k, n)
    AIC = -2.0*n*log(k) + 2.0*k
    @inbounds for i in eachindex(N)
        if N[i] > 0
            AIC = AIC - 2.0*N[i] * log(N[i]/n)
        end
    end
    return AIC
end

# Optimization criterion due to Birge and Rozenholc
# Note that the sign is flipped relative to AIC to stay consistent with their original paper
function compute_BR(N, k, n)
    BR = n*log(k) - k - log(k)^(2.5)
    @inbounds for i in eachindex(N)
        if N[i] > 0
            BR = BR + N[i] * log(N[i]/n)
        end
    end
    return BR
end

# Objective (maximization) for regular histograms based on minimum description length (Hall and Hannan)
function compute_MDL(N, k, n)
    MDL = -(n-0.5*k)*log(n-0.5*k) + n*log(k) - 0.5*k*log(n)
    @inbounds for i in eachindex(N)
        if N[i] > 0
            MDL = MDL + (N[i]-0.5) * log(N[i]-0.5)
        end
    end
    return MDL
end

# Objective (maximization) for regular histograms based on the stochastic complexity criterion ()
function compute_SC(N, k, n)
    SC = n*log(k) + loggamma(k) - loggamma(k+n)
    @inbounds for i in eachindex(N)
        if N[i] > 0
            SC = SC + loggamma(N[i]+1)
        end
    end
    return SC
end

# Objective (maximization) for regular histograms based on Kullback-Leibler Cross Validation
function compute_KLCV(N, k, n)
    KLCV = n*log(k)
    @inbounds for i in eachindex(N)
        if N[i] > 0
            KLCV = KLCV + N[i] * log(N[i]-1.0)
        end
    end
    return KLCV
end

# Objective (maximization) for regular histograms based on Normalized Maximum Likelihood
function compute_NML(N, K, n)
    NML = n*log(k) + 
          0.5*(k-1.0)*log(0.5*n) + log(sqrt(π)/gamma(0.5*k)) + 
          sqrt(2.0)*k*gamma(0.5*k)/(3.0*sqrt(n)*gamma(0.5*(k-1))) + 
          1.0/n*( ( 3.0 + k*(k-2.0)*(2.0*k+1.0) )/ 36.0 - k^2 * gamma(0.5*k)^2/(9.0*gamma(0.5*(k-1.0))^2) )
    for i in eachindex(N)
        NML = NML + N[i]*log(N[i]/n)
    end
    return NML
end