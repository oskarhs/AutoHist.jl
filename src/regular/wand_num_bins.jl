# Compute the number of bins using Wand's (1997) criterion
# Based on code from Julia's StatsPlots library: https://github.com/JuliaPlots/StatsPlots.jl/blob/master/src/hist.jl 
# Originally Ported from R code located here https://github.com/cran/KernSmooth/tree/master/R
function wand_num_bins(x, level=2, scalest = :minim, gridsize = 401, range_x = extrema(x))
    n = length(x)
    minx, maxx = range_x
    gpoints = range(minx, stop = maxx, length = gridsize)
    gcounts = linbin(x, gpoints)

    scalest = if scalest === :stdev
        sqrt(var(x))
    elseif scalest === :iqr
        (quantile(x, 0.75) - quantile(x, 0.25)) / 1.349
    elseif scalest === :minim
        min((quantile(x, 0.75) - quantile(x, 0.25)) / 1.349, sqrt(var(x)))
    else
        error("scalest must be one of :stdev, :iqr or :minim (default)")
    end

    scalest == 0 && error("scale estimate is zero for input data")
    sx = (x .- mean(x)) ./ scalest
    sa = (minx - mean(x)) / scalest
    sb = (maxx - mean(x)) / scalest

    gpoints = range(sa, stop = sb, length = gridsize)
    gcounts = linbin(sx, gpoints)

    hpi = if level == 0
        (24.0*sqrt(pi)/n)^(1/3)
    elseif level == 1
        alpha = (2.0/(3.0*n))^(1/5)*sqrt(2.0)
        psi2hat = bkfe(gcounts, 2, alpha, [sa, sb])
        (6.0 / (-psi2hat * n))^(1 / 3)
    elseif level == 2
        alpha = ((2.0/(5.0*n))^(1/7))*sqrt(2.0)
        psi4hat = bkfe(gcounts, 4, alpha, [sa, sb])
        alpha = (sqrt(2.0 / pi) / (psi4hat * n))^(1 / 5)
        psi2hat = bkfe(gcounts, 2, alpha, [sa, sb])
        (6.0 / (-psi2hat * n))^(1 / 3)
    elseif level == 3
        alpha = ((2.0/(7.0*n))^(1/9))*sqrt(2.0)
        psi6hat = bkfe(gcounts, 6, alpha, [sa, sb])
        alpha = (-3.0 * sqrt(2.0 / pi) / (psi6hat * n))^(1 / 7)
        psi4hat = bkfe(gcounts, 4, alpha, [sa, sb])
        alpha = (sqrt(2.0 / pi) / (psi4hat * n))^(1 / 5)
        psi2hat = bkfe(gcounts, 2, alpha, [sa, sb])
        (6.0 / (-psi2hat * n))^(1 / 3)
    elseif level == 4
        alpha = ((2.0/(9.0*n))^(1/11))*sqrt(2.0)
        psi8hat = bkfe(gcounts, 8, alpha, [sa, sb])
        alpha = (15.0 * sqrt(2.0 / pi) / (psi8hat * n))^(1 / 9)
        psi6hat = bkfe(gcounts, 6, alpha, [sa, sb])
        alpha = (-3.0 * sqrt(2.0 / pi) / (psi6hat * n))^(1 / 7)
        psi4hat = bkfe(gcounts, 4, alpha, [sa, sb])
        alpha = (sqrt(2.0 / pi) / (psi4hat * n))^(1 / 5)
        psi2hat = bkfe(gcounts, 2, alpha, [sa, sb])
        (6.0 / (-psi2hat * n))^(1 / 3)
    elseif level == 5
        alpha = ((2.0 / (11.0 * n))^(1 / 13)) * sqrt(2)
        psi10hat = bkfe(gcounts, 10, alpha, [sa, sb])
        alpha = (-105 * sqrt(2.0 / pi) / (psi10hat * n))^(1 // 11)
        psi8hat = bkfe(gcounts, 8, alpha, [sa, sb])
        alpha = (15.0 * sqrt(2.0 / pi) / (psi8hat * n))^(1 / 9)
        psi6hat = bkfe(gcounts, 6, alpha, [sa, sb])
        alpha = (-3.0 * sqrt(2.0 / pi) / (psi6hat * n))^(1 / 7)
        psi4hat = bkfe(gcounts, 4, alpha, [sa, sb])
        alpha = (sqrt(2.0 / pi) / (psi4hat * n))^(1 / 5)
        psi2hat = bkfe(gcounts, 2, alpha, [sa, sb])
        (6.0 / (-psi2hat * n))^(1 / 3)
    end
    h = scalest * hpi
    k = ceil(Int64, (range_x[2]-range_x[1])/h)
    return k
end

# Compute linear binning of data
function linbin(X, gpoints)
    n, M = length(X), length(gpoints)

    a, b = gpoints[1], gpoints[M]
    gcnts = zeros(M)
    delta = (b - a) / (M - 1)

    @inbounds for i = 1:n
        lxi = ((X[i] - a) / delta) + 1
        li = floor(Int, lxi)
        rem = lxi - li

        if 1 <= li < M
            gcnts[li] += 1 - rem
            gcnts[li + 1] += rem
        end

#=         if !trun
            if lt < 1
                gcnts[1] += 1
            end

            if li >= M
                gcnts[M] += 1
            end
        end =#
    end
    return gcnts
end

# Kernel estimator of the functionals needed for Wand's criterion
function bkfe(gcounts, drv, bandwidth, range_x)
    bandwidth <= 0 && error("'bandwidth' must be strictly positive")

    a, b = range_x
    h = bandwidth
    M = length(gcounts)
    gpoints = range(a, stop = b, length = M)

    ## Set the sample size and bin width

    n = sum(gcounts)
    delta = (b - a) / (M - 1)

    ## Obtain kernel weights

    tau = 4 + drv
    L = min(Int(fld(tau * h, delta)), M)

    lvec = 0:L
    arg = lvec .* delta / h

    # kappam = pdf.(Normal(), arg) ./ h^(drv + 1)
    kappam = @. 1.0/sqrt(2.0*pi) * exp(-0.5*arg ^ 2) / h^(drv + 1)
    hmold0, hmnew = ones(length(arg)), ones(length(arg))
    hmold1 = arg

    if drv >= 2
        for i in (2:drv)
            hmnew = arg .* hmold1 .- (i - 1) .* hmold0
            hmold0 = hmold1       # Compute mth degree Hermite polynomial
            hmold1 = hmnew        # by recurrence.
        end
    end
    kappam = hmnew .* kappam

    ## Now combine weights and counts to obtain estimate
    ## we need P >= 2L+1L, M: L <= M.
    P = nextpow(2, M + L + 1)
    kappam = [kappam; zeros(P - 2 * L - 1); reverse(kappam[2:end])]
    Gcounts = [gcounts; zeros(P - M)]
    kappam = fft(kappam)
    Gcounts = fft(Gcounts)

    return sum(gcounts .* (real(ifft(kappam .* Gcounts)))[1:M]) / (n^2)
end