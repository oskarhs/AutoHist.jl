module AutoHist

export AutomaticHistogram
export fit, convert, loglikelihood, logmarginallikelihood, minimum, maximum, extrema, peaks, pdf, cdf, insupport, length, distance
export histogram_regular, histogram_irregular
export DP
export RRH, Knuth, AIC, BIC, BR, MDL, NML_R, L2CV_R, KLCV_R, Sturges, FD, Scott, Wand

using StatsBase, Base.Threads, Distributions
import StatsBase: modes
import Statistics: quantile
import SpecialFunctions: loggamma, logabsbinomial, gamma
import FFTW: fft, ifft
import StatsAPI: fit, loglikelihood
import Distributions: pdf, insupport

abstract type AbstractRule end

abstract type AbstractRegularRule <: AbstractRule end

abstract type AbstractIrregularRule <: AbstractRule end

include("AutomaticHistogram.jl")
include("utils.jl")

include(joinpath("regular", "objective_functions.jl"))
include(joinpath("regular", "wand_num_bins.jl"))

include(joinpath("irregular", "algorithms.jl"))
include(joinpath("irregular" ,"greedy_grid.jl"))
include(joinpath("irregular", "objective_functions_irregular.jl"))
include(joinpath("irregular", "dynprog.jl"))
include(joinpath("irregular", "compute_bounds.jl"))

include(joinpath("regular", "regular_histogram.jl"))
include(joinpath("regular", "regular_rules.jl"))
include(joinpath("irregular", "irregular_histogram.jl"))

function __init__()
    return nothing
end

end
