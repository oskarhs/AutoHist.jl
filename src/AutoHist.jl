module AutoHist

export AutomaticHistogram
export fit, convert, loglikelihood, logmarginallikelihood, minimum, maximum, extrema, peaks, pdf, insupport, length, distance
export histogram_regular, histogram_irregular
export DP, GPDP

using StatsBase, Base.Threads, Distributions
import StatsBase: modes
import Statistics: quantile
import SpecialFunctions: loggamma, logabsbinomial, gamma
import FFTW: fft, ifft
import StatsAPI: fit, loglikelihood
import Distributions: pdf, insupport
if !isdefined(Base, :get_extension)
    using Requires
end

include("AutomaticHistogram.jl")
include("utils.jl")

include(joinpath("regular", "objective_functions.jl"))
include(joinpath("regular", "wand_num_bins.jl"))

include(joinpath("irregular", "algorithms.jl"))
include(joinpath("irregular" ,"greedy_grid.jl"))
include(joinpath("irregular", "objective_functions_irregular.jl"))
include(joinpath("irregular", "dynprog.jl"))
include(joinpath("irregular", "dynprog_greedy.jl"))
include(joinpath("irregular", "compute_bounds.jl"))

include(joinpath("regular", "regular_histogram.jl"))
include(joinpath("irregular", "irregular_histogram.jl"))

function __init__()
    @static if !isdefined(Base, :get_extension) # Latest makie version requires Julia >= v1.10
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include(joinpath("..", "ext", "PlotsExt.jl"))
    end
end

end
