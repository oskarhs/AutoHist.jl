module AutoHist

export AutomaticHistogram
export fit
export convert, loglikelihood, logmarginallikelihood, minimum, maximum, extrema, modes
export histogram_regular, histogram_irregular

import StatsBase: modes
using StatsBase, Base.Threads
import Statistics: quantile
import SpecialFunctions: loggamma, logabsbinomial, gamma
import FFTW: fft, ifft
import StatsAPI: fit, loglikelihood
if !isdefined(Base, :get_extension)
    using Requires
end

include("utils.jl")

include(joinpath("regular", "objective_functions.jl"))
include(joinpath("regular", "wand_num_bins.jl"))

include(joinpath("irregular" ,"greedy_grid.jl"))
include(joinpath("irregular", "dynamic_algorithm.jl"))

include(joinpath("regular", "regular_histogram.jl"))
include(joinpath("irregular", "irregular_histogram.jl"))

include("AutomaticHistogram.jl")

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80" include(joinpath("..", "ext", "PlotsExt.jl"))
    end
end

end
