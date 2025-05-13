module AutoHist

export histogram_regular, histogram_irregular

using StatsBase, Base.Threads, LoopVectorization
import Statistics: quantile
import SpecialFunctions: loggamma, logabsbinomial

include("utils.jl")

include(joinpath("regular", "objective_functions.jl"))

include(joinpath("irregular" ,"greedy_grid.jl"))
include(joinpath("irregular", "dynamic_algorithm.jl"))

include(joinpath("regular", "regular_histogram.jl"))
include(joinpath("irregular", "irregular_histogram.jl"))

end
