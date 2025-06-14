module PlotsExt

isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

using AutoHist, RecipesBase, StatsBase

# Import package based on version
#isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

# Plot the corresponding StatsBase.Histogram
@recipe f(::Type{AutomaticHistogram}, h::AutomaticHistogram) = convert(Histogram, h)

end