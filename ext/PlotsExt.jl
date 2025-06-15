module PlotsExt

isdefined(Base, :get_extension) ? (using Plots) : (using ..Plots)

import AutoHist: AutomaticHistogram
import RecipesBase: @recipe
import StatsBase: Histogram

# Plot the corresponding StatsBase.Histogram
@recipe f(::Type{AutomaticHistogram}, h::AutomaticHistogram) = convert(Histogram, h)

end