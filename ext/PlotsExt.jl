module PlotsExt
using Plots

import AutoHist: AutomaticHistogram
import StatsBase: Histogram

# Plot the corresponding StatsBase.Histogram
@recipe f(h::AutomaticHistogram) = convert(Histogram, h)

end