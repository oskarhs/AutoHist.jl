module AutoHistPlotsExt
using Plots

import AutoHist: AutomaticHistogram, AbstractRule, fit
import StatsBase: Histogram

@recipe function f(h::AutomaticHistogram)
    seriestype --> :barbins

    st_map = Dict(
        :histogram => :barbins,
        :barbins => :barbins,
        :barhist => :barbins,
        :bar => :barbins,
        :stephist => :stepbins,
    )
    seriestype := get(st_map, plotattributes[:seriestype], plotattributes[:seriestype])
    if !(plotattributes[:seriestype] in [:barbins, :stepbins])
        throw(ArgumentError("Seriestype :$(plotattributes[:seriestype]) not supported for automatic histograms."))
    end

    x = h.breaks
    y = h.density
    x, y
end

# Fit histogram, plot the result.
@recipe function f(x::AbstractVector{<:Real}, rule::AbstractRule)
    closed = get(plotattributes, :closed, :right)
    support = get(plotattributes, :support, (-Inf, Inf))
    fit(AutomaticHistogram, x, rule; support=support, closed=closed)
end

end