module AutoHistMakieExt
using Makie

import AutoHist: AutomaticHistogram, AbstractRule, fit

function Makie.convert_arguments(P::Type{<:BarPlot}, h::AutomaticHistogram)
    centers = 0.5 * (h.breaks[1:end-1] + h.breaks[2:end] )
    widths = diff(h.breaks)
    kwargs = (;width = widths, gap = 0, dodge_gap = 0)
    return Makie.to_plotspec(BarPlot, Makie.convert_arguments(P, centers, h.density); kwargs...)
end
Makie.plottype(::AutomaticHistogram) = BarPlot

function Makie.plot!(plot::Hist{<:Tuple{<:AutomaticHistogram}}) # enables hist(h)
    attributes = Makie.shared_attributes(plot, BarPlot)
    barplot!(plot, attributes, plot[1][])
end

Makie.convert_arguments(P::Type{<:Stairs}, h::AutomaticHistogram) = convert_arguments(P, vcat(h.breaks, h.breaks[end]), vcat(0.0, h.density, 0.0))

function Makie.plot!(plot::StepHist{<:Tuple{<:AutomaticHistogram}}) # enables stephist(h)
    attributes = Makie.shared_attributes(plot, BarPlot)
    stairs!(plot, attributes, plot[1][])
end

function Makie.convert_arguments(P::Type{<:BarPlot}, x::AbstractVector{<:Real}, rule::AbstractRule) # enables e.g. `plot(x, AIC())`
    h = fit(AutomaticHistogram, x, rule)
    return Makie.convert_arguments(P, h)
end

Makie.plottype(::AbstractVector{<:Real}, ::AbstractRule) = BarPlot

function Makie.convert_arguments(P::Type{<:Hist}, x::AbstractVector{<:Real}, rule::AbstractRule) # enables e.g. `hist(x, AIC())`
    h = fit(AutomaticHistogram, x, rule)
    return Makie.convert_arguments(P, h)
end

function Makie.plot!(plot::Hist{<:Tuple{<:AbstractVector{<:Real}, <:AbstractRule}}) # enables e.g. `hist(x, AIC())`
    attributes = Makie.shared_attributes(plot, BarPlot)
    barplot!(plot, attributes, plot[1][])
end

function Makie.convert_arguments(P::Type{<:StepHist}, x::AbstractVector{<:Real}, rule::AbstractRule)
    h = fit(AutomaticHistogram, x, rule)
    Makie.convert_arguments(P, h)
end


function Makie.plot!(plot::StepHist{<:Tuple{<:AbstractVector{<:Real}, <:AbstractRule}}) # enables e.g. `stephist(x, AIC())`
    attributes = Makie.shared_attributes(plot, BarPlot)
    stairs!(plot, attributes, plot[1][])
end

end # module