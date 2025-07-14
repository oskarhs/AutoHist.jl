module MakieExt
isdefined(Base, :get_extension) ? (using Makie) : (using ..Makie)

import AutoHist: AutomaticHistogram

function Makie.convert_arguments(P::Type{<:BarPlot}, h::AutomaticHistogram)
    centers = h.breaks[1:end-1] .+ 0.5 .* diff(h.breaks)
    widths = diff(h.breaks)
    kwargs = (;width = widths, gap = 0, dodge_gap = 0)
    return Makie.to_plotspec(BarPlot, Makie.convert_arguments(P, centers, h.density); kwargs...)
end
Makie.plottype(::AutomaticHistogram) = BarPlot

function Makie.plot!(plot::Hist{<:Tuple{<:AutomaticHistogram}}) # enables hist(h)
    attributes = Makie.Attributes(plot)
    barplot!(plot, attributes, plot[1][])
end

Makie.convert_arguments(P::Type{<:Stairs}, h::AutomaticHistogram) = convert_arguments(P, vcat(h.breaks, h.breaks[end]), vcat(0.0, h.density, 0.0))

function Makie.plot!(plot::StepHist{<:Tuple{<:AutomaticHistogram}}) # enables stephist(h)
    valid_attributes = Makie.Attributes(plot)
    stairs!(plot, valid_attributes, plot[1])
end

end # module