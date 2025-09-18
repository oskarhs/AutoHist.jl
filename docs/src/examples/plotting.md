# Plotting tutorial

The following markdown page gives a brief, example-driven introduction to plotting histograms with AutoHist.jl. All the plotting functionality is provided through extensions for [Makie.jl](https://docs.makie.org/stable/) and [Plots.jl](https://docs.juliaplots.org/stable/). Both extensions provide recipes for plotting [`AutomaticHistogram`](@ref) objects, and additionally lets the user utilize the automatic bin selection rules from AutoHist when calling the respective built-in histogram plotting functions from the two libraries.

## Plots.jl
In order to use the AutoHist.jl plotting functionality with Plots.jl, we will need to load both libraries first.
```@example Plots; continued = true
using AutoHist; import Plots
```

We start by showing how to plot the result from calling `fit` or `autohist`. The Plots extension provides a type recipe for objects of the `AutomaticHistogram` type, so that `Plots.plot(h)` will draw the histogram `h`, as illustrated in the code snippet below.
```@example Plots
import Random
x = Random.randexp(Random.Xoshiro(1812), 10^5)

h1 = autohist(x, AIC())
h2 = autohist(x, BIC())

# Set axis limits
xlims = [-0.5, 6.0]
ylims = [-0.05, 1.0]

# Plot a filled histogram
p1 = Plots.plot(h1, color=:red, alpha=0.4, label="", title="AIC", xlims=xlims, ylims=ylims)

# Plot a stephist
p2 = Plots.stephist(h2, color=:green, label="", title="BIC", xlims=xlims, ylims=ylims)

# Display the two plots side-by-side
Plots.plot(p1, p2, layout=(1, 2), size=(670, 320))
```
As illustrated here, the call to `Plots.plot` also accepts additional keyword arguments which allow us to customize the look of the histogram. `AutomaticHistogram` objects are compatible with most histogram-like seriestypes[^1], so the above call to `Plots.plot` is equivalent to calling `Plots.histogram` with the same input. By default, the plotted histogram is a filled bar chart. If one desires, a non-filled histogram can be drawn by using the `stephist` seriestype, as we did above.

[^1]:
    In particular, the seriestypes `histogram`, `barbins`, `barhist` and `bar` will draw the histogram as a filled bar chart, while `stephist` draws a non-filled histogram. An attempt to plot an `AutomaticHistogram` with another [seriestype](https://docs.juliaplots.org/latest/generated/attributes_series/) throws an error.

!!! note
    Some of the arguments normally supported by the `histogram` and `stephist` seriestypes are deprecated when used together with `AutomaticHistogram`. In particular, the `weights` keyword is not supported, as `fit`/`autohist` currently only handles equally weighted samples. The `bins` keyword is ignored, as the histogram partition has already been chosen by `autohist`. Additionally, the `normalize` kwarg is deprecated, so that e.g. `Plots.histogram(h, normalize=norm)` will always draw a histogram on the density scale no matter the value of `norm`.

Automatic histograms can also be plotted by accessing Plots.jl functions directly. The general syntax for fitting a histogram in this manner is `Plots.plot(x, rule)`, where `rule` is any of the bin selection rules provided by the AutoHist.jl package, see the [methods page](../methods.md). Note that the above syntax supports any histogram-like seriestype, so the above inline code snippet is equivalent to `Plots.histogram(x, rule)`.
When using Plots.jl functions directly, one can pass additional keyword arguments directly to the plotting function. In the code snippet below, we show how to directly plot a right-closed histogram with support constrained to be positive.
```@example Plots
# Plot a filled histogram
p3 = Plots.plot(x, BayesBlocks(), color=:black, alpha=0.4, label="",
                title="BayesBlocks", xlims=xlims, ylims=ylims)

# Plot a stephist, supported on [0, maximum(x)] with left-closed intervals
p4 = Plots.stephist(x, RRH(), support=(0.0, Inf), closed=:left, color=:blue,
                    label="", title="RRH", xlims=xlims, ylims=ylims)

# Display the two plots
Plots.plot(p3, p4, layout=(1, 2), size=(670, 320))
```

## Makie.jl
To showcase the Makie.jl extension, we start by loading the required packages:
```@example Makie
using AutoHist; import Makie, CairoMakie;
```

Quickly plotting a pre-fitted `AutomaticHistogram` object can be done easily via `Makie.plot(h)`.
`AutomaticHistogram` objects also supports other histogram-like recipes such as [`hist`](https://docs.makie.org/stable/reference/plots/hist), [`barplot`](https://docs.makie.org/stable/reference/plots/barplot) and [`stephist`](https://docs.makie.org/stable/reference/plots/stephist), so e.g. `Makie.hist(h)` can also be used to draw a fitted histogram.
Below, we show a slightly more involved example where we plot two histograms fitted using different rules side-by-side, using a `stephist` for the last panel.
```@example Makie
import Random
x = Random.randexp(Random.Xoshiro(1812), 10^5)

h1 = autohist(x, AIC())
h2 = autohist(x, BIC())

# Specify axis lims, ticks
limits = ((-0.5, 6.0), (-0.05, 1.0))
yticks = 0.0:0.25:1.0

# Make a figure with a 2x1 Layout
fig = Makie.Figure(size=(670, 320))
ax1 = Makie.Axis(fig[1, 1], title="AIC", xlabel="x",
                 ylabel="Density", limits=limits, yticks = yticks)
ax2 = Makie.Axis(fig[1, 2], title="BIC", xlabel="x", limits=limits, yticks=yticks)

# Draw a filled AIC histogram in the left panel
Makie.plot!(ax1, h1, alpha=0.4, color="red")

# Draw a BIC stephist in the right panel
Makie.stephist!(ax2, h2, color="green")

fig # Display figure
```
As shown in the above example, we can pass additional plot attributes to the Makie plotting functions in order to customize our histogram plots. We refer the interested reader to the [Makie documentation](https://docs.makie.org/stable/reference/plots/hist) for further details on the supported attributes.

!!! note
    Some of the arguments normally supported by the `barplot`, `hist` and `stephist` seriestypes are deprecated when used together with `AutomaticHistogram`. In particular, the `weights` keyword is not supported, as `fit`/`autohist` currently only handle equally weigthed samples. The `bins` keyword is ignored, as the histogram partition has already been chosen by `autohist`. Additionally, the `normalization` kwarg is deprecated, so that e.g. `Makie.hist(h, normalization=norm)` will always draw a histogram on the density scale no matter the value of `norm`.

We can also plot automatically selected histograms directly without needing to construct an `AutomaticHistogram` object directly first. This is achieved by calling e.g. `Makie.plot(x, rule)`, where `rule` is any of the bin selection rules provided as part of AutoHist.jl, see the [methods page](../methods.md). Alternatively, we can use another histogram-like series type like `stephist` instead. In the following code snippet, we illustrate the use of these methods on the simulated exponential dataset considered previously:
```@example Makie
# Make a figure with a 2x1 Layout
fig = Makie.Figure(size=(670, 320))
ax1 = Makie.Axis(fig[1, 1], title="BayesBlocks", xlabel="x",
                 ylabel="Density", limits=limits, yticks=yticks)
ax2 = Makie.Axis(fig[1, 2], title="RRH", xlabel="x", limits=limits, yticks=yticks)

# Draw a filled BayesBlocks histogram in the left panel
Makie.plot!(ax1, x, BayesBlocks(), alpha=0.4, color="black")

# Draw a RRH stephist in the right panel
Makie.stephist!(ax2, x, RRH(), color="blue")

fig # Display figure
```

!!! note
    Currently, the Makie.jl plotting functions do not support passing additional keyword arguments to `fit`/`autohist`.
    In order to draw an automatic histogram using Makie with non-default values for the `support` and `closed` kwargs, it is therefore necessary to manually fit an `AutomaticHistogram` object first, and then plot this using the syntax outlined previously.