abstract type AbstractAlgorithm end

struct OptPart <: AbstractAlgorithm
    greedy::Bool
    gr_maxbins::Union{Int, Symbol}
end

"""
    OptPart(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)

The optimal partitioning algorithm for constructing an irregular histogram.

# Keyword arguments
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most ``\\max \\{ 3000, n^{1/2} \\}+1`` candidate cutpoints (including edges).

!!! note
    This algorithm can be quite slow for large datasets when the `greedy` keyword is set to `false`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = LinRange(eps(), 1.0-eps(), 1000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x, L2CV_I(alg = OptPart(greedy=true, gr_maxbins=100)))
AutomaticHistogram
breaks: [0.0001220703125, 0.16676839192708331, 0.2977047874813988, 0.4048345656622024, 0.5119643438430059, 0.6071908133370535, 0.6786106654575893, 0.7619338262648809, 0.8452569870721726, 0.9285801478794643, 1.0]
density: [0.006000732511292883, 0.05346107146424569, 0.17735498311154516, 0.3920478574044684, 0.7035858869490909, 1.0641298986692698, 1.5001831278232214, 2.0762534489073383, 2.7963413502624808, 3.5984392626052997]
counts: [1, 7, 19, 42, 67, 76, 125, 173, 233, 257]
type: irregular
closed: right
a: NaN

julia> h = fit(AutomaticHistogram, x, L2CV_I(alg = OptPart(greedy=false)))
AutomaticHistogram
breaks: [0.0001220703125, 0.16676839192708331, 0.2977047874813988, 0.4048345656622024, 0.5119643438430059, 0.6071908133370535, 0.6786106654575893, 0.7619338262648809, 0.8452569870721726, 0.9285801478794643, 1.0]
density: [0.006000732511292883, 0.05346107146424569, 0.17735498311154516, 0.3920478574044684, 0.7035858869490909, 1.0641298986692698, 1.5001831278232214, 2.0762534489073383, 2.7963413502624808, 3.5984392626052997]
counts: [1, 7, 19, 42, 67, 76, 125, 173, 233, 257]
type: irregular
closed: right
a: NaN
```
"""
function OptPart(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)
    if typeof(gr_maxbins) == Symbol && gr_maxbins != :default
        throw(ArgumentError("Keyword argument gr_maxbins must be an integer or :default."))
    elseif typeof(gr_maxbins) == Int && gr_maxbins ≤ 0
        throw(DomainError("Keyword argument gr_maxbins must be positive."))
    end
    return OptPart(greedy, gr_maxbins)
end

struct SegNeig <: AbstractAlgorithm
    greedy::Bool
    gr_maxbins::Union{Int, Symbol}
end

"""
    SegNeig(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)

The segment neighbourhood algorithm for constructing an irregular histogram.

# Keyword arguments
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most ``\\max\\{ 500, n^{1/3} \\}+1`` candidate cutpoints (including edges).

!!! note
    This algorithm can be quite slow for large datasets when the `greedy` keyword is set to `false`.

# Examples
```jldoctest; setup = :(using AutoHist)
julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x, RIH(alg = SegNeig(greedy=true, gr_maxbins=200)))
AutomaticHistogram
breaks: [0.0001220703125, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9239223314368207, 1.0]
density: [0.006626835974128547, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.7513848134529346, 3.5648421829491657]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1147, 1357]
type: irregular
closed: right
a: 5.0

julia> h = fit(AutomaticHistogram, x, RIH(alg = SegNeig(greedy=false)))
AutomaticHistogram
breaks: [0.0001220703125, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9202995853147645, 1.0]
density: [0.006626835974128547, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.733509595364622, 3.545742066060377]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1090, 1414]
type: irregular
closed: right
a: 5.0
```
"""
function SegNeig(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)
    if typeof(gr_maxbins) == Symbol && gr_maxbins != :default
        throw(ArgumentError("Keyword argument gr_maxbins must be an integer or :default."))
    elseif typeof(gr_maxbins) == Int && gr_maxbins ≤ 0
        throw(DomainError("Keyword argument gr_maxbins must be positive."))
    end
    return SegNeig(greedy, gr_maxbins)
end