abstract type AbstractAlgorithm end

# Corresponds to using the dynamic programming algorithm of Kanazawa (1988)
struct DP <: AbstractAlgorithm
    greedy::Bool
    gr_maxbins::Union{Int, Symbol}
end

"""
    DP(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)

Dynamic programming algorithm for constructing an irregular histogram.

# Keyword arguments
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most max(500, n^(1/3))+1 candidate cutpoints (including edges).

!!! note
    This algorithm can be quite slow for large datasets when the `greedy` keyword is set to `false`.

# Examples
```jldoctest 
julia> using AutoHist

julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x; alg = DP(greedy=true, gr_maxbins=200))
AutomaticHistogram
breaks: [0.00012207031249977798, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9239223314368207, 1.0000000000000002]
density: [0.00662683597412854, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.7513848134529346, 3.564842182949155]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1147, 1357]
type: irregular
closed: right
a: 5.0

julia> h = fit(AutomaticHistogram, x; alg = DP(greedy=false))
AutomaticHistogram
breaks: [0.00012207031249977798, 0.17763663029325183, 0.29718725232110504, 0.4022468898607337, 0.4928155429121377, 0.5797614498414855, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9202995853147645, 1.0000000000000002]
density: [0.00662683597412854, 0.057821970706400425, 0.17596277991076312, 0.36279353706969375, 0.6214544825215076, 0.9730458529384184, 1.4481767793920146, 2.0440057561776532, 2.733509595364622, 3.545742066060367]
counts: [5, 34, 92, 164, 270, 423, 656, 852, 1090, 1414]
type: irregular
closed: right
a: 5.0
```
"""
function DP(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default)
    if typeof(gr_maxbins) == Symbol && gr_maxbins != :default
        throw(ArgumentError("Keyword argument gr_maxbins must be an integer or :default."))
    elseif typeof(gr_maxbins) == Int && gr_maxbins ≤ 0
        throw(DomainError("Keyword argument gr_maxbins must be positive."))
    end
    return DP(greedy, gr_maxbins)
end

# Corresponds to using the greedy pruned dynamic programming algorithm of Simensen et al. (2025)
struct GPDP <: AbstractAlgorithm
    greedy::Bool
    gr_maxbins::Union{Int, Symbol}
    max_cand::Union{Int, Symbol}
end

"""
    GPDP(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default, max_cand::Union{Int, Symbol}=:default)

Greedy pruned dynamic programming algorithm for constructing an irregular histogram.  Only supported for rules `:bayes`, `:penr`, `:penb`, `:pena`, and `:nml`.

# Keyword arguments
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of [Rozenholc et al. (2010)](https://doi.org/10.1016/j.csda.2010.04.021) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most max(1000, n^(1/3))+1 cutpoint candidates (including edges).
- `max_cand`: Number of cutpoints to consider as possible candidates for the maximum in each step of the greedy algorithm. Default value is `15` for the greedy segment neighborhood algorithm.

!!! note
    This algorithm can be quite slow for large datasets when the `greedy` keyword is set to `false`.

# Examples
```jldoctest 
julia> using AutoHist

julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x; alg = GPDP())
AutomaticHistogram
breaks: [0.00012207031249977798, 0.1703911380491395, 0.30080999844316125, 0.3841331592504529, 0.4928155429121377, 0.5978751804517662, 0.6775755951370018, 0.7645215020663496, 0.8478446628736412, 0.9202995853147645, 1.0000000000000002]
density: [0.006866313123570044, 0.056150710478631356, 0.16405599040198, 0.3429389655637976, 0.6552110576392909, 1.0413586147484761, 1.5038798827835582, 2.101555238803439, 2.764093657532795, 3.545742066060367]
counts: [5, 36, 68, 186, 344, 415, 654, 876, 1002, 1414]
type: irregular
closed: right
a: 5.0

julia> h = fit(AutomaticHistogram, x; alg = GPDP(greedy=false, max_cand=30))
AutomaticHistogram
breaks: [0.00012207031249977798, 0.1160499462182971, 0.2573370449784873, 0.3768876670063406, 0.4747018123018569, 0.5761387037194293, 0.6667073567708333, 0.7572760098222373, 0.8405991706295289, 0.9202995853147645, 1.0000000000000002]
density: [0.0027226100721400408, 0.030696131624918386, 0.13302868095600692, 0.3135247855550627, 0.5859998905466929, 0.9650488727485436, 1.4481767793920146, 2.0440057561776532, 2.733509595364622, 3.545742066060367]
counts: [1, 21, 79, 153, 297, 437, 656, 852, 1090, 1414]
type: irregular
closed: right
a: 5.0
```
"""
function GPDP(; greedy::Bool=true, gr_maxbins::Union{Int, Symbol}=:default, max_cand::Union{Int, Symbol}=:default)
    if typeof(gr_maxbins) == Symbol && gr_maxbins != :default
        throw(ArgumentError("Keyword argument gr_maxbins must be a positive integer or :default."))
    elseif typeof(gr_maxbins) == Int && gr_maxbins ≤ 0
        throw(DomainError("Keyword argument gr_maxbins must be positive."))
    end
    if typeof(max_cand) == Symbol && max_cand != :default
        throw(ArgumentError("Keyword argument max_cand must be a positive integer or :default."))
    elseif typeof(max_cand) == Int && max_cand ≤ 0
        throw(DomainError("Keyword argument max_cand must be positive."))
    end
    return GPDP(greedy, gr_maxbins, max_cand)
end
