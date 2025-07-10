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
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of Rozenholc et al. (2010) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`. The algorithm can be quite slow for large datasets when this keyword is set to `false`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most max(500, n^(1/3))+1 candidate cutpoints (including edges).

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

Greedy pruned dynamic programming algorithm for constructing an irregular histogram.

# Keyword arguments
- `greedy`: Boolean indicating whether or not the greedy cutpoint selection strategy of Rozenholc et al. (2010) should be used to select a smaller number of candidate cutpoints prior to running the dynamic programming algorithm. Defaults to `true`. The algorithm can be quite slow for large datasets when this keyword is set to `false`.
- `gr_maxbins`: Number of candidate cutpoints chosen by the greedy algorithm. Supplying `gr_maxbins=:default` results in the selection of at most max(1000, n^(1/3))+1 cutpoint candidates (including edges).
- `max_cand`: Number of cutpoints to consider as possible candidates for the maximum in each step of the greedy algorithm. Default value is `15` for the greedy segment neighborhood algorithm and `45` for the greedy optimal partitioning algorithm.

# Examples
```jldoctest 
julia> using AutoHist

julia> x = LinRange(eps(), 1.0-eps(), 5000) .^(1.0/4.0);

julia> h = fit(AutomaticHistogram, x; alg = GPDP())
AutomaticHistogram
breaks: [0.00012207031249977798, 0.19575036090353262, 0.30080999844316125, 0.3841331592504529, 0.4928155429121377, 0.5978751804517662, 0.6667073567708333, 0.7645215020663496, 0.8514674089956975, 0.9239223314368207, 1.0000000000000002]
density: [0.009169728520235164, 0.0637578259981318, 0.16405599040198, 0.3429389655637976, 0.6552110576392909, 1.0140467041841466, 1.4717081233990028, 2.11514331109008, 2.7999421894184446, 3.564842182949155]
counts: [8, 33, 68, 186, 344, 349, 720, 920, 1015, 1357]
type: irregular
closed: right
a: 5.0

julia> h = fit(AutomaticHistogram, x; alg = GPDP(greedy=false, max_cand=30))
AutomaticHistogram
breaks: [0.00012207031249977798, 0.13778642295063406, 0.2682052833446558, 0.3841331592504529, 0.4928155429121377, 0.5870069420855979, 0.6775755951370018, 0.7572760098222373, 0.8405991706295289, 0.9202995853147645, 1.0000000000000002]
density: [0.0039018380946941765, 0.03776684797317198, 0.14404855308285908, 0.3429389655637976, 0.6352423794006015, 1.0157883332636115, 1.4800644245378491, 2.0440057561776532, 2.733509595364622, 3.545742066060367]
counts: [2, 24, 83, 186, 299, 460, 590, 852, 1090, 1414]
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
