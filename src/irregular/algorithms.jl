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
breaks: [0.00012207031249977798, 0.20299585314764493, 0.30080999844316125, 0.36239668251811596, 0.4783245584239131, 0.5725159575973732, 0.6594618645267211, 0.750030517578125, 0.8297309322633605, 0.9130540930706521, 1.0000000000000002]
density: [0.009862770955956847, 0.0663639674261088, 0.1502328303595117, 0.30260936719244613, 0.5864544365976523, 0.9385761107406808, 1.4040555093789138, 1.9764287121852526, 2.650673218857813, 3.5054229130654364]
counts: [9, 32, 46, 175, 276, 408, 636, 788, 1105, 1525]
type: irregular
closed: right
a: 5.0

julia> h = fit(AutomaticHistogram, x; alg = GPDP(greedy=false, max_cand=30))
AutomaticHistogram
breaks: [0.00012207031249977798, 0.20299585314764493, 0.3044327445652174, 0.40586963598278986, 0.49643828903419385, 0.5906296882076539, 0.681198341259058, 0.7608987559442936, 0.8442219167515851, 0.9239223314368207, 1.0000000000000002]
density: [0.009862770955956847, 0.06796890780356953, 0.18418118149879492, 0.37161779107231424, 0.6479696688274588, 1.0334368412688513, 1.502626437612729, 2.0727804974905486, 2.766099169806118, 3.564842182949155]
counts: [9, 34, 93, 168, 305, 468, 599, 864, 1103, 1357]
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
