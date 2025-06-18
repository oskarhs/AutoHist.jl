# API

## The AutomaticHistogram type
```@docs
AutomaticHistogram
```

### Fitting an automatic histogram to data
An automatic histogram based on regular or irregular partitions can be fitted to the data by calling the `fit` method.
```@docs
fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}; rule::Symbol=:bayes, type::Symbol=:irregular, kwargs...)
```

Alternatively, an automatic histogram can be fitted to the data by the following methods:
```@docs
histogram_irregular
histogram_regular
```

### Additional methods for AutomaticHist

```@docs
modes(::AutomaticHistogram)
minimum(::AutomaticHistogram)
maximum(::AutomaticHistogram)
extrema(::AutomaticHistogram)
AutoHist.insupport(::AutomaticHistogram, ::Real)
AutoHist.pdf(::AutomaticHistogram, ::Real)
loglikelihood(::AutomaticHistogram)
logmarginallikelihood
convert(::Type{Histogram}, h::AutomaticHistogram)
```


```@index
Pages = ["api.md"]
```