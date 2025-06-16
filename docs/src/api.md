# API

## The AutomaticHistogram type
```@docs
AutomaticHistogram
```

An automatic histogram based on regular or irregular partitions can be fitted to the data by calling the `fit` method.
```@docs
fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}; rule::Symbol=:bayes, type::Symbol=:irregular, kwargs...)
```

Alternatively, an automatic histogram can be fitted to the data by the following methods:
```@docs
histogram_irregular
histogram_regular
```

Additional methods

```@docs
loglikelihood(::AutomaticHistogram)
logmarginallikelihood
convert(::Histogram, h::AutomaticHistogram)
minimum(::AutomaticHistogram)
maximum(::AutomaticHistogram)
extrema(::AutomaticHistogram)
```
