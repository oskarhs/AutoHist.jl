# API

## The AutomaticHistogram type
```@docs
AutomaticHistogram
```

An automatic histogram based on regular or irregular partitions can be fitted to the data by calling the `fit` method.
```@docs
fit(::Type{AutomaticHistogram}, x::AbstractVector{<:Real}; rule=:bayes, type=:irregular, kwargs...)
```

Additional methods

```@docs
loglikelihood
logmarginallikelihood
convert
```

## Automatic histogram construction

```@docs
histogram_irregular
histogram_regular
```