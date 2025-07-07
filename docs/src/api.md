# API

## The AutomaticHistogram type
```@docs
AutomaticHistogram
```

## Fitting an automatic histogram to data
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
peaks(::AutomaticHistogram)
minimum(::AutomaticHistogram)
maximum(::AutomaticHistogram)
extrema(::AutomaticHistogram)
insupport(::AutomaticHistogram, ::Real)
pdf(::AutomaticHistogram, ::Real)
length(::AutomaticHistogram)
loglikelihood(::AutomaticHistogram)
logmarginallikelihood
convert(::Type{Histogram}, h::AutomaticHistogram)
distance
```

## Index

```@index
Pages = ["api.md"]
```