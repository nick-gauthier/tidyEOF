# Calculate teleconnections between a patterns object and another layer

Computes pixel-wise correlations between a spatiotemporal field and each
PC amplitude time series. Checks for overlapping time steps between the
raster field and the PC amplitudes.

## Usage

``` r
get_correlation(dat, patterns, amplitudes = NULL)
```

## Arguments

- dat:

  A stars object with a time dimension. For multivariate analyses,
  correlate one variable at a time (e.g.
  \`get_correlation(dat\["tmean"\], pat)\`).

- patterns:

  A patterns object from patterns()

- amplitudes:

  Optional amplitudes tibble (defaults to patterns\$amplitudes)

## Value

A stars object with correlation values for each PC
