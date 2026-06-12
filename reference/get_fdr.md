# Calculate FDR-corrected significance contours for teleconnections

Computes pixel-wise correlations with FDR correction and returns
significance contour lines as sf polygons.

## Usage

``` r
get_fdr(dat, patterns, fdr = 0.1, amplitudes = NULL)
```

## Arguments

- dat:

  A stars object with a time dimension. For multivariate analyses,
  correlate one variable at a time (e.g. \`get_fdr(dat\["tmean"\],
  pat)\`).

- patterns:

  A patterns object from patterns()

- fdr:

  False discovery rate threshold (default 0.1)

- amplitudes:

  Optional amplitudes tibble (defaults to patterns\$amplitudes)

## Value

An sf object with significance contour polygons for each PC
