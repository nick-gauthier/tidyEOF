# Restore original field from anomalies and climatology

Reverses the operation of
[`get_anomalies()`](https://nick-gauthier.github.io/tidyEOF/reference/get_anomalies.md),
adding the climatological mean (and optionally multiplying by standard
deviation) back to anomaly fields.

## Usage

``` r
restore_climatology(anomalies, clim, scale = FALSE, monthly = FALSE)
```

## Arguments

- anomalies:

  A stars object containing anomalies (from
  [`get_anomalies()`](https://nick-gauthier.github.io/tidyEOF/reference/get_anomalies.md))

- clim:

  A climatology list with `mean` and `sd` stars objects (from
  [`get_climatology()`](https://nick-gauthier.github.io/tidyEOF/reference/get_climatology.md))

- scale:

  Logical. If TRUE, multiply by standard deviation before adding mean
  (use when anomalies were standardized)

- monthly:

  Logical. If TRUE, restore using monthly climatology. Time steps are
  matched to the climatology by calendar month, so the anomalies may
  start in any month or span partial years.

## Value

A stars object with the original field restored

## Examples

``` r
if (FALSE) { # \dontrun{
clim <- get_climatology(dat)
anom <- get_anomalies(dat, clim)
restored <- restore_climatology(anom, clim)
} # }
```
