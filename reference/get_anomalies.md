# Calculate anomalies from a climatological mean

Calculate anomalies from a climatological mean

## Usage

``` r
get_anomalies(dat, clim = NULL, scale = FALSE, monthly = FALSE)
```

## Arguments

- dat:

  A stars object with dimensions (x, y, time)

- clim:

  Optional climatology from get_climatology(). If NULL, computed
  internally

- scale:

  Logical. If TRUE, divide by standard deviation

- monthly:

  Logical. If TRUE, compute monthly anomalies. Each time step is matched
  to its calendar month, so the data may start in any month, span
  partial years, or cover a single year. Aborts if the data contains a
  month absent from the climatology.

## Value

A stars object with anomalies
