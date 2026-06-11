# Apply or remove a monthly climatology by calendar month index

Shared engine for \[get_anomalies()\] and \[restore_climatology()\] with
\`monthly = TRUE\`. Because the climatology's month dimension is always
1:12 in calendar order, each time step's climatology is looked up by
direct indexing with its month number – no ordering convention to
maintain, no complete-years requirement, and identical code for raster
and geometry data.

## Usage

``` r
apply_monthly_climatology(
  dat,
  clim,
  scale,
  direction = c("anomalize", "restore")
)
```

## Arguments

- dat:

  A stars object with a time dimension (data or anomalies)

- clim:

  Climatology list from \[get_climatology()\] with \`monthly = TRUE\`

- scale:

  Logical, whether to divide/multiply by the climatological sd

- direction:

  "anomalize" (subtract climatology) or "restore" (add it back)

## Value

A stars object with the same structure as \`dat\`, units dropped
