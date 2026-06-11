# Extract calendar month numbers from a time vector

Locale-independent month extraction: \`" \`format()\` dispatches on the
time class, so \`Date\`, \`POSIXct\` (respecting its tzone attribute),
and CF-calendar classes with format methods all work.

## Usage

``` r
month_index(times)
```

## Arguments

- times:

  A vector of time values

## Value

Integer vector of calendar months (1-12)
