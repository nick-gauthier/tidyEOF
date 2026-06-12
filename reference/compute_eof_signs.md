# Compute sign vector for consistent EOF orientation

Determines the sign (+1 or -1) needed for each EOF component so that the
dominant loading is positive across ALL variable blocks. Supports both
single-attribute (univariate) and multi-attribute (multivariate) stars
objects, and both raster (x, y, PC) and sf geometry (geometry, PC)
layouts.

## Usage

``` r
compute_eof_signs(eofs)
```

## Arguments

- eofs:

  A stars object with EOF spatial patterns (must have a PC dimension as
  the last dimension). May contain one or more attributes (variables).

## Value

Named numeric vector of +1/-1 values, one per PC
