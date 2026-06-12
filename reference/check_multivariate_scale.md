# Check that multivariate input is standardized

PCA is variance-driven, so joint EOFs of variables with different units
require per-pixel standardization to contribute comparably.

## Usage

``` r
check_multivariate_scale(dat, scale, call = rlang::caller_env())
```

## Arguments

- dat:

  A stars object

- scale:

  The scale argument passed to the caller

- call:

  Calling environment for error messages
