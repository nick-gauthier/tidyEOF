# Apply sign flips to a patterns object

Flips the sign of EOFs, amplitudes, and projection matrix according to
the supplied sign vector. This keeps all components synchronized.
Supports both single-attribute (univariate) and multi-attribute
(multivariate) EOF stars objects.

## Usage

``` r
apply_sign_flips(patterns, signs)
```

## Arguments

- patterns:

  A patterns object

- signs:

  Named numeric vector of +1/-1 values (from compute_eof_signs)

## Value

The patterns object with signs applied
