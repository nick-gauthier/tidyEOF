# Project New Data onto Existing Patterns

This function projects new spatial-temporal data onto existing EOF
patterns, returning the corresponding principal component time series.
This is a core function used in pattern-based downscaling and
reconstruction.

## Usage

``` r
project_patterns(patterns, newdata)
```

## Arguments

- patterns:

  A patterns object containing EOFs, climatology, and other metadata

- newdata:

  A stars object with new spatial-temporal data to project

## Value

A tibble with time column and PC amplitude columns

## Details

For multivariate patterns, \`newdata\` must contain the same variables
as the training data (any order); univariate patterns accept any
single-attribute object regardless of name.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get patterns from training data
pat <- patterns(training_data, k = 5)

# Project new data onto these patterns
new_amplitudes <- project_patterns(pat, new_data)
} # }
```
