# Couple Pattern Relationships Using CCA

This function couples predictor and response patterns using Canonical
Correlation Analysis (CCA) as the primary method. CCA finds linear
combinations of predictor and response patterns that maximize
correlation between them.

## Usage

``` r
couple(
  predictor_patterns,
  response_patterns,
  k = NULL,
  method = "cca",
  center = TRUE,
  validate = TRUE
)
```

## Arguments

- predictor_patterns:

  A patterns object containing predictor patterns (e.g., from
  patterns())

- response_patterns:

  A patterns object containing response patterns (e.g., from patterns())

- k:

  Number of CCA modes to retain. If NULL, uses min(ncol(predictor),
  ncol(response))

- method:

  Coupling method. Currently only "cca" is supported

- center:

  Logical, whether to center the amplitudes before CCA (default: TRUE).
  Centering is the statistically standard choice and makes retaining all
  modes equivalent to multivariate regression with an intercept. It is a
  no-op when the amplitudes are already zero-mean over the coupled
  period (the usual case), but is essential when the predictor and
  response patterns were fit on different periods and then filtered to a
  common one, which leaves the common-period amplitudes with a nonzero
  mean.

- validate:

  Logical, whether to validate input patterns compatibility

## Value

A coupled_patterns object containing:

- cca:

  The CCA results from cancor()

- predictor_patterns:

  The original predictor patterns

- response_patterns:

  The original response patterns

- k:

  Number of CCA modes retained

- method:

  Coupling method used

## Examples

``` r
if (FALSE) { # \dontrun{
# Get patterns from your data
pred_patterns <- patterns(predictor_data, k = 5)
resp_patterns <- patterns(response_data, k = 5)

# Couple the patterns
coupled <- couple(pred_patterns, resp_patterns, k = 3)

# Make predictions
predictions <- predict(coupled, new_predictor_data)
} # }
```
