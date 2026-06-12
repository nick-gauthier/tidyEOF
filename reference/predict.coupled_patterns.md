# Predict Method for Coupled Patterns

This function makes predictions using a coupled_patterns object created
by couple(). It applies the learned CCA relationship to new predictor
data to predict response patterns.

## Usage

``` r
# S3 method for class 'coupled_patterns'
predict(
  object,
  newdata,
  k = NULL,
  reconstruct = TRUE,
  predictor_patterns = NULL,
  ...
)
```

## Arguments

- object:

  A coupled_patterns object from couple()

- newdata:

  New predictor data (stars object) for making predictions

- k:

  Number of CCA modes to use for prediction (CCA only; if NULL, uses all
  available modes). Ignored for \`method = "pcr"\`.

- reconstruct:

  Logical, whether to reconstruct the full spatial field (default: TRUE)

- predictor_patterns:

  Optional patterns object to use instead of the one stored in the
  coupled object. Useful for cross-source prediction with common EOFs:
  the override patterns share the same EOF space but carry a different
  climatology.

- ...:

  Additional arguments (currently unused)

## Value

If reconstruct=TRUE, returns a stars object with reconstructed spatial
fields. If reconstruct=FALSE, returns a tibble with predicted
amplitudes.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create coupled patterns
coupled <- couple(pred_patterns, resp_patterns, k = 3)

# Make predictions on new data
predictions <- predict(coupled, new_predictor_data)

# Just get predicted amplitudes without spatial reconstruction
amplitudes <- predict(coupled, new_predictor_data, reconstruct = FALSE)

# Cross-source prediction with common EOFs
cpat <- common_patterns(list(era = era, phyda = phyda), k = 5)
coupled <- couple(cpat$era, fine_patterns, k = 3)
predict(coupled, phyda_new, predictor_patterns = cpat$phyda)
} # }
```
