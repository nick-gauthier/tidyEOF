# Cross-validate EOF truncation for a single field

Evaluates reconstruction skill for different numbers of EOFs using
k-fold cross-validation with a speckled holdout. For each held-out fold
a random scatter of grid cells is hidden, mode amplitudes are estimated
from the visible cells, and the hidden cells are predicted. Because the
hidden cells are not used to estimate the amplitudes, prediction error
is genuinely out-of-sample and stops improving once \`k\` exceeds the
field's effective rank – so the RMSE-minimising \`k\` is a meaningful
estimate of how many modes the data support (Bro et al. 2008).

## Usage

``` r
tune_eof(
  data,
  k = 1:10,
  kfolds = 5,
  max_k = max(k),
  metrics = c("rmse", "cor_spatial", "cor_temporal"),
  scale = FALSE,
  monthly = FALSE,
  weight = TRUE,
  hidden_fraction = 0.2,
  n_reps = 5,
  seed = 1L
)
```

## Arguments

- data:

  A stars object with spatial-temporal data

- k:

  Vector of EOF counts to evaluate (default 1:10)

- kfolds:

  Number of cross-validation folds (default 5)

- max_k:

  Maximum EOFs to compute per fold (default max(k))

- metrics:

  Character vector of metrics to compute. Options: "rmse",
  "cor_spatial", "cor_temporal" (default: all three). Metrics are
  computed on the hidden cells only.

- scale:

  Logical, whether to scale data before EOF extraction (default FALSE)

- monthly:

  Logical, whether to compute monthly climatology (default FALSE)

- weight:

  Logical, whether to apply area weighting (default TRUE)

- hidden_fraction:

  Fraction of grid cells to hide in each held-out fold (default 0.2).
  Hidden cells are predicted from the visible ones.

- n_reps:

  Number of random hidden-cell masks to average over per fold (default
  5). More replicates give smoother, more stable estimates.

- seed:

  Base random seed for hidden-cell masks (default 1). Masks depend only
  on the fold and replicate, not on \`k\`, so all \`k\` are compared on
  the same hidden cells. The global RNG is left undisturbed.

## Value

A tibble with columns: k, fold, and one column per metric.

## Details

Naively projecting the full held-out field onto the EOFs and scoring the
reconstruction (the approach used before tidyeof 0.1.0) does NOT work:
the projection is least-squares optimal for the very data being scored,
so error decreases monotonically with \`k\` and the "best" \`k\` is
always the largest.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find optimal k for precipitation field
results <- tune_eof(precip_data, k = 1:15, kfolds = 5)
summary <- summarize_eof_cv(results, metric = "rmse")

# Plot reconstruction skill vs k
library(ggplot2)
results %>%
  group_by(k) %>%
  summarize(rmse = mean(rmse)) %>%
  ggplot(aes(k, rmse)) + geom_line() + geom_point()
} # }
```
