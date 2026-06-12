# Evaluate EOF reconstruction for a single fold via speckled holdout

Hides a random scatter of grid cells in the held-out data, estimates
mode amplitudes from the visible cells by least squares, and scores the
prediction on the hidden cells. This makes reconstruction skill a
genuine out-of-sample quantity, so over-fitting with too many EOFs is
penalised (Bro et al. 2008). Masks depend only on the fold and
replicate, not on \`k\`, so every \`k\` is scored on the same hidden
cells.

## Usage

``` r
evaluate_eof_fold(
  fold,
  k,
  metrics,
  hidden_fraction = 0.2,
  n_reps = 5,
  seed = 1L
)
```

## Arguments

- fold:

  A fold list containing train_patterns and test_data

- k:

  Number of EOFs to use

- metrics:

  Metrics to compute

- hidden_fraction:

  Fraction of valid grid cells to hide per replicate

- n_reps:

  Number of random hidden-cell masks to average over

- seed:

  Base seed; combined with the fold id so masks are reproducible

## Value

Tibble with fold_id and metric values. For multivariate fields the
metrics include pooled scores plus per-variable scores (e.g.
\`rmse_tmean\`).
