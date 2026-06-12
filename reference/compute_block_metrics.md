# Compute per-variable and pooled metrics over block-structured matrices

The plain metric name is always the pooled score. For a single block it
is the metric itself, preserving univariate behavior exactly. For
multiple blocks, RMSE is normalized by each block's observed standard
deviation and RMS-combined (raw pooling across different units is
meaningless), and correlations (already unitless) are averaged.
Per-variable values get suffixed names (e.g. rmse_tmean) only when there
is more than one block.

## Usage

``` r
compute_block_metrics(
  pred_matrix,
  obs_matrix,
  block_map,
  metrics = c("rmse", "cor_spatial", "cor_temporal")
)
```

## Arguments

- pred_matrix:

  Predicted values matrix (time x space)

- obs_matrix:

  Observed values matrix (time x space)

- block_map:

  Named list of column indices per variable

- metrics:

  Character vector of metric names

## Value

Named list of metric values

## Details

Pooled RMSE assumes each block's observed standard deviation is \> 0; a
constant (zero-variance) observed variable yields a non-finite pooled
score.
