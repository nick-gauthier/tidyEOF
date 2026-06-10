# Calculate spatial anomaly correlation (averaged over time)

For each time step, compute correlation across spatial locations, then
average across all time steps. Each cell's temporal mean is removed
first so the correlation measures agreement of the anomaly patterns
rather than the shared climatology, which is constant in time and would
otherwise inflate the correlation toward one regardless of skill.

## Usage

``` r
calc_cor_spatial(pred_matrix, obs_matrix)
```

## Arguments

- pred_matrix:

  Predicted values matrix (time x space)

- obs_matrix:

  Observed values matrix (time x space)

## Value

Mean spatial anomaly correlation across time steps
