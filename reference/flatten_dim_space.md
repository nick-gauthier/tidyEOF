# Flatten a stars object to a (dim x space) numeric matrix

Moves \`dim_name\` to the first dimension and flattens the remaining
(spatial) dimensions, matching the column ordering used by
\[flatten_time_space()\]. Units are dropped.

## Usage

``` r
flatten_dim_space(x, dim_name)
```

## Arguments

- x:

  A single-attribute stars object

- dim_name:

  Name of the dimension to keep as rows

## Value

A matrix with rows = \`dim_name\`, columns = flattened space
