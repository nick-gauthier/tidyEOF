# Flatten a stars object to a (dim x space) numeric matrix

Moves \`dim_name\` to the first dimension and flattens the remaining
(spatial) dimensions. For multi-attribute objects, attributes are
concatenated variable-major — all spatial columns for the first
attribute, then all columns for the second, etc. — matching the column
ordering of \[flatten_time_space()\].

## Usage

``` r
flatten_dim_space(x, dim_name)
```

## Arguments

- x:

  A stars object (one or more attributes)

- dim_name:

  Name of the dimension to keep as rows

## Value

A plain numeric matrix with rows = \`dim_name\` and columns = flattened
space (concatenated variable-major for multi-attribute objects). Units
are stripped; no block map is attached.

## Details

Units are dropped: the result is always a plain \`double\` matrix with
no \`"units"\` class. No block map is returned; callers that need
per-variable column ranges should derive them from the spatial shape or
use \[flatten_time_space()\].
