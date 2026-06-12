# Extract the EOF loading matrix (concatenated variable-space x PC) from a patterns object

Stacks every variable's loadings into the concatenated variable-major
layout used by \[flatten_time_space()\], so rows align with
\`patterns\$valid_pixels\` and \`patterns\$block_map\`.

## Usage

``` r
eof_loading_matrix(patterns)
```

## Arguments

- patterns:

  A patterns object

## Value

A numeric matrix with prod(spatial) \* n_vars rows and k columns
