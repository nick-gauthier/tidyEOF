# Smart PCA Selection with IRLBA Support

Automatically selects between base \`prcomp()\` and \`prcomp_irlba()\`
based on data size. For large datasets, IRLBA provides significant
computational savings when only the first few components are needed.

## Usage

``` r
perform_pca_smart(
  x,
  k = NULL,
  center = TRUE,
  scale. = FALSE,
  size_threshold,
  ...
)
```

## Arguments

- x:

  A numeric matrix for PCA computation

- k:

  Number of components to compute

- center:

  Logical, whether to center the data

- scale.:

  Logical, whether to scale the data

- size_threshold:

  Minimum number of elements to trigger IRLBA

- ...:

  Additional arguments passed to the PCA function

## Value

A PCA result object compatible with \`prcomp()\` output
