# Tidy PCA standard deviations into a tibble

Replaces \`broom::tidy(pca, matrix = "pcs")\` with a dependency-free
version. Returns a tibble with columns PC, std.dev, percent, cumulative.

## Usage

``` r
tidy_pca_sdev(pca_result, total_var = NULL)
```

## Arguments

- pca_result:

  A prcomp (or prcomp_irlba) result

- total_var:

  Total variance of the data (sum of all eigenvalues). Defaults to
  \`sum(pca_result\$sdev^2)\`, which is only correct when the full
  spectrum is present; pass \`pca_result\$totalvar\` for truncated
  (IRLBA) results.

## Value

Tibble with PC, std.dev, percent, cumulative
