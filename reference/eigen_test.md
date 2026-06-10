# Test EOF significance using modified Rule N

Tests whether the k-th eigenvalue is significantly different from noise
using a modified Rule N approach based on Tracy-Widom distribution.

## Usage

``` r
eigen_test(lambdas, k, M, n, p = 0.05, total_var = NULL)
```

## Arguments

- lambdas:

  Vector of eigenvalues from PCA

- k:

  Index of eigenvalue to test

- M:

  Number of spatial points (grid cells)

- n:

  Number of time steps

- p:

  Significance level (default 0.05)

- total_var:

  Total variance of the data (sum of all eigenvalues). Required when
  \`lambdas\` contains only the leading modes (e.g., from IRLBA); the
  remaining noise variance is then \`total_var\` minus the eigenvalues
  above \`k\`. If NULL (default), all eigenvalues must be present.

## Value

Logical, TRUE if eigenvalue is significant at level p

## Details

Based on the gamma approximation to the Tracy-Widom distribution
described in Cheng & Wallace (1993) and implemented following Overland &
Preisendorfer (1982). Constants (shape = 46.4, scale factor = 0.186,
location = 9.85) derive from fitting the gamma CDF to the Tracy-Widom
Type 1 distribution.
