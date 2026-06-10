# Test EOF significance using modified Rule N

Tests whether the k-th eigenvalue is significantly different from noise
using a modified Rule N approach based on the Tracy-Widom distribution.

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

Follows the Rule N significance framework for EOFs of Overland &
Preisendorfer (1982), comparing each eigenvalue against the noise null
given by the Tracy-Widom Type-1 law for the largest eigenvalue of a
white-noise covariance matrix (Johnstone 2001). The Tracy-Widom CDF is
evaluated through the gamma approximation of Chiani (2014); its
constants (shape = 46.4, scale factor = 0.186, location = 9.85) come
from fitting a gamma CDF to the Tracy-Widom Type-1 distribution. Because
the null assumes independent noise, the test tends to be liberal for
spatially correlated geophysical fields.

## References

Overland, J.E. & Preisendorfer, R.W. (1982). A significance test for
principal components applied to a cyclone climatology. *Monthly Weather
Review*, 110(1), 1-4.

Johnstone, I.M. (2001). On the distribution of the largest eigenvalue in
principal components analysis. *Annals of Statistics*, 29(2), 295-327.

Chiani, M. (2014). Distribution of the largest eigenvalue for real
Wishart and Gaussian random matrices and a simple approximation for the
Tracy-Widom distribution. *Journal of Multivariate Analysis*, 129,
69-81.
