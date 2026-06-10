# Matrix-centric rotation helper that keeps all components synchronized

Standard REOF convention (Hannachi et al. 2007): varimax operates on the
sqrt(eigenvalue)-scaled loadings, and the rotated scaled loadings ARE
the rotated patterns. They are stored unit-norm with the variance
carried by the scores, mirroring the unrotated convention (score sd =
sqrt(eigenvalue)); rotated patterns are not mutually orthogonal, but
scores stay uncorrelated.

## Usage

``` r
rotate_pca_components(loadings_matrix, scores_matrix, sdev_vector)
```
