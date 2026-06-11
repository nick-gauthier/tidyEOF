# Get EOFs and PCs from spatiotemporal data

This function performs Empirical Orthogonal Function (EOF) analysis on
spatial-temporal data. For large datasets, it automatically uses IRLBA
(Implicitly Restarted Lanczos Bidiagonalization Algorithm) for efficient
computation when the irlba package is available.

## Usage

``` r
patterns(
  dat,
  k = 4,
  scale = FALSE,
  rotate = FALSE,
  monthly = FALSE,
  weight = TRUE,
  irlba_threshold = 5e+05
)
```

## Arguments

- dat:

  A \`stars\` object containing spatial and temporal dimensions

- k:

  The number of PC/EOF modes to retain

- scale:

  Logical, whether to scale before PCA

- rotate:

  Logical, whether to apply Varimax rotation. Rotation follows the
  standard REOF convention (Hannachi et al. 2007): varimax operates on
  sqrt(eigenvalue)-scaled EOFs without Kaiser row-normalization
  (\`normalize = FALSE\`), the stored patterns are the rotated loadings
  (unit norm, not mutually orthogonal), and amplitudes remain
  uncorrelated with sd = sqrt(rotated eigenvalue). Requires k \> 1.

- monthly:

  Logical, whether to use monthly climatology

- weight:

  Logical, whether to apply area weighting

- irlba_threshold:

  Minimum number of data elements to trigger IRLBA usage (default:
  500000). Set to Inf to always use base prcomp().

## Value

A \`patterns\` object containing EOFs, amplitudes, and metadata

## Details

The eigenvalue table includes North et al. (1982) sampling-error bars
(\`low\`/\`hi\`), computed with a relative error of sqrt(2/n) where n is
the number of time steps. This assumes temporally independent samples:
for autocorrelated data (e.g., monthly anomalies) the effective sample
size is smaller and the bars are too narrow, so modes that appear
well-separated may not be.
