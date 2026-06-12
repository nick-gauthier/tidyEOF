# Changelog

## tidyeof (development version)

- [`patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/patterns.md)
  accepts multi-attribute `stars` objects for combined (multivariate)
  EOF analysis — e.g. joint temperature + precipitation downscaling with
  physically consistent predictions. Requires `scale = TRUE`.
- The whole stack is multivariate-aware:
  [`project_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/project_patterns.md),
  [`reconstruct()`](https://nick-gauthier.github.io/tidyEOF/reference/reconstruct.md),
  [`couple()`](https://nick-gauthier.github.io/tidyEOF/reference/couple.md)/[`predict()`](https://rdrr.io/r/stats/predict.html),
  [`tune_eof()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_eof.md),
  [`tune_cca()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_cca.md),
  and plot methods. CV metrics gain per-variable columns plus a pooled
  score.
- Breaking change: EOF attributes in `patterns$eofs` are now named after
  the variable (previously `"weight"`).
- Cells with non-finite standardized anomalies (climatological sd of 0)
  are now excluded like NA cells.

## tidyeof 0.1.0

Initial release.

### Core Features

- [`patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/patterns.md)
  for EOF extraction from `stars` spatiotemporal data, with area
  weighting, standardization, varimax rotation, and monthly
  anomalization
- [`couple()`](https://nick-gauthier.github.io/tidyEOF/reference/couple.md)
  and [`predict()`](https://rdrr.io/r/stats/predict.html) for CCA-based
  statistical downscaling
- [`common_patterns()`](https://nick-gauthier.github.io/tidyEOF/reference/common_patterns.md)
  for joint EOF decomposition across multiple datasets
- [`tune_eof()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_eof.md)
  and
  [`tune_cca()`](https://nick-gauthier.github.io/tidyEOF/reference/tune_cca.md)
  for cross-validated hyperparameter selection with contiguous temporal
  folds
- Scree plots with North et al. (1982) error bars and modified Rule N
  significance testing
- Teleconnection maps with FDR-corrected significance contours
- Full integration with `stars`, `sf`, and the tidyverse

### Architecture

- S3 classes: `patterns`, `coupled_patterns`, `common_patterns`,
  `cv_folds`
- Plot, print, summary, and predict methods for all major classes
- Automatic IRLBA truncated SVD for large datasets
- Unit preservation throughout the pipeline via the `units` package
