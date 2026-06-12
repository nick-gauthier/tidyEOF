# tidyeof (development version)

* `couple()` and `tune_cca()` gain `method = "pcr"` for principal components
  regression — an OLS alternative to CCA that maps predictor PC amplitudes to
  predictand PC amplitudes. Regularization is the predictor truncation
  (`k_pred`); the CCA-style coupling `k`/`k_cca` is inert for PCR (a warning is
  emitted if passed). CCA-specific diagnostics (`get_canonical_*`) are not
  defined for PCR couplings.
* `patterns()` accepts multi-attribute `stars` objects for combined
  (multivariate) EOF analysis — e.g. joint temperature + precipitation
  downscaling with physically consistent predictions. Requires `scale = TRUE`.
* The whole stack is multivariate-aware: `project_patterns()`,
  `reconstruct()`, `couple()`/`predict()`, `tune_eof()`, `tune_cca()`, and
  plot methods. CV metrics gain per-variable columns plus a pooled score.
* Breaking change: EOF attributes in `patterns$eofs` are now named after the
  variable (previously `"weight"`).
* Cells with non-finite standardized anomalies (climatological sd of 0) are
  now excluded like NA cells.

# tidyeof 0.1.0

Initial release.

## Core Features

* `patterns()` for EOF extraction from `stars` spatiotemporal data, with area weighting, standardization, varimax rotation, and monthly anomalization
* `couple()` and `predict()` for CCA-based statistical downscaling
* `common_patterns()` for joint EOF decomposition across multiple datasets
* `tune_eof()` and `tune_cca()` for cross-validated hyperparameter selection with contiguous temporal folds
* Scree plots with North et al. (1982) error bars and modified Rule N significance testing
* Teleconnection maps with FDR-corrected significance contours
* Full integration with `stars`, `sf`, and the tidyverse

## Architecture

* S3 classes: `patterns`, `coupled_patterns`, `common_patterns`, `cv_folds`
* Plot, print, summary, and predict methods for all major classes
* Automatic IRLBA truncated SVD for large datasets
* Unit preservation throughout the pipeline via the `units` package
