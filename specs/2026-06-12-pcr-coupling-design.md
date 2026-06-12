# PCR Coupling Method — Design

**Date:** 2026-06-12
**Status:** Approved
**Branch:** `pcr-coupling` (off `multivariate-eof-cca`)

## Motivation

CCA is currently the only coupling method linking predictor EOF amplitudes to
predictand EOF amplitudes. Principal Components Regression (PCR) is a natural
alternative: because `patterns()` has already reduced each field to orthogonal,
truncated principal components, the coupling is a small linear map from
predictor PC amplitudes (X, `n_times × k_pred`) to predictand PC amplitudes
(Y, `n_times × k_resp`). PCR realises that map as ordinary least squares —
the canonical PCR formulation, where the EOF truncation *is* the
regularization. It is a transparent, dependency-free baseline that, unlike CCA,
does not chase correlation into low-variance directions, so it is often more
robust for short records or strongly collinear predictors.

This is the textbook trio for EOF-based downscaling: PCR (X-variance), PLS
(X–Y covariance), CCA (X–Y correlation). This design adds the first.

## Decisions

| Decision | Choice |
|---|---|
| Estimator | **OLS only.** Plain least squares of predictand PCs on predictor PCs. No ridge/lasso/elastic-net in v1. |
| Regularization knob | The predictor truncation `k_pred`. The CCA-style coupling `k` (canonical-mode count) is **inert** for PCR. |
| Dependency | None. Base R linear algebra (`qr`). No `glmnet`. |
| Diagnostics | Minimal first cut: `couple`/`predict`/`reconstruct`/CV work; the CCA-only accessors abort with a clear message; `print`/`summary` become method-aware. PCR-specific diagnostic maps are deferred. |
| Class | Reuse `coupled_patterns`; store a method-specific `pcr` slot alongside the existing `cca` slot. |
| Multivariate predictands | Work for free — PCR predicts all `k_resp` predictand PCs, and reconstruction/metrics are already multivariate-aware. |

## 1. Architecture and dispatch

All changes are localised to the coupling layer (`R/couple.R`) plus a
`method` pass-through in `R/tune_cca.R`. The downstream stack —
`project_patterns()`, `reconstruct()`, `compute_spatial_metrics()`,
`prep_cv_folds()` — is method-agnostic and unchanged, because PCR is only a
different X→Y map in amplitude space.

`couple()` already takes `method = "cca"` and aborts on anything else; that
stub becomes a real branch. `predict.coupled_patterns()` already has the
matching guard. The `coupled_patterns` S3 class is reused; a PCR object stores
a `pcr` slot and leaves the `cca` slot absent (CCA objects are unchanged).

## 2. `couple(method = "pcr")`

Shared front matter (identical to the CCA path): validate compatibility,
resolve `common_times`, extract `pred_amps`/`resp_amps` filtered to those
times. Then, for PCR:

- **Centering** (honoring the existing `center` argument): when
  `center = TRUE`, `xcenter <- colMeans(pred_amps)`, `ycenter <- colMeans(resp_amps)`,
  and `Xc`/`Yc` are the column-centered matrices. When `center = FALSE`,
  `xcenter <- FALSE`, `ycenter <- FALSE` and no centering is applied (mirroring
  `cancor()`'s sentinel convention so the apply path can branch identically).
- **Fit** via a rank-safe QR: `qrX <- qr(Xc)`; if `qrX$rank < ncol(Xc)`, abort
  with class `tidyeof_rank_deficient` (message points at reducing `k_pred` /
  the `k_pred < n_times` requirement). Otherwise
  `B <- qr.coef(qrX, Yc)` — a `k_pred × k_resp` coefficient matrix (`qr.coef`
  handles the multivariate Y directly).
- **Object**: `method = "pcr"`; `pcr = list(coefficients = B, xcenter, ycenter)`;
  `predictor_patterns`, `response_patterns`, `center` as for CCA; `k <- ncol(pred_amps)`
  (the number of predictor PCs used — for reporting only).
- The CCA-style coupling `k` argument is ignored for PCR. The CCA `k ≤ max_k`
  canonical-mode validation is skipped.

Helper: factor the OLS fit into `fit_pcr(pred_amps, resp_amps, center)` returning
`list(coefficients, xcenter, ycenter)`, keeping `couple()` readable.

## 3. `predict.coupled_patterns()` for PCR

Dispatch on `object$method`. The projection step is unchanged
(`new_amplitudes <- project_patterns(proj_patterns, newdata)`, including the
cross-source `predictor_patterns` override). For PCR, a new internal
`apply_pcr_prediction(new_amplitudes, object$pcr)`:

- `pred_matrix <- as.matrix(new_amplitudes[, -time])`.
- If `xcenter` is not `FALSE`: `pred_matrix <- sweep(pred_matrix, 2, xcenter, "-")`.
- `response <- pred_matrix %*% coefficients`.
- If `ycenter` is not `FALSE`: `response <- sweep(response, 2, ycenter, "+")`.
- Return a tibble with `time` and `PC1..PC{k_resp}` columns (same shape the CCA
  path returns, so `reconstruct()` consumes it unchanged).

`ncol(B) == k_resp == response_patterns$k`, so the predicted amplitude tibble
matches what `reconstruct()` expects. The `k` argument to `predict()` is inert
for PCR (no canonical-mode truncation); `reconstruct = TRUE/FALSE` behaves as
for CCA.

## 4. Accessors, print, summary

- `get_canonical_correlations()`, `get_canonical_variables()`,
  `get_canonical_patterns()`: add a guard at the top — if
  `object$method != "cca"`, abort with class `tidyeof_cca_only` and a message
  explaining these are correlation-specific CCA diagnostics not defined for PCR.
- `print.coupled_patterns()` / `summary.coupled_patterns()`: branch on
  `method`. For PCR, report the method, the number of predictor PCs, and the
  number of response PCs (no canonical correlations). The CCA output is
  unchanged.

## 5. Cross-validation

- `tune_cca()` gains a `method = "cca"` argument, threaded through
  `evaluate_fold()` into `couple(..., method = method)`. `evaluate_fold` and
  the metric flow are otherwise unchanged (PCR's predicted amplitudes go
  through the same `predict()` → `compute_spatial_metrics()` path, so
  per-variable + pooled metrics work for multivariate predictands for free).
- `prep_cv_folds()` is unchanged (it only precomputes patterns).
- For PCR the `k_cca` axis is inert. Left at its default (`k_cca = NULL`),
  `tune_cca` assigns one `k_cca = min(k_pred, k_resp)` per `(k_pred, k_resp)`
  combination, so the grid is effectively `k_pred × k_resp` with no redundant
  rows. Documented: for PCR, leave `k_cca = NULL`; an explicit `k_cca` vector
  would produce duplicate rows that all evaluate identically.

## 6. Error handling

| Condition | Behavior |
|---|---|
| Rank-deficient OLS (`k_pred ≥ n_times` or collinear amplitudes) | `tidyeof_rank_deficient` abort with a hint to reduce predictor modes |
| CCA accessor called on a PCR object | `tidyeof_cca_only` abort |
| Unknown `method` | existing `tidyeof_unsupported_method` abort (extend the allowed set to `c("cca", "pcr")`) |

## 7. Testing

New tests in `tests/testthat/test-pcr.R` (loading `prism` per house convention):

1. **OLS equivalence**: on the orthogonal predictor PCs, the fitted `B` matches
   independent column-wise OLS (and, for a self-coupling `couple(pat, pat)`,
   `B` is ~identity), pinning the fit math.
2. **End-to-end univariate**: `couple(method="pcr") → predict → reconstruct`
   returns a correctly named/unitted stars field with the right dims; amplitude
   prediction (`reconstruct = FALSE`) returns a tibble.
3. **End-to-end multivariate predictand**: univariate predictor → multivariate
   (`scale = TRUE`) response; output has both variables with units restored.
4. **Cross-source predict**: the `predictor_patterns` override path works under
   PCR.
5. **CV**: `tune_cca(method="pcr")` runs, emits the metric columns (incl.
   per-variable for a multivariate response), and `summarize_cv` selects
   sensible `k_pred`/`k_resp`.
6. **Centering**: `center = FALSE` path fits through the origin and round-trips.
7. **Guards**: CCA accessors abort (`tidyeof_cca_only`) on a PCR object;
   rank-deficient fit aborts (`tidyeof_rank_deficient`).
8. **Skill sanity**: on the test data with small `k`, PCR and CCA give similar
   reconstruction skill (loose tolerance) — catches a wiring bug without
   asserting they're identical.

## Out of scope (future work)

- Regularized estimators: ridge (closed-form on orthogonal PCs, no dependency),
  lasso / elastic-net (via `glmnet`), with their `lambda`/`alpha` tuning.
- PCR-specific diagnostic maps (regression coefficient matrix as stars,
  predictor→predictand "pattern" maps, per-response variance explained).
- PLS coupling (`method = "pls"`) — the covariance-criterion sibling.
