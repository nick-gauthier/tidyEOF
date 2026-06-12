# PCR Coupling Method Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Principal Components Regression (`method = "pcr"`) as an OLS-based alternative to CCA in the `couple()` downscaling stack, fully integrated with `predict()`, the diagnostics, and cross-validation.

**Architecture:** PCR is a different linear map in amplitude space; all work is localized to `R/couple.R` (a `fit_pcr()` helper, an `apply_pcr_prediction()` applier, method dispatch in `couple()`/`predict()`, method-aware `print`/`summary`, and CCA-only guards on the canonical accessors) plus a `method` pass-through in `R/tune_cca.R`. The patterns/projection/reconstruction/metrics layers are unchanged. Regularization is the predictor truncation `k_pred`; the CCA-style coupling `k` is inert for PCR.

**Tech Stack:** R package; base-R linear algebra (`qr`/`qr.coef`), tidyverse, testthat. No new dependency. Test data: `inst/testdata/prism_test.RDS` (single attribute `tmean`, °C, 51×51 grid, 36 monthly Dates 2017-01-01..2019-12-01). The helper `make_multivar()` (in `tests/testthat/helper-multivariate.R`, auto-loaded by testthat for all test files) builds a two-attribute `tmean`+`ppt` object.

**Spec:** `specs/2026-06-12-pcr-coupling-design.md`. **Branch:** `pcr-coupling`.

**Verified facts (do not re-derive):**
- `couple()` (`R/couple.R:44-100`) currently: validates, extracts `pred_amps`/`resp_amps` via `extract_amplitudes_matrix(..., common_times)`, then `cancor()`, and builds a `coupled_patterns` list with a `cca` slot. A `method != "cca"` abort stub exists at lines 56-59.
- `extract_amplitudes_matrix(x, times)` returns a numeric matrix (rows = times, cols = `PC1..PCk`), time-filtered and row-sorted. Column names are `PC1`, `PC2`, … for both predictor and response.
- `predict.coupled_patterns()` (`R/couple.R:197-257`) aborts unless `method == "cca"` (lines 206-209), defaults/validates `k`, projects newdata, calls `apply_cca_prediction()`, then `reconstruct()` (unless `reconstruct = FALSE`).
- `apply_cca_prediction()` (`R/couple.R:271-325`) shows the centering convention: `cca_result$xcenter`/`ycenter` are either a named numeric vector or the sentinel `FALSE`; the applier re-aligns by `colnames` and `sweep`s. Mirror this for PCR.
- Accessors `get_canonical_variables` (`:341`), `get_canonical_patterns` (`:414`), `get_canonical_correlations` (`:476`) all read `object$cca$...`.
- `print.coupled_patterns` (`:106`) reads `x$cca$cor`; `summary.coupled_patterns` (`:120`) calls `get_canonical_correlations`.
- `tune_cca()` (`R/tune_cca.R:221`) calls `evaluate_fold()` (`:320`) which calls `couple(pred, resp, k = k_cca, validate = FALSE)` then `predict()` then `compute_spatial_metrics()`. `evaluate_fold`'s metric loop iterates `names(metric_values)`, so per-variable + pooled columns flow through automatically.
- NAMESPACE is unchanged by this work (new functions are `@keywords internal`; `couple`/`predict`/`tune_cca` are already exported).

**Run tests with:** `cd "/Users/nick/UF Dropbox/Nick Gauthier/Projects/tidyeof" && Rscript -e 'devtools::test(filter = "pcr")'` for the new file, and `Rscript -e 'devtools::test()'` for the full suite. (Quote the path — it has spaces.)

---

### Task 1: `fit_pcr()` + `couple(method = "pcr")` + method-aware print/summary

**Files:**
- Modify: `R/couple.R` — add `fit_pcr()`, branch `couple()` on `method`, make `print`/`summary` method-aware
- Create: `tests/testthat/test-pcr.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-pcr.R`:

```r
prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()

test_that("couple(method = 'pcr') builds a PCR coupled object", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")

  expect_s3_class(cpl, "coupled_patterns")
  expect_equal(cpl$method, "pcr")
  expect_true(is.matrix(cpl$pcr$coefficients))
  expect_equal(dim(cpl$pcr$coefficients), c(4, 3))  # k_pred x k_resp
  expect_null(cpl$cca)
  expect_equal(cpl$k, 4)  # predictor PCs used (reporting only)
})

test_that("PCR OLS coefficients match a direct lm fit", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")

  X <- as.matrix(pred$amplitudes[, -1])
  Y <- as.matrix(resp$amplitudes[, -1])
  B_lm <- coef(lm(Y ~ X))[-1, , drop = FALSE]  # drop intercept row
  expect_equal(unname(cpl$pcr$coefficients), unname(B_lm), tolerance = 1e-8)
})

test_that("PCR aborts on rank-deficient predictors", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  pred_bad <- pred
  pred_bad$amplitudes$PC4 <- pred_bad$amplitudes$PC2  # collinear -> rank 3 < 4
  expect_error(couple(pred_bad, resp, method = "pcr"),
               class = "tidyeof_rank_deficient")
})

test_that("couple rejects an unknown method", {
  pred <- patterns(prism, k = 3, weight = FALSE)
  resp <- patterns(prism, k = 2, weight = FALSE)
  expect_error(couple(pred, resp, method = "banana"),
               class = "tidyeof_unsupported_method")
})

test_that("print and summary do not error for a PCR coupled object", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")
  expect_no_error(print(cpl))
  expect_no_error(summary(cpl))
})

test_that("CCA coupling is unchanged", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "cca", k = 2)
  expect_equal(cpl$method, "cca")
  expect_false(is.null(cpl$cca))
  expect_equal(cpl$k, 2)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "pcr")'`
Expected: FAIL — `couple(method = "pcr")` hits the `tidyeof_unsupported_method` abort stub, so no `pcr` slot exists.

- [ ] **Step 3: Add the `fit_pcr()` helper**

In `R/couple.R`, immediately BEFORE the `couple <- function(...)` definition (currently line 44), insert:

```r
#' Fit the OLS coefficient matrix for PCR coupling
#'
#' Ordinary least squares of the predictand PC amplitudes on the predictor PC
#' amplitudes. The EOF truncation already performed the dimension reduction, so
#' this is plain (optionally centered) multivariate regression. The fit is
#' rank-safe: a rank-deficient predictor matrix (e.g. `k_pred >= n_times` or
#' collinear amplitudes) aborts rather than returning `NA` coefficients.
#'
#' @param pred_amps Predictor amplitude matrix (time x k_pred)
#' @param resp_amps Predictand amplitude matrix (time x k_resp)
#' @param center Logical; center both sides before fitting (default TRUE)
#' @return List with `coefficients` (k_pred x k_resp), `xcenter`, `ycenter`
#'   (each a named numeric vector when centered, or `FALSE`)
#' @keywords internal
fit_pcr <- function(pred_amps, resp_amps, center = TRUE) {
  if (isTRUE(center)) {
    xcenter <- colMeans(pred_amps)
    ycenter <- colMeans(resp_amps)
    Xc <- sweep(pred_amps, 2, xcenter, "-")
    Yc <- sweep(resp_amps, 2, ycenter, "-")
  } else {
    xcenter <- FALSE
    ycenter <- FALSE
    Xc <- pred_amps
    Yc <- resp_amps
  }

  qrX <- qr(Xc)
  if (qrX$rank < ncol(Xc)) {
    cli::cli_abort(
      c(
        "Predictor amplitudes are rank-deficient ({qrX$rank} < {ncol(Xc)} columns); OLS is not identifiable.",
        "i" = "Reduce the number of predictor EOFs so it is below the number of time steps and free of collinearity."
      ),
      class = "tidyeof_rank_deficient"
    )
  }

  coefficients <- qr.coef(qrX, Yc)
  list(coefficients = coefficients, xcenter = xcenter, ycenter = ycenter)
}
```

- [ ] **Step 4: Branch `couple()` on `method`**

Replace the body of `couple()` (`R/couple.R:44-100`) with the following (keep the existing roxygen block above it; the only `@param method` doc change happens in Task 5):

```r
couple <- function(predictor_patterns, response_patterns, k = NULL,
                  method = "cca", center = TRUE, validate = TRUE) {

  if (!method %in% c("cca", "pcr")) {
    cli::cli_abort(
      "Unsupported coupling method {.val {method}}. Use {.val cca} or {.val pcr}.",
      class = "tidyeof_unsupported_method"
    )
  }

  # Validate inputs and get common times
  common_times <- if (validate) {
    validate_patterns_compatibility(predictor_patterns, response_patterns)
  } else {
    pred_times <- get_times(predictor_patterns)
    resp_times <- get_times(response_patterns)
    pred_times[pred_times %in% resp_times]
  }

  # Extract amplitude matrices filtered to common times
  pred_amps <- extract_amplitudes_matrix(predictor_patterns, common_times)
  resp_amps <- extract_amplitudes_matrix(response_patterns, common_times)

  # Notify user if filtering occurred
  pred_n <- length(get_times(predictor_patterns))
  resp_n <- length(get_times(response_patterns))
  if (length(common_times) < pred_n || length(common_times) < resp_n) {
    cli::cli_inform("Filtered to {length(common_times)} common time steps (predictor had {pred_n}, response had {resp_n}).")
  }

  coupled <- list(
    predictor_patterns = predictor_patterns,
    response_patterns = response_patterns,
    method = method,
    center = center
  )

  if (method == "cca") {
    # Determine and validate k (number of canonical modes)
    if (is.null(k)) {
      k <- min(ncol(pred_amps), ncol(resp_amps))
    }
    max_k <- min(ncol(pred_amps), ncol(resp_amps))
    if (k > max_k) {
      warning("k = ", k, " exceeds maximum possible (", max_k, "). Setting k = ", max_k)
      k <- max_k
    }
    coupled$cca <- cancor(pred_amps, resp_amps, xcenter = center, ycenter = center)
    coupled$k <- k
  } else {  # method == "pcr"
    # The coupling-k is inert for PCR; regularization is the predictor
    # truncation (k_pred). k is stored only for reporting.
    coupled$pcr <- fit_pcr(pred_amps, resp_amps, center = center)
    coupled$k <- ncol(pred_amps)
  }

  class(coupled) <- "coupled_patterns"
  return(coupled)
}
```

- [ ] **Step 5: Make `print`/`summary` method-aware**

Replace `print.coupled_patterns` (`R/couple.R:106-114`) with:

```r
print.coupled_patterns <- function(x, ...) {
  cli::cli_h1("Coupled Patterns Object")
  cli::cli_text("Method: {.field {x$method}}")
  if (x$method == "cca") {
    cli::cli_text("CCA modes retained: {.field {x$k}}")
    cli::cli_text("Canonical correlations: {.val {round(x$cca$cor[1:x$k], 3)}}")
  } else {
    cli::cli_text("Predictor PCs used: {.field {x$k}}")
  }
  cli::cli_text("Predictor patterns: {.field {ncol(extract_amplitudes_matrix(x$predictor_patterns))}} PCs")
  cli::cli_text("Response patterns: {.field {ncol(extract_amplitudes_matrix(x$response_patterns))}} PCs")
  invisible(x)
}
```

Replace `summary.coupled_patterns` (`R/couple.R:120-128`) with:

```r
summary.coupled_patterns <- function(object, ...) {
  cli::cli_h1("Coupled Patterns Summary")
  cli::cli_text("Method: {.field {object$method}}")
  cli::cli_text("Centered: {.field {object$center}}")
  if (object$method == "cca") {
    cli::cli_text("CCA modes retained: {.field {object$k}}")
    cli::cli_h2("Canonical Correlations")
    print(get_canonical_correlations(object))
  } else {
    cli::cli_text("Predictor PCs used: {.field {object$k}}")
    cli::cli_text("Response PCs predicted: {.field {ncol(object$pcr$coefficients)}}")
  }
  invisible(object)
}
```

- [ ] **Step 6: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "pcr")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (CCA path unchanged; `test-refactored_coupling.R`, `test-tune_cca.R` green)

- [ ] **Step 7: Commit**

```bash
git add R/couple.R tests/testthat/test-pcr.R
git commit -m "PCR coupling: fit_pcr and couple(method = 'pcr')"
```

---

### Task 2: PCR prediction path

**Files:**
- Modify: `R/couple.R` — add `apply_pcr_prediction()`, dispatch `predict.coupled_patterns()` on method
- Test: `tests/testthat/test-pcr.R`

- [ ] **Step 1: Write the failing tests**

APPEND to `tests/testthat/test-pcr.R`:

```r
test_that("PCR predict + reconstruct returns a named field (univariate)", {
  set.seed(1)
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))
  pred <- patterns(filter(coarse, time <= as.Date("2018-12-01")), k = 4, weight = FALSE)
  resp <- patterns(filter(prism, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")
  test <- filter(coarse, time > as.Date("2018-12-01"))

  rec <- predict(cpl, test)
  expect_s3_class(rec, "stars")
  expect_named(rec, "tmean")
  expect_equal(unname(dim(rec)), c(51, 51, 12))

  amps <- predict(cpl, test, reconstruct = FALSE)
  expect_s3_class(amps, "tbl_df")
  expect_named(amps, c("time", "PC1", "PC2", "PC3"))
})

test_that("PCR downscales to a multivariate response", {
  set.seed(2)
  prism_mv <- make_multivar(prism)
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))
  pred <- patterns(filter(coarse, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  resp <- patterns(filter(prism_mv, time <= as.Date("2018-12-01")),
                   k = 3, scale = TRUE, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")

  rec <- predict(cpl, filter(coarse, time > as.Date("2018-12-01")))
  expect_named(rec, c("tmean", "ppt"))
  expect_equal(units(rec[["ppt"]]), units(prism_mv[["ppt"]]))
})

test_that("PCR honors center = FALSE", {
  pred <- patterns(prism, k = 3, weight = FALSE)
  resp <- patterns(prism, k = 2, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr", center = FALSE)
  expect_identical(cpl$pcr$xcenter, FALSE)
  expect_identical(cpl$pcr$ycenter, FALSE)
  amps <- predict(cpl, prism, reconstruct = FALSE)
  expect_named(amps, c("time", "PC1", "PC2"))
})

test_that("PCR supports the cross-source predictor_patterns override", {
  cpat <- common_patterns(list(a = prism, b = setNames(prism, "tmean")),
                          k = 4, scale = TRUE, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(cpat$a, resp, method = "pcr")
  pr <- predict(cpl, prism, predictor_patterns = cpat$b, reconstruct = FALSE)
  expect_named(pr, c("time", "PC1", "PC2", "PC3"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "pcr")'`
Expected: FAIL — `predict()` aborts with `tidyeof_unsupported_method` for a PCR object (the guard at `R/couple.R:206-209` only allows `"cca"`).

- [ ] **Step 3: Add the `apply_pcr_prediction()` applier**

In `R/couple.R`, immediately AFTER `apply_cca_prediction()` (it ends at line 325), insert:

```r
#' Apply PCR Prediction Transform
#'
#' Internal function that maps predictor PC amplitudes to predicted predictand
#' PC amplitudes via the fitted OLS coefficient matrix.
#'
#' @param new_amplitudes Tibble with `time` and predictor PC amplitudes
#' @param pcr The `pcr` slot of a coupled object: `coefficients`, `xcenter`, `ycenter`
#' @return Tibble with `time` and predicted response PC amplitudes
#' @keywords internal
apply_pcr_prediction <- function(new_amplitudes, pcr) {
  new_times <- new_amplitudes$time

  pred_matrix <- new_amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  # Apply training centering if the fit was centered (mirrors apply_cca_prediction)
  if (!identical(pcr$xcenter, FALSE)) {
    xcenter <- pcr$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(pred_matrix))) {
      xcenter <- xcenter[colnames(pred_matrix)]
    }
    pred_matrix <- sweep(pred_matrix, 2, xcenter, "-")
  }

  response_amplitudes <- pred_matrix %*% pcr$coefficients

  if (!identical(pcr$ycenter, FALSE)) {
    ycenter <- pcr$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(response_amplitudes))) {
      ycenter <- ycenter[colnames(response_amplitudes)]
    }
    response_amplitudes <- sweep(response_amplitudes, 2, ycenter, "+")
  }

  n_response_pcs <- ncol(response_amplitudes)
  pc_names <- paste0("PC", 1:n_response_pcs)

  response_amplitudes %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(pc_names) %>%
    mutate(time = new_times, .before = 1)
}
```

- [ ] **Step 4: Dispatch `predict.coupled_patterns()` on method**

Replace `predict.coupled_patterns()` (`R/couple.R:197-257`) with (keeps the CCA path numerically identical; routes PCR to the new applier; `k` is inert for PCR):

```r
predict.coupled_patterns <- function(object, newdata, k = NULL, reconstruct = TRUE,
                                   predictor_patterns = NULL, ...) {

  if (!inherits(object, "coupled_patterns")) {
    cli::cli_abort("object must be a coupled_patterns object from couple()",
                   class = "tidyeof_invalid_input")
  }

  if (!object$method %in% c("cca", "pcr")) {
    cli::cli_abort("Unsupported coupling method {.val {object$method}}.",
                   class = "tidyeof_unsupported_method")
  }

  # Use override patterns if provided, otherwise use stored patterns
  proj_patterns <- predictor_patterns %||% object$predictor_patterns

  # Validate compatibility if overriding
  if (!is.null(predictor_patterns)) {
    if (!identical(dim(predictor_patterns$proj_matrix),
                   dim(object$predictor_patterns$proj_matrix))) {
      cli::cli_abort(
        "{.arg predictor_patterns} must share the same EOF space as the training predictor patterns.",
        class = "tidyeof_incompatible_patterns"
      )
    }
  }

  # Project new data onto predictor patterns to get PC amplitudes
  new_amplitudes <- project_patterns(proj_patterns, newdata)

  predicted_amplitudes <- if (object$method == "cca") {
    if (is.null(k)) {
      k <- object$k
    }
    if (k > object$k) {
      warning("Requested k (", k, ") exceeds available modes (", object$k, "). Using k = ", object$k)
      k <- object$k
    }
    apply_cca_prediction(new_amplitudes = new_amplitudes, cca_result = object$cca, k = k)
  } else {
    # k is inert for PCR
    apply_pcr_prediction(new_amplitudes = new_amplitudes, pcr = object$pcr)
  }

  if (!reconstruct) {
    return(predicted_amplitudes)
  }

  reconstruct(target_patterns = object$response_patterns,
              amplitudes = predicted_amplitudes)
}
```

- [ ] **Step 5: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "pcr")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-refactored_coupling.R`, `test-multivariate.R` exercise the CCA predict path and must stay green)

- [ ] **Step 6: Commit**

```bash
git add R/couple.R tests/testthat/test-pcr.R
git commit -m "PCR coupling: prediction path"
```

---

### Task 3: CCA-only guards on the canonical accessors

**Files:**
- Modify: `R/couple.R` — guard `get_canonical_variables`, `get_canonical_patterns`, `get_canonical_correlations`
- Test: `tests/testthat/test-pcr.R`

- [ ] **Step 1: Write the failing tests**

APPEND to `tests/testthat/test-pcr.R`:

```r
test_that("CCA accessors abort on a PCR coupled object", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")
  expect_error(get_canonical_correlations(cpl), class = "tidyeof_cca_only")
  expect_error(get_canonical_patterns(cpl, type = "response"), class = "tidyeof_cca_only")
  expect_error(get_canonical_variables(cpl, resp, type = "response"), class = "tidyeof_cca_only")
})

test_that("CCA accessors still work for a CCA coupled object", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "cca", k = 2)
  expect_s3_class(get_canonical_correlations(cpl), "data.frame")
  expect_s3_class(get_canonical_patterns(cpl, type = "response"), "stars")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "pcr")'`
Expected: FAIL — the accessors read `object$cca$...` on a PCR object, producing a subscript/`NULL` error rather than the typed `tidyeof_cca_only` class.

- [ ] **Step 3: Add a shared guard and call it from each accessor**

In `R/couple.R`, immediately BEFORE `get_canonical_variables` (currently line 341), insert the helper:

```r
#' Guard a CCA-only accessor against non-CCA coupled objects
#' @keywords internal
check_cca_method <- function(object, fn, call = rlang::caller_env()) {
  if (object$method != "cca") {
    cli::cli_abort(
      c(
        "{.fn {fn}} is a CCA-specific diagnostic, not defined for a {.val {object$method}} coupling.",
        "i" = "Canonical correlations, variates, and patterns exist only for {.code method = \"cca\"}."
      ),
      class = "tidyeof_cca_only",
      call = call
    )
  }
}
```

Then add a call as the FIRST statement inside each accessor body:

- In `get_canonical_variables` (after the `function(...) {` line at `:341`):
  ```r
    check_cca_method(object, "get_canonical_variables")
  ```
- In `get_canonical_patterns` (after the `function(...) {` line at `:414`):
  ```r
    check_cca_method(object, "get_canonical_patterns")
  ```
- In `get_canonical_correlations` (after the `function(...) {` line at `:476`):
  ```r
    check_cca_method(object, "get_canonical_correlations")
  ```

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "pcr")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-refactored_coupling.R` exercises CCA accessors)

- [ ] **Step 5: Commit**

```bash
git add R/couple.R tests/testthat/test-pcr.R
git commit -m "PCR coupling: guard CCA-only accessors"
```

---

### Task 4: Cross-validation method pass-through

**Files:**
- Modify: `R/tune_cca.R` — add `method` to `tune_cca()` and `evaluate_fold()`, thread to `couple()`
- Test: `tests/testthat/test-pcr.R`

- [ ] **Step 1: Write the failing tests**

APPEND to `tests/testthat/test-pcr.R`:

```r
test_that("tune_cca(method = 'pcr') runs and selects sensible k", {
  set.seed(3)
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))
  cv <- prep_cv_folds(coarse, prism, kfolds = 3,
                      max_k_pred = 4, max_k_resp = 4, weight = FALSE)
  res <- tune_cca(cv, k_pred = 2:3, k_resp = 2:3, method = "pcr")
  expect_true(all(c("rmse", "cor_spatial", "cor_temporal") %in% names(res)))
  s <- summarize_cv(res, metric = "rmse")
  expect_true(all(c("k_pred", "k_resp") %in% names(s)))
})

test_that("tune_cca(method = 'pcr') handles a multivariate response", {
  set.seed(4)
  prism_mv <- make_multivar(prism)
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))
  cv <- prep_cv_folds(coarse, prism_mv, kfolds = 3,
                      max_k_pred = 4, max_k_resp = 4,
                      scale_resp = TRUE, weight = FALSE)
  res <- tune_cca(cv, k_pred = 2:3, k_resp = 2:3, method = "pcr")
  expect_true(all(c("rmse", "rmse_tmean", "rmse_ppt") %in% names(res)))
})

test_that("PCR and CCA give comparable skill at full rank (wiring sanity)", {
  set.seed(5)
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))
  pred <- patterns(filter(coarse, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  resp <- patterns(filter(prism, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  test_pred <- filter(coarse, time > as.Date("2018-12-01"))
  test_resp <- filter(prism, time > as.Date("2018-12-01"))

  cca <- couple(pred, resp, method = "cca", k = 3)
  pcr <- couple(pred, resp, method = "pcr")
  m_cca <- tidyeof:::compute_spatial_metrics(predict(cca, test_pred), test_resp, "rmse")$rmse
  m_pcr <- tidyeof:::compute_spatial_metrics(predict(pcr, test_pred), test_resp, "rmse")$rmse
  # All canonical modes retained => CCA equals multivariate OLS, so PCR matches closely
  expect_lt(abs(m_cca - m_pcr) / m_cca, 0.5)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "pcr")'`
Expected: FAIL — `tune_cca()` has no `method` argument, so `method = "pcr"` is an unused argument error.

- [ ] **Step 3: Add `method` to `tune_cca()` and thread it through**

In `R/tune_cca.R`, change the `tune_cca` signature (`:221-226`) to add `method`:

```r
tune_cca <- function(cv_folds,
                     k_pred = 1:10,
                     k_resp = 1:10,
                     k_cca = NULL,
                     method = "cca",
                     metrics = c("rmse", "cor_spatial", "cor_temporal"),
                     parallel = FALSE) {
```

In the `evaluate_fold()` call inside `tune_cca` (`:288-294`), add `method = method`:

```r
      evaluate_fold(
        fold = fold,
        k_pred = params$k_pred,
        k_resp = params$k_resp,
        k_cca = params$k_cca,
        method = method,
        metrics = metrics
      )
```

Change the `evaluate_fold` signature and the `couple()` call inside it (`:320-326`):

```r
evaluate_fold <- function(fold, k_pred, k_resp, k_cca, method = "cca", metrics) {
  # Truncate patterns to requested k (cheap operation using [.patterns)
  pred_patterns <- fold$train_pred_patterns[1:k_pred]
  resp_patterns <- fold$train_resp_patterns[1:k_resp]

  # Couple patterns (k_cca is inert for method = "pcr")
  coupled <- couple(pred_patterns, resp_patterns, k = k_cca,
                    method = method, validate = FALSE)
```

(The rest of `evaluate_fold` — `predict`, `compute_spatial_metrics`, the metric loop — is unchanged.)

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "pcr")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-tune_cca.R` uses the default `method = "cca"` and must stay green)

- [ ] **Step 5: Commit**

```bash
git add R/tune_cca.R tests/testthat/test-pcr.R
git commit -m "PCR coupling: method pass-through in tune_cca"
```

---

### Task 5: Documentation, NEWS, regeneration

**Files:**
- Modify: roxygen in `R/couple.R` (`couple` `@param method`, `predict.coupled_patterns` `@param k`), `R/tune_cca.R` (`tune_cca` `@param method`)
- Modify: `NEWS.md`
- Regenerate: `man/`, `NAMESPACE` via `devtools::document()`

- [ ] **Step 1: Update roxygen for `couple()`**

In `R/couple.R`, replace the `@param method` line of `couple`'s roxygen (currently `#' @param method Coupling method. Currently only "cca" is supported`) with:

```
#' @param method Coupling method: `"cca"` (canonical correlation, default) or
#'   `"pcr"` (principal components regression — OLS of the predictand PC
#'   amplitudes on the predictor PC amplitudes). For PCR the regularization is
#'   the predictor truncation `k_pred`, and the canonical-mode argument `k` is
#'   ignored.
```

In `predict.coupled_patterns`'s roxygen, replace the `@param k` line (currently `#' @param k Number of CCA modes to use for prediction. If NULL, uses all available modes`) with:

```
#' @param k Number of CCA modes to use for prediction (CCA only; if NULL, uses
#'   all available modes). Ignored for `method = "pcr"`.
```

- [ ] **Step 2: Update roxygen for `tune_cca()`**

In `R/tune_cca.R`, add a `@param method` line to `tune_cca`'s roxygen (place it just after the `@param k_cca` block):

```
#' @param method Coupling method passed to [couple()]: `"cca"` (default) or
#'   `"pcr"`. For `"pcr"` the `k_cca` axis is inert — leave `k_cca = NULL` so
#'   the grid is effectively `k_pred` x `k_resp` (an explicit `k_cca` vector
#'   would produce duplicate rows that all evaluate identically).
```

- [ ] **Step 3: Add a NEWS entry**

In `NEWS.md`, under the existing `# tidyeof (development version)` heading (add the heading if it is not present), add a bullet:

```markdown
* `couple()` and `tune_cca()` gain `method = "pcr"` for principal components
  regression — an OLS alternative to CCA that maps predictor PC amplitudes to
  predictand PC amplitudes. Regularization is the predictor truncation
  (`k_pred`); the CCA-style coupling `k` is inert for PCR. CCA-specific
  diagnostics (`get_canonical_*`) are not defined for PCR couplings.
```

- [ ] **Step 4: Regenerate docs and run the suite**

Run: `Rscript -e 'devtools::document()'` — Expected: `man/couple.Rd`, `man/predict.coupled_patterns.Rd`, `man/tune_cca.Rd` regenerate; NAMESPACE unchanged (no new exports — `fit_pcr`, `apply_pcr_prediction`, `check_cca_method` are `@keywords internal`).
Run: `Rscript -e 'devtools::test()'` — Expected: full suite green.

- [ ] **Step 5: Confirm NAMESPACE is unchanged and commit**

```bash
git status --short        # confirm NAMESPACE not modified; vignettes/eof-analysis.qmd NOT staged
git add R/couple.R R/tune_cca.R man NEWS.md
git commit -m "document PCR coupling method"
```

---

## Self-review notes

- **Spec coverage:** estimator = OLS (Task 1 `fit_pcr`); `k` inert (Tasks 1, 2, 4); rank-safe QR + `tidyeof_rank_deficient` (Task 1); `pcr` slot on reused class (Task 1); `predict` dispatch (Task 2); cross-source override (Task 2); multivariate predictand (Tasks 2, 4); CCA-only accessor guards + `tidyeof_cca_only` (Task 3); method-aware print/summary (Task 1); `tune_cca(method=)` with inert `k_cca` (Task 4); docs/NEWS (Task 5). All spec sections map to a task.
- **No new dependency, no NAMESPACE change** (internal helpers).
- **Univariate/CCA preservation:** the CCA branches of `couple`/`predict` and all accessors are byte-identical for `method = "cca"`; only the dispatch wrapping is added.
