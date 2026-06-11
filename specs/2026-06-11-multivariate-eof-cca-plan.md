# Multivariate EOF/CCA Support Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `patterns()` and the full CCA downscaling stack accept multi-attribute stars objects (e.g., temperature + precipitation), producing joint EOFs with per-variable spatial loadings and shared amplitudes.

**Architecture:** Attribute-aware core (Approach A from `specs/2026-06-11-multivariate-eof-cca-design.md`). Each variable is flattened to a `time × p` block and cbind'd variable-major into one `time × (V·p)` matrix; PCA/rotation/projection are unchanged. A `block_map` on the `patterns` object records which columns belong to which variable. Univariate is the single-block special case — one code path everywhere.

**Tech Stack:** R package; stars/sf/units, tidyverse, testthat. Test data: `inst/testdata/prism_test.RDS` — single attribute `tmean` (units °C), dims x=51, y=51, time=36 (monthly Dates 2017-01-01 to 2019-12-01), 2601 cells, no NA cells.

**Verified facts (do not re-derive):**
- `Ops.stars` FAILS for two multi-attribute operands (`dat - clim$mean` errors with `length(e2) == 1 is not TRUE`) → all stars-arithmetic must loop attributes.
- These DO work on multi-attribute stars: `aperm()` (permutes every attribute), `dplyr::filter()`, `dplyr::slice()`, `units::drop_units()`, `st_apply(..., rename = FALSE)` (returns all attributes), `stars * numeric_vector` (applies per attribute), `do.call(c, <list of single-attribute stars>)` (recombines attributes).
- No test currently references the `"weight"` EOF attribute name; renaming it (approved) only touches `R/patterns_core.R`.
- `tests/testthat.R` defines `prism` for `R CMD check`; each test file ALSO loads `prism` itself at the top (house convention — follow it). `helper-*.R` files are auto-sourced by testthat.

**Run tests with:** `Rscript -e 'devtools::test(filter = "<name>")'` from the package root. Quote the path when using cd: `cd "/Users/nick/UF Dropbox/Nick Gauthier/Projects/tidyeof"`.

---

### Task 1: Multivariate flattening helpers

**Files:**
- Modify: `R/stars_utils.R` (`flatten_time_space`, lines 11-37)
- Modify: `R/get_climatology.R` (`flatten_dim_space`, lines 24-28)
- Create: `tests/testthat/helper-multivariate.R`
- Create: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Create the test helper**

`tests/testthat/helper-multivariate.R`:

```r
# Build a two-variable stars object from the single-variable PRISM test data.
# The second attribute is a deterministic nonlinear transform of temperature
# (no RNG, so helpers stay reproducible without touching the seed).
make_multivar <- function(x) {
  vals <- units::drop_units(x[[1]])
  idx <- array(seq_along(vals), dim = dim(vals))
  ppt <- exp(0.15 * vals) + 3 * sin(idx / 7)
  out <- x
  out$ppt <- units::set_units(ppt, "mm")
  out
}
```

- [ ] **Step 2: Write the failing tests**

`tests/testthat/test-multivariate.R` (start the file with the house-convention prism load):

```r
prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()
prism_mv <- make_multivar(prism)

test_that("flatten_time_space concatenates attributes variable-major with a block map", {
  flat <- tidyeof:::flatten_time_space(units::drop_units(prism_mv))
  expect_equal(dim(flat$matrix), c(36, 2 * 2601))
  expect_named(flat$block_map, c("tmean", "ppt"))
  expect_equal(flat$block_map$tmean, 1:2601)
  expect_equal(flat$block_map$ppt, 2602:5202)
  expect_equal(flat$n_space, 2601)

  # Each block equals flattening that variable alone
  flat_t <- tidyeof:::flatten_time_space(units::drop_units(prism_mv["tmean"]))
  flat_p <- tidyeof:::flatten_time_space(units::drop_units(prism_mv["ppt"]))
  expect_equal(flat$matrix[, flat$block_map$tmean], flat_t$matrix)
  expect_equal(flat$matrix[, flat$block_map$ppt], flat_p$matrix)
})

test_that("flatten_time_space is unchanged for univariate input", {
  flat <- tidyeof:::flatten_time_space(units::drop_units(prism))
  expect_equal(dim(flat$matrix), c(36, 2601))
  expect_equal(flat$block_map, list(tmean = 1:2601))
  expect_equal(flat$spatial_dims, c("x", "y"))
  expect_equal(unname(flat$spatial_shape), c(51L, 51L))
})

test_that("flatten_dim_space concatenates attributes", {
  m <- tidyeof:::flatten_dim_space(prism_mv, "time")
  expect_equal(dim(m), c(36, 2 * 2601))
  m1 <- tidyeof:::flatten_dim_space(prism_mv["tmean"], "time")
  expect_equal(m[, 1:2601], m1)
})
```

- [ ] **Step 3: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — `flat$block_map` is NULL (old function has no block_map), and `flatten_dim_space` returns 2601 columns (it reads only `[[1]]`).

- [ ] **Step 4: Implement `flatten_time_space`**

Replace the whole function body in `R/stars_utils.R` (keep the roxygen block, update its text to mention multiple attributes and the returned `block_map`/`n_space`):

```r
flatten_time_space <- function(dat) {
  check_stars_object(dat)

  dims <- stars::st_dimensions(dat)
  if (!"time" %in% names(dims)) {
    cli::cli_abort("Object must contain a {.field time} dimension to be flattened.")
  }

  spatial_dims <- setdiff(names(dims), "time")
  permuted <- aperm(dat, c("time", spatial_dims))

  var_names <- names(dat)
  blocks <- purrr::map(seq_along(var_names), function(i) {
    arr <- permuted[[i]]
    matrix(arr, nrow = dim(arr)[1], ncol = prod(dim(arr)[-1]))
  })
  mat <- do.call(cbind, blocks)

  n_space <- ncol(blocks[[1]])
  block_map <- setNames(
    purrr::map(seq_along(var_names), ~((.x - 1L) * n_space + 1L):(.x * n_space)),
    var_names
  )

  spatial_values <- purrr::map(spatial_dims, ~stars::st_get_dimension_values(permuted, .x))
  names(spatial_values) <- spatial_dims

  list(
    matrix = mat,
    block_map = block_map,
    n_space = n_space,
    spatial_dims = spatial_dims,
    spatial_shape = dim(permuted[[1]])[-1],
    spatial_values = spatial_values
  )
}
```

Note the old `if (length(dat) != 1) dat <- dat[1]` is gone. Existing callers that pass `dat[1]` keep working (single block).

- [ ] **Step 5: Implement `flatten_dim_space`**

Replace the function body in `R/get_climatology.R` (update roxygen to say attributes are concatenated variable-major):

```r
flatten_dim_space <- function(x, dim_name) {
  spatial <- setdiff(names(stars::st_dimensions(x)), dim_name)
  permuted <- aperm(units::drop_units(x), c(dim_name, spatial))
  blocks <- purrr::map(seq_along(names(x)), function(i) {
    arr <- permuted[[i]]
    matrix(arr, nrow = dim(arr)[1])
  })
  do.call(cbind, blocks)
}
```

- [ ] **Step 6: Run the new tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (univariate behavior unchanged)

- [ ] **Step 7: Commit**

```bash
git add R/stars_utils.R R/get_climatology.R tests/testthat/helper-multivariate.R tests/testthat/test-multivariate.R
git commit -m "multivariate flattening helpers with block map"
```

---

### Task 2: Multivariate climatology and anomalies

**Files:**
- Modify: `R/get_climatology.R` (`get_climatology` monthly branch, `apply_monthly_climatology`, `get_anomalies` annual branch, `restore_climatology` annual branch)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("multivariate anomalies equal per-variable anomalies (annual and monthly)", {
  for (m in c(FALSE, TRUE)) {
    anom_mv <- get_anomalies(prism_mv, scale = TRUE, monthly = m)
    anom_t <- get_anomalies(prism_mv["tmean"], scale = TRUE, monthly = m)
    anom_p <- get_anomalies(prism_mv["ppt"], scale = TRUE, monthly = m)
    expect_named(anom_mv, c("tmean", "ppt"))
    expect_equal(units::drop_units(anom_mv)[["tmean"]],
                 units::drop_units(anom_t)[[1]], tolerance = 1e-12)
    expect_equal(units::drop_units(anom_mv)[["ppt"]],
                 units::drop_units(anom_p)[[1]], tolerance = 1e-12)
  }
})

test_that("multivariate climatology round-trips through restore_climatology", {
  for (m in c(FALSE, TRUE)) {
    clim <- get_climatology(prism_mv, monthly = m)
    expect_named(clim$mean, c("tmean", "ppt"))
    anom <- get_anomalies(prism_mv, clim, scale = TRUE, monthly = m)
    restored <- restore_climatology(anom, clim, scale = TRUE, monthly = m)
    expect_equal(units::drop_units(restored[["tmean"]]),
                 units::drop_units(prism_mv[["tmean"]]), tolerance = 1e-8)
    expect_equal(units::drop_units(restored[["ppt"]]),
                 units::drop_units(prism_mv[["ppt"]]), tolerance = 1e-8)
  }
})

test_that("get_anomalies rejects climatology with mismatched attribute count", {
  clim <- get_climatology(prism_mv["tmean"])
  expect_error(get_anomalies(prism_mv, clim), class = "tidyeof_attribute_mismatch")
})
```

Note on `units::drop_units(anom_mv)`: annual scaled anomalies carry dimensionless `units` objects while monthly anomalies are bare numeric — dropping units before comparison makes the test path-independent. `units::drop_units` on a unitless attribute is a no-op for these comparisons because we apply it to the whole object on both sides.

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — annual path errors with `length(e2) == 1 is not TRUE` (Ops.stars), monthly path silently uses only `tmean`.

- [ ] **Step 3: Fix `get_anomalies` annual branch**

In `R/get_climatology.R`, replace the end of `get_anomalies` (currently `out <- dat - clim$mean; if (scale) out <- out / clim$sd; out`) with a positional per-attribute loop, and add the count guard before the `if (monthly)` branch (so both paths get it):

```r
  if (length(clim$mean) != length(dat)) {
    cli::cli_abort(
      "Climatology has {length(clim$mean)} attribute{?s} but the data has {length(dat)}.",
      class = "tidyeof_attribute_mismatch"
    )
  }
```

and at the end of the function:

```r
  # Ops.stars cannot subtract multi-attribute objects, so anomalize each
  # attribute separately (paired by position, preserving the univariate
  # behavior of ignoring attribute names) and recombine.
  out <- do.call(c, purrr::map(seq_along(names(dat)), function(i) {
    out_i <- dat[i] - clim$mean[i]
    if (scale) out_i <- out_i / clim$sd[i]
    out_i
  }))
  out
```

- [ ] **Step 4: Fix `restore_climatology` annual branch**

Replace the block after the units-mismatch handling (currently `if (scale) anomalies <- anomalies * target_sd; out <- anomalies + target_mean`) with:

```r
  out <- do.call(c, purrr::map(seq_along(names(anomalies)), function(i) {
    out_i <- anomalies[i]
    if (scale) out_i <- out_i * target_sd[i]
    out_i + target_mean[i]
  }))
```

(Keep the existing `restore_units(out, clim$mean)` call after it.)

- [ ] **Step 5: Fix the monthly climatology**

In `get_climatology`'s monthly branch, change `mat <- flatten_dim_space(dat[1], "time")` to `mat <- flatten_dim_space(dat, "time")` (now `n_space` in that branch is V·p), and replace the `to_stars` helper so it splits blocks back into per-variable attributes:

```r
    new_dims <- st_dimensions(dat)[spatial]
    new_dims$month <- month_dimension()
    class(new_dims) <- "dimensions"
    spatial_shape <- dim(dat)[spatial]
    n_cells <- prod(spatial_shape)

    to_stars <- function(values) {
      per_var <- purrr::map(seq_along(names(dat)), function(i) {
        block <- values[((i - 1L) * n_cells + 1L):(i * n_cells), , drop = FALSE]
        stars::st_as_stars(
          array(block, dim = c(spatial_shape, month = 12L)),
          dimensions = new_dims
        ) %>%
          setNames(names(dat)[i])
      })
      do.call(c, per_var)
    }
```

(`month_stat`'s `n_space` local stays as `ncol(mat)` — the full V·p width — no change needed there.)

- [ ] **Step 6: Fix `apply_monthly_climatology`**

Replace the flatten/rebuild section (the lines from `spatial <- setdiff(...)` through `aperm(permuted, ...)`) with a per-attribute version. The arithmetic in the middle is unchanged:

```r
  spatial <- setdiff(names(st_dimensions(dat)), "time")
  permuted <- aperm(units::drop_units(dat), c("time", spatial))
  blocks <- purrr::map(seq_along(names(dat)), function(i) {
    arr <- permuted[[i]]
    matrix(arr, nrow = length(times))
  })
  mat <- do.call(cbind, blocks)
  n_cells <- ncol(blocks[[1]])

  if (ncol(mat) != ncol(mn_mat)) {
    cli::cli_abort(
      "Spatial size mismatch: data has {ncol(mat)} cells but the climatology has {ncol(mn_mat)}.",
      class = "tidyeof_grid_mismatch"
    )
  }

  # ... (existing missing-month check and anomalize/restore arithmetic on `mat`,
  #      unchanged) ...

  for (i in seq_along(names(dat))) {
    cols <- ((i - 1L) * n_cells + 1L):(i * n_cells)
    permuted[[i]] <- array(mat[, cols, drop = FALSE], dim = dim(permuted[[i]]))
  }
  aperm(permuted, names(st_dimensions(dat)))
```

(`mn_mat`/`sd_mat` already come from the now-multivariate `flatten_dim_space`, so their width is V·p too. The `ncol` equality check therefore also catches attribute-count mismatches in the monthly path.)

- [ ] **Step 7: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (watch `test-get_climatology.R` and `test-monthly-alignment.R` especially)

- [ ] **Step 8: Commit**

```bash
git add R/get_climatology.R tests/testthat/test-multivariate.R
git commit -m "multivariate climatology, anomalies, and restoration"
```

---

### Task 3: Sign alignment over concatenated blocks

**Files:**
- Modify: `R/align_patterns.R` (`compute_eof_signs`, `apply_sign_flips`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing test**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("compute_eof_signs works on multi-attribute EOF objects", {
  # Hand-build a 2-attribute (x, y, PC) stars object with known column sums
  arr_pos <- array(1, dim = c(3, 3, 2))      # both PCs sum positive
  arr_neg <- array(-1, dim = c(3, 3, 2))     # both PCs sum negative
  eofs <- stars::st_as_stars(list(a = arr_pos, b = arr_neg)) %>%
    stars::st_set_dimensions(names = c("x", "y", "PC"))
  # Block sums cancel: 9 + (-9) = 0 per PC -> orient() maps 0 to +1
  expect_equal(unname(tidyeof:::compute_eof_signs(eofs)), c(1, 1))

  # Make b dominate negatively
  eofs$b <- arr_neg * 3
  expect_equal(unname(tidyeof:::compute_eof_signs(eofs)), c(-1, -1))
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — old `compute_eof_signs` reads `split(eofs)` / `eofs[[1]]`, ignoring attribute `b` (returns `c(1, 1)` for the second case).

- [ ] **Step 3: Implement**

Replace both functions in `R/align_patterns.R` (keep roxygen, note multi-attribute support):

```r
compute_eof_signs <- function(eofs) {
  # Orient so the dominant loading across ALL variable blocks is positive.
  # A zero-sum pattern maps to +1 rather than sign(0) = 0, which would zero
  # out the mode.
  orient <- function(x) if (sum(x, na.rm = TRUE) >= 0) 1 else -1
  mat <- do.call(rbind, purrr::map(names(eofs), function(v) {
    arr <- eofs[[v]]
    matrix(arr, nrow = prod(dim(arr)[-length(dim(arr))]))
  }))
  signs <- apply(mat, 2, orient)
  names(signs) <- paste0("PC", seq_along(signs))
  signs
}
```

(This single matrix-based path replaces the old raster/geometry branches: for raster the attribute array is `(x, y, PC)` so `matrix()` collapses x·y rows; for geometry it is `(geometry, PC)` and is already a matrix.)

In `apply_sign_flips`, replace the single `sweep` line on `patterns$eofs` with a per-attribute loop:

```r
  for (v in names(patterns$eofs)) {
    arr <- patterns$eofs[[v]]
    patterns$eofs[[v]] <- sweep(arr, length(dim(arr)), signs, `*`)
  }
```

(The amplitude and proj_matrix flips below it are unchanged.)

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add R/align_patterns.R tests/testthat/test-multivariate.R
git commit -m "sign alignment over concatenated variable blocks"
```

---

### Task 4: Multivariate `patterns()` core

**Files:**
- Modify: `R/patterns_core.R` (new `check_multivariate_scale`, `patterns`, `get_eofs`)
- Modify: `R/patterns_class.R` (`new_patterns` gains `block_map`; `print.patterns` lists variables)
- Modify: `R/common_patterns.R` (pass `block_map` through `new_patterns`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("patterns() requires scale = TRUE for multivariate input", {
  expect_error(patterns(prism_mv, k = 3), class = "tidyeof_multivariate_scale")
  expect_error(patterns(prism_mv, k = 3, scale = FALSE), class = "tidyeof_multivariate_scale")
})

test_that("multivariate patterns have per-variable EOFs and shared amplitudes", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  expect_s3_class(pat, "patterns")
  expect_named(pat$eofs, c("tmean", "ppt"))
  expect_equal(pat$names, c("tmean", "ppt"))
  expect_equal(pat$block_map, list(tmean = 1:2601, ppt = 2602:5202))
  expect_equal(unname(dim(pat$eofs[["tmean"]])), c(51, 51, 3))
  expect_named(pat$amplitudes, c("time", "PC1", "PC2", "PC3"))
  expect_equal(nrow(pat$proj_matrix), length(pat$valid_pixels))
  expect_equal(pat$units$ppt, units(prism_mv[["ppt"]]))
})

test_that("duplicated variable yields identical loading blocks", {
  dup <- prism
  dup$tmean2 <- prism[[1]]
  pat <- patterns(dup, k = 3, scale = TRUE, weight = FALSE)
  expect_equal(pat$eofs[["tmean"]], pat$eofs[["tmean2"]], tolerance = 1e-8)
})

test_that("univariate EOF attribute is named after the variable", {
  pat <- patterns(prism, k = 2)
  expect_named(pat$eofs, "tmean")
  expect_equal(pat$block_map, list(tmean = 1:2601))
})

test_that("non-finite cells from sd = 0 standardization are dropped", {
  z <- units::drop_units(prism)
  arr <- z[[1]]
  arr[1, 1, ] <- 5  # constant cell -> sd 0 -> NaN anomalies under scaling
  z[[1]] <- arr
  pat <- patterns(z, k = 2, scale = TRUE, weight = FALSE)
  expect_false(1 %in% pat$valid_pixels)
  expect_true(all(is.finite(pat$proj_matrix)))
})

test_that("monthly multivariate patterns run end to end", {
  pat <- patterns(prism_mv, k = 2, scale = TRUE, monthly = TRUE)
  expect_named(pat$eofs, c("tmean", "ppt"))
  expect_equal(nrow(pat$amplitudes), 36)
})

test_that("per-variable NA masks are handled independently", {
  z <- prism_mv
  arr <- units::drop_units(z[["ppt"]])
  arr[2, 1, ] <- NA
  z$ppt <- units::set_units(arr, "mm")
  pat <- patterns(z, k = 2, scale = TRUE, weight = FALSE)
  expect_false((2601 + 2) %in% pat$valid_pixels)  # ppt block, cell 2
  expect_true(2 %in% pat$valid_pixels)            # tmean cell 2 still valid
})

test_that("multivariate rotation runs and reorders blocks together", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE, rotate = TRUE)
  expect_named(pat$eofs, c("tmean", "ppt"))
  expect_equal(ncol(pat$rotation), 3)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — `patterns(prism_mv, ...)` aborts with `tidyeof_multiple_attributes` (and no `tidyeof_multivariate_scale` class exists yet).

- [ ] **Step 3: Add the scale gate and rewire `patterns()`**

In `R/patterns_core.R`, add below `check_single_attribute` (which STAYS — `common_patterns()` still uses it):

```r
#' Check that multivariate input is standardized
#'
#' PCA is variance-driven, so joint EOFs of variables with different units
#' require per-pixel standardization to contribute comparably.
#' @param dat A stars object
#' @param scale The scale argument passed to the caller
#' @param call Calling environment for error messages
#' @keywords internal
check_multivariate_scale <- function(dat, scale, call = rlang::caller_env()) {
  if (length(dat) > 1 && !isTRUE(scale)) {
    cli::cli_abort(
      c(
        "Multivariate input ({length(dat)} attributes: {.field {names(dat)}}) requires {.code scale = TRUE}.",
        "i" = "PCA is variance-driven; variables with different units must be standardized to contribute comparably."
      ),
      class = "tidyeof_multivariate_scale",
      call = call
    )
  }
}
```

In `patterns()`, replace `check_single_attribute(dat)` with `check_multivariate_scale(dat, scale)`, and pass the block map to the constructor: add `block_map = eofs$block_map,` to the `new_patterns(...)` call.

- [ ] **Step 4: Rewrite the multivariate parts of `get_eofs`**

Replace the top of `get_eofs` (from `flattened <- ...` through the weights block) with:

```r
  var_names <- names(dat)
  n_vars <- length(var_names)

  flattened <- flatten_time_space(units::drop_units(dat))
  anomaly_matrix_full <- flattened$matrix
  n_pixels <- ncol(anomaly_matrix_full)   # n_vars * n_space
  n_space <- flattened$n_space
  block_map <- flattened$block_map

  # Valid pixels must be finite everywhere: this drops NA-masked cells and
  # also Inf/NaN cells produced by sd ~ 0 standardization (e.g. arid cells
  # for precipitation)
  valid_pixels <- which(apply(anomaly_matrix_full, 2,
                              function(col) all(is.finite(col))))

  # Validate k value
  max_k <- min(length(times) - 1, length(valid_pixels))
  check_k_valid(k, max_k)

  # Extract matrix for valid pixels (time x space)
  anomaly_matrix <- anomaly_matrix_full[, valid_pixels, drop = FALSE]

  # Apply spatial weights column-wise if provided; one weight per grid cell,
  # replicated across variable blocks
  if (!is.null(weights)) {
    if (length(weights) != n_space) {
      cli::cli_abort(
        "Length of {.arg weights} ({length(weights)}) must match number of spatial points per variable ({n_space}).",
        class = "tidyeof_weight_mismatch"
      )
    }
    weights_valid <- rep(weights, n_vars)[valid_pixels]
    anomaly_matrix <- sweep(anomaly_matrix, 2, weights_valid, `*`)
  } else {
    weights_valid <- rep(1, length(valid_pixels))
  }
```

Then replace the EOF-stars construction section (from `full_patterns <- ...` through the end of the raster/geometry `if/else`) with:

```r
  full_patterns <- array(NA, dim = c(k, n_pixels))
  full_patterns[, valid_pixels] <- t(loadings)

  # Build a multi-attribute template (one attribute per variable) with the
  # time dimension relabeled as PC, then fill each attribute with its block.
  # Attribute names come from the data, so univariate EOFs are named after
  # their variable (previously "weight").
  if (has_geometry_dimension(dat)) {
    template <- dat[, , 1:k, drop = FALSE] %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC')
    for (i in seq_along(var_names)) {
      template[[i]] <- t(full_patterns[, block_map[[i]], drop = FALSE])  # geometry x PC
    }
  } else {
    template <- dat[, , , 1:k, drop = FALSE] %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC')
    for (i in seq_along(var_names)) {
      pattern_array <- array(full_patterns[, block_map[[i]], drop = FALSE],
                             dim = c(k, dims[[1]], dims[[2]]))
      template[[i]] <- aperm(pattern_array, c(2, 3, 1))  # x, y, PC
    }
  }
  spatial_patterns <- template
```

Finally add `block_map = block_map,` to the returned list at the bottom of `get_eofs`.

- [ ] **Step 5: Thread `block_map` through the constructor**

In `R/patterns_class.R`, add the parameter `block_map = NULL` to `new_patterns()` (document with `@param block_map Named list mapping each variable to its column range in the concatenated space-time matrix`), and store `block_map = block_map,` in the structure list (next to `valid_pixels`).

In `print.patterns`, after the `Modes:` line add:

```r
  cli::cli_text("Variables: {.field {x$names}}")
```

In `R/common_patterns.R`, add `block_map = eofs$block_map,` to its `new_patterns(...)` call (sources are univariate, so this is a single-entry map).

- [ ] **Step 6: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass. If anything fails it will be in plotting or reconstruction reading `eofs[[1]]` — those still work because `[[1]]` is positional, but investigate any failure before proceeding.

- [ ] **Step 7: Commit**

```bash
git add R/patterns_core.R R/patterns_class.R R/common_patterns.R tests/testthat/test-multivariate.R
git commit -m "multivariate EOF extraction in patterns()"
```

---

### Task 5: `eof_loading_matrix` + multivariate `reconstruct()`

**Files:**
- Modify: `R/stars_utils.R` (add generalized `eof_loading_matrix`)
- Modify: `R/tune_cca.R` (delete the old `eof_loading_matrix`, lines ~560-569)
- Modify: `R/reconstruct_field.R` (`reconstruct`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("multivariate reconstruction round-trips at full rank", {
  pat <- patterns(prism_mv, k = 35, scale = TRUE, weight = FALSE)
  rec <- reconstruct(pat)
  expect_named(rec, c("tmean", "ppt"))
  expect_equal(units(rec[["ppt"]]), units(prism_mv[["ppt"]]))
  expect_equal(units::drop_units(rec[["tmean"]]),
               units::drop_units(prism_mv[["tmean"]]), tolerance = 1e-6)
  expect_equal(units::drop_units(rec[["ppt"]]),
               units::drop_units(prism_mv[["ppt"]]), tolerance = 1e-6)
})

test_that("truncated multivariate reconstruction returns both variables with weighting", {
  pat <- patterns(prism_mv, k = 4, scale = TRUE)
  rec <- reconstruct(pat)
  expect_named(rec, c("tmean", "ppt"))
  expect_equal(unname(dim(rec)), c(51, 51, 36))
})

test_that("eof_loading_matrix stacks variable blocks", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  m <- tidyeof:::eof_loading_matrix(pat)
  expect_equal(dim(m), c(2 * 2601, 3))
  expect_equal(m[1:2601, ], matrix(pat$eofs[["tmean"]], nrow = 2601, ncol = 3))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — `reconstruct` reads `eofs[[1]]` only, so its matrix has 2601 rows vs 5202-row `valid_pixels` indexing (subscript out of bounds or wrong values), and the round-trip returns only one variable's worth of data.

- [ ] **Step 3: Move and generalize `eof_loading_matrix`**

Delete the existing `eof_loading_matrix` from `R/tune_cca.R`. Add to `R/stars_utils.R`:

```r
#' Extract the EOF loading matrix (space x PC) from a patterns object
#'
#' Stacks every variable's loadings into the concatenated variable-major
#' layout used by [flatten_time_space()], so rows align with
#' `patterns$valid_pixels` and `patterns$block_map`.
#' @param patterns A patterns object
#' @return A numeric matrix with prod(spatial) * n_vars rows and k columns
#' @keywords internal
eof_loading_matrix <- function(patterns) {
  do.call(rbind, purrr::map(names(patterns$eofs), function(v) {
    arr <- patterns$eofs[[v]]
    matrix(arr, nrow = prod(dim(arr)[-length(dim(arr))]), ncol = patterns$k)
  }))
}
```

- [ ] **Step 4: Rewrite `reconstruct()`'s spatial rebuild**

In `R/reconstruct_field.R`, replace the section from `eof_array <- target_patterns$eofs[[1]]` through the `matrix_to_spacetime(...)` call with:

```r
  eof_matrix <- eof_loading_matrix(target_patterns)
  valid_pixels <- target_patterns$valid_pixels
  eof_valid <- eof_matrix[valid_pixels, , drop = FALSE]

  anomalies_valid <- amps_matrix %*% t(eof_valid)

  # Re-insert into the full concatenated space, then split per variable
  n_total <- nrow(eof_matrix)
  full_mat <- matrix(NA_real_, nrow = nrow(anomalies_valid), ncol = n_total)
  full_mat[, valid_pixels] <- anomalies_valid

  block_map <- target_patterns$block_map
  if (is.null(block_map)) {
    # Patterns objects from before block_map existed are univariate
    block_map <- setNames(list(seq_len(n_total)), target_patterns$names[[1]])
  }

  var_list <- purrr::map(seq_along(block_map), function(i) {
    matrix_to_spacetime(
      full_mat[, block_map[[i]], drop = FALSE],
      template_eofs = target_patterns$eofs[i],
      spatial_template = target_patterns$climatology$mean[i],
      valid_pixels = seq_along(block_map[[i]]),
      times = amplitudes$time,
      var_names = names(block_map)[[i]]
    )
  })
  anomalies <- do.call(c, var_list)
```

(Everything below — `restore_climatology` and the per-variable units loop — is unchanged; both already handle multiple attributes after Task 2.)

- [ ] **Step 5: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-get_patterns.R` and `test-tune_eof.R` exercise the moved helper)

- [ ] **Step 6: Commit**

```bash
git add R/stars_utils.R R/tune_cca.R R/reconstruct_field.R tests/testthat/test-multivariate.R
git commit -m "multivariate reconstruction via shared loading-matrix helper"
```

---

### Task 6: Multivariate `project_patterns()`

**Files:**
- Modify: `R/project_patterns.R` (lines 22-62)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("projecting training data reproduces stored amplitudes (multivariate)", {
  pat <- patterns(prism_mv, k = 4, scale = TRUE)
  proj <- project_patterns(pat, prism_mv)
  expect_equal(as.matrix(proj[-1]), as.matrix(pat$amplitudes[-1]),
               tolerance = 1e-6, ignore_attr = TRUE)
})

test_that("project_patterns reorders newdata attributes to training order", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  expect_equal(project_patterns(pat, prism_mv[c("ppt", "tmean")]),
               project_patterns(pat, prism_mv))
})

test_that("project_patterns rejects mismatched attribute sets", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  bad <- setNames(prism_mv, c("tmean", "precip"))
  expect_error(project_patterns(pat, bad), class = "tidyeof_attribute_mismatch")
})

test_that("univariate projection stays name-agnostic", {
  pat <- patterns(prism, k = 3)
  renamed <- setNames(prism, "tas")
  expect_equal(project_patterns(pat, renamed), project_patterns(pat, prism))
})

test_that("rotated multivariate projection reproduces stored amplitudes", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE, rotate = TRUE)
  proj <- project_patterns(pat, prism_mv)
  expect_equal(as.matrix(proj[-1]), as.matrix(pat$amplitudes[-1]),
               tolerance = 1e-6, ignore_attr = TRUE)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — `project_patterns` aborts with `tidyeof_multiple_attributes` on multivariate newdata.

- [ ] **Step 3: Implement**

In `R/project_patterns.R`, replace `check_single_attribute(newdata)` with:

```r
  # Univariate patterns accept any single-attribute newdata (name-agnostic,
  # e.g. cross-source prediction where the variable is named differently).
  # Multivariate projection needs the same variable set, reordered to match.
  if (length(patterns$names) > 1 || length(newdata) > 1) {
    if (!setequal(names(newdata), patterns$names)) {
      cli::cli_abort(
        c(
          "Attributes of {.arg newdata} ({.field {names(newdata)}}) must match the training variables ({.field {patterns$names}}).",
          "i" = "Multivariate patterns require the same set of variables."
        ),
        class = "tidyeof_attribute_mismatch"
      )
    }
    newdata <- newdata[patterns$names]
  }
```

and change `flattened <- flatten_time_space(anomalies[1])` to `flattened <- flatten_time_space(anomalies)`.

(No other changes: `anomalies * area_weights(newdata)` multiplies each attribute by the per-cell weights — verified stars behavior — and `valid_pixels`/`proj_matrix` already work on concatenated flat indices.)

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-projection_consistency.R` covers the univariate path)

- [ ] **Step 5: Commit**

```bash
git add R/project_patterns.R tests/testthat/test-multivariate.R
git commit -m "multivariate projection with attribute matching"
```

---

### Task 7: Per-variable + pooled metrics

**Files:**
- Modify: `R/metrics.R` (`compute_spatial_metrics`; add `compute_block_metrics`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("multivariate metrics report pooled plus per-variable scores", {
  pat <- patterns(prism_mv, k = 4, scale = TRUE)
  rec <- reconstruct(pat)
  m <- tidyeof:::compute_spatial_metrics(rec, prism_mv)

  expect_true(all(c("rmse", "rmse_tmean", "rmse_ppt",
                    "cor_spatial", "cor_spatial_tmean", "cor_spatial_ppt",
                    "cor_temporal", "cor_temporal_tmean", "cor_temporal_ppt")
                  %in% names(m)))

  # Pooled rmse = RMS of per-variable rmse normalized by observed sd
  sd_t <- sd(units::drop_units(prism_mv[["tmean"]]), na.rm = TRUE)
  sd_p <- sd(units::drop_units(prism_mv[["ppt"]]), na.rm = TRUE)
  expect_equal(m$rmse,
               sqrt(mean(c((m$rmse_tmean / sd_t)^2, (m$rmse_ppt / sd_p)^2))))

  # Pooled correlations are means of per-variable correlations
  expect_equal(m$cor_spatial, mean(c(m$cor_spatial_tmean, m$cor_spatial_ppt)))
  expect_equal(m$cor_temporal, mean(c(m$cor_temporal_tmean, m$cor_temporal_ppt)))
})

test_that("univariate metrics are unchanged", {
  pat <- patterns(prism, k = 4)
  rec <- reconstruct(pat)
  m <- tidyeof:::compute_spatial_metrics(rec, prism)
  expect_named(m, c("rmse", "cor_spatial", "cor_temporal"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — no per-variable columns; pooled rmse is computed over raw mixed-unit cells.

- [ ] **Step 3: Implement `compute_block_metrics`**

Add to `R/metrics.R`:

```r
#' Compute per-variable and pooled metrics over block-structured matrices
#'
#' The plain metric name is always the pooled score. For a single block it is
#' the metric itself, preserving univariate behavior exactly. For multiple
#' blocks, RMSE is normalized by each block's observed standard deviation and
#' RMS-combined (raw pooling across different units is meaningless), and
#' correlations (already unitless) are averaged. Per-variable values get
#' suffixed names (e.g. rmse_tmean) only when there is more than one block.
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#' @param block_map Named list of column indices per variable
#' @param metrics Character vector of metric names
#' @return Named list of metric values
#' @keywords internal
compute_block_metrics <- function(pred_matrix, obs_matrix, block_map,
                                  metrics = c("rmse", "cor_spatial", "cor_temporal")) {
  calc <- list(rmse = calc_rmse, cor_spatial = calc_cor_spatial,
               cor_temporal = calc_cor_temporal)
  metrics <- intersect(metrics, names(calc))

  per_var <- purrr::map(block_map, function(cols) {
    p <- pred_matrix[, cols, drop = FALSE]
    o <- obs_matrix[, cols, drop = FALSE]
    vals <- purrr::map(calc[metrics], ~.x(p, o))
    vals$.obs_sd <- stats::sd(o, na.rm = TRUE)
    vals
  })

  results <- list()
  for (m in metrics) {
    vals <- purrr::map_dbl(per_var, m)
    if (length(per_var) == 1) {
      results[[m]] <- vals[[1]]
    } else {
      if (m == "rmse") {
        nrmse <- vals / purrr::map_dbl(per_var, ".obs_sd")
        results[[m]] <- sqrt(mean(nrmse^2))
      } else {
        results[[m]] <- mean(vals)
      }
      for (v in names(per_var)) {
        results[[paste0(m, "_", v)]] <- per_var[[v]][[m]]
      }
    }
  }
  results
}
```

- [ ] **Step 4: Rewire `compute_spatial_metrics`**

In `compute_spatial_metrics`, add attribute alignment at the top (before time alignment):

```r
  if (length(predicted) > 1 || length(observed) > 1) {
    if (!setequal(names(predicted), names(observed))) {
      cli::cli_abort(
        "Predicted attributes ({.field {names(predicted)}}) must match observed attributes ({.field {names(observed)}}).",
        class = "tidyeof_attribute_mismatch"
      )
    }
    observed <- observed[names(predicted)]
  }
```

and replace everything from `results <- list()` to the end with:

```r
  compute_block_metrics(pred_flat$matrix, obs_flat$matrix,
                        pred_flat$block_map, metrics)
```

- [ ] **Step 5: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-metrics.R`, `test-tune_cca.R` cover univariate equivalence)

- [ ] **Step 6: Commit**

```bash
git add R/metrics.R tests/testthat/test-multivariate.R
git commit -m "per-variable and pooled spatial metrics"
```

---

### Task 8: Cross-validation with multivariate fields

**Files:**
- Modify: `R/tune_cca.R` (`evaluate_eof_fold`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("tune_eof handles multivariate input with per-variable metrics", {
  res <- tune_eof(prism_mv, k = 1:3, kfolds = 3, scale = TRUE,
                  weight = FALSE, n_reps = 2)
  expect_true(all(c("rmse", "rmse_tmean", "rmse_ppt") %in% names(res)))
  expect_equal(nrow(res), 9)
  s <- summarize_eof_cv(res)
  expect_true(attr(s, "best_k") %in% 1:3)
})

test_that("tune_cca downscales to a multivariate response", {
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))

  cv <- prep_cv_folds(coarse, prism_mv,
                      kfolds = 3, max_k_pred = 4, max_k_resp = 4,
                      scale_resp = TRUE, weight = FALSE)
  res <- tune_cca(cv, k_pred = 2:3, k_resp = 2:3)
  expect_true(all(c("rmse", "rmse_tmean", "rmse_ppt") %in% names(res)))
  s <- summarize_cv(res)
  expect_true(all(c("k_pred", "k_resp", "k_cca") %in% names(s)))
  expect_true("rmse_mean" %in% names(s))
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: `tune_cca` test may already PASS (its pieces landed in Tasks 4-7); `tune_eof` FAILS — `evaluate_eof_fold` flattens `[1]` only and indexes per-cell weights with concatenated `valid` indices.

- [ ] **Step 3: Rewrite `evaluate_eof_fold`**

Replace the function body in `R/tune_cca.R` (roxygen unchanged apart from noting per-variable metrics):

```r
evaluate_eof_fold <- function(fold, k, metrics, hidden_fraction = 0.2,
                              n_reps = 5, seed = 1L) {
  patterns_k <- fold$train_patterns[1:k]

  # Held-out anomalies in the space the patterns were fit in (own climatology)
  anomalies <- get_anomalies(fold$test_data, patterns_k$climatology,
                             scale = patterns_k$scaled, monthly = patterns_k$monthly)
  anom_mat <- flatten_time_space(units::drop_units(anomalies))$matrix

  valid <- patterns_k$valid_pixels
  n_valid <- length(valid)
  obs <- anom_mat[, valid, drop = FALSE]               # time x valid cells

  # Area weights applied during fitting must also weight the LS estimate;
  # one weight per grid cell, replicated across variable blocks
  n_vars <- length(patterns_k$names)
  w <- if (isTRUE(patterns_k$weight)) {
    rep(area_weights(fold$test_data), n_vars)[valid]
  } else {
    rep(1, n_valid)
  }

  block_map <- patterns_k$block_map
  if (is.null(block_map)) {
    block_map <- setNames(list(seq_len(ncol(anom_mat))), patterns_k$names[[1]])
  }

  weighted_loadings <- eof_loading_matrix(patterns_k)[valid, , drop = FALSE] * w
  weighted_obs <- sweep(obs, 2, w, `*`)                # time x valid cells

  n_hidden <- max(2L, round(hidden_fraction * n_valid))
  if (n_valid - n_hidden < k) {
    cli::cli_abort(
      "hidden_fraction = {hidden_fraction} leaves fewer than k = {k} visible cells.",
      class = "tidyeof_insufficient_cells"
    )
  }

  rep_metrics <- purrr::map(seq_len(n_reps), function(rep) {
    hidden <- with_seed(seed + fold$fold_id * 1000L + rep,
                        sample.int(n_valid, n_hidden))
    visible <- setdiff(seq_len(n_valid), hidden)

    # Estimate amplitudes from visible cells, predict the hidden ones
    amps <- t(qr.solve(weighted_loadings[visible, , drop = FALSE],
                       t(weighted_obs[, visible, drop = FALSE])))
    pred_hidden <- sweep(amps %*% t(weighted_loadings[hidden, , drop = FALSE]),
                         2, w[hidden], `/`)
    obs_hidden <- obs[, hidden, drop = FALSE]

    # Map hidden columns back to variable blocks for per-variable metrics
    hidden_cols <- valid[hidden]
    hidden_blocks <- purrr::map(block_map, ~which(hidden_cols %in% .x))
    hidden_blocks <- hidden_blocks[lengths(hidden_blocks) > 0]

    compute_block_metrics(pred_hidden, obs_hidden, hidden_blocks, metrics)
  })

  result <- tibble::tibble(fold = fold$fold_id)
  metric_names <- unique(unlist(purrr::map(rep_metrics, names)))
  for (m in metric_names) {
    vals <- vapply(rep_metrics,
                   function(v) if (is.null(v[[m]])) NA_real_ else v[[m]],
                   numeric(1))
    result[[m]] <- mean(vals, na.rm = TRUE)
  }
  result
}
```

(`evaluate_fold` for `tune_cca` needs NO change — `predict` and `compute_spatial_metrics` already handle multivariate, and its metric loop iterates `names(metric_values)`.)

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass (`test-tune_eof.R` covers univariate equivalence)

- [ ] **Step 5: Commit**

```bash
git add R/tune_cca.R tests/testthat/test-multivariate.R
git commit -m "multivariate cross-validation in tune_eof and tune_cca"
```

---

### Task 9: Coupling end-to-end + canonical patterns

**Files:**
- Modify: `R/couple.R` (`get_canonical_patterns`, lines ~414-488)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("CCA downscaling predicts a multivariate response end to end", {
  coarse <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.5), "°C"))

  pred_pat <- patterns(filter(coarse, time <= as.Date("2018-12-01")),
                       k = 3, weight = FALSE)
  resp_pat <- patterns(filter(prism_mv, time <= as.Date("2018-12-01")),
                       k = 3, scale = TRUE, weight = FALSE)

  coupled <- couple(pred_pat, resp_pat, k = 2)
  prediction <- predict(coupled, filter(coarse, time > as.Date("2018-12-01")))

  expect_s3_class(prediction, "stars")
  expect_named(prediction, c("tmean", "ppt"))
  expect_equal(units(prediction[["ppt"]]), units(prism_mv[["ppt"]]))
  expect_equal(unname(dim(prediction)), c(51, 51, 12))

  amps <- predict(coupled, filter(coarse, time > as.Date("2018-12-01")),
                  reconstruct = FALSE)
  expect_s3_class(amps, "tbl_df")
})

test_that("get_canonical_patterns returns all response variables", {
  coarse <- prism %>% mutate(tmean = tmean * 0.9)
  pred_pat <- patterns(coarse, k = 3, weight = FALSE)
  resp_pat <- patterns(prism_mv, k = 3, scale = TRUE, weight = FALSE)
  coupled <- couple(pred_pat, resp_pat, k = 2)

  cp <- get_canonical_patterns(coupled, type = "response")
  expect_named(cp, c("tmean", "ppt"))
  expect_equal(stars::st_get_dimension_values(cp, "CV"), c("CV1", "CV2"))

  cp_pred <- get_canonical_patterns(coupled, type = "predictor")
  expect_named(cp_pred, "tmean")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: the end-to-end test should PASS already (couple/predict are amplitude-space); `get_canonical_patterns` FAILS — it reads `eof_stars[[1]]` and returns one attribute.

- [ ] **Step 3: Rewrite the array section of `get_canonical_patterns`**

Replace everything in `get_canonical_patterns` from `# Extract EOF array (spatial dims + PC)` to the final `result` with:

```r
  # EOF attribute arrays are constructed spatial-first with PC last
  # (see get_eofs), one attribute per variable
  eof_stars <- patterns$eofs
  eof_dims <- stars::st_dimensions(eof_stars)
  spatial_dim_names <- setdiff(names(eof_dims), "PC")
  n_pcs <- length(stars::st_get_dimension_values(eof_stars, "PC"))

  new_dims <- eof_dims[spatial_dim_names]
  cv_dim <- list(
    from = 1L,
    to = k,
    offset = NA_real_,
    delta = NA_real_,
    refsys = NA_character_,
    point = FALSE,
    values = paste0("CV", 1:k)
  )
  class(cv_dim) <- "dimension"
  new_dims$CV <- cv_dim
  class(new_dims) <- "dimensions"

  result_list <- purrr::map(names(eof_stars), function(v) {
    eof_array <- eof_stars[[v]]
    spatial_shape <- dim(eof_array)[-length(dim(eof_array))]
    eof_matrix <- matrix(eof_array, nrow = prod(spatial_shape), ncol = n_pcs)

    # canonical_pattern[i] = sum_j EOF[j] * coef[j, i]
    canonical_array <- array(eof_matrix %*% coef_matrix,
                             dim = c(spatial_shape, k))
    setNames(stars::st_as_stars(canonical_array, dimensions = new_dims), v)
  })

  do.call(c, result_list)
```

- [ ] **Step 4: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add R/couple.R tests/testthat/test-multivariate.R
git commit -m "multivariate canonical patterns and coupling test"
```

---

### Task 10: Plot methods and teleconnection guards

**Files:**
- Modify: `R/plot_patterns.R` (`plot.patterns`, `.plot_eofs_internal`, `.plot_amplitudes_internal`)
- Modify: `R/teleconnections.R` (`get_correlation`, `get_fdr`: add `check_single_attribute(dat)`)
- Test: `tests/testthat/test-multivariate.R`

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-multivariate.R`:

```r
test_that("multivariate patterns plot with one panel row per variable", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  expect_s3_class(plot(pat), "patchwork")
  expect_s3_class(plot(pat, type = "eofs"), "patchwork")
  expect_s3_class(plot(pat, type = "amplitudes"), "ggplot")
})

test_that("unsupported multivariate plot modes error clearly", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  expect_error(plot(pat, type = "eofs", scaled = TRUE, rawdata = prism_mv),
               class = "tidyeof_multivariate_unsupported")
  expect_error(plot(pat, type = "amplitudes", scale = "raw"),
               class = "tidyeof_multivariate_unsupported")
})

test_that("teleconnection functions require a single-attribute field", {
  pat <- patterns(prism, k = 2)
  expect_error(get_correlation(prism_mv, pat),
               class = "tidyeof_multiple_attributes")
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'`
Expected: FAIL — `plot(pat, type = "eofs")` returns a single ggplot faceting only PC1..k of the first attribute; no `tidyeof_multivariate_unsupported` class exists.

- [ ] **Step 3: Rewrite `.plot_eofs_internal`**

Replace the function in `R/plot_patterns.R`:

```r
.plot_eofs_internal <- function(x, scaled = FALSE, rawdata = NULL, layout = NULL,
                                 overlay = NULL, overlay_color = "grey30", overlay_fill = NA) {
  facet_args <- if(!is.null(layout)) layout else list()

  overlay_layer <- if(!is.null(overlay)) {
    ggplot2::geom_sf(data = overlay, fill = overlay_fill, color = overlay_color, inherit.aes = FALSE)
  } else {
    NULL
  }

  if(scaled) {
    if (length(x$eofs) > 1) {
      cli::cli_abort(
        c(
          "Correlation (scaled) EOF maps are only supported for single-variable patterns.",
          "i" = "Fit patterns on one variable, or plot with {.code scaled = FALSE}."
        ),
        class = "tidyeof_multivariate_unsupported"
      )
    }
    if(is.null(rawdata)) {
      rlang::abort("rawdata must be provided when scaled = TRUE for correlation calculation", class = "tidyeof_missing_rawdata")
    }
    return(
      ggplot2::ggplot() +
        stars::geom_stars(data = get_correlation(rawdata, x)) +
        overlay_layer +
        do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
        ggplot2::scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1)) +
        ggplot2::coord_sf() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "right") +
        ggplot2::labs(fill = "Correlation")
    )
  }

  vars <- names(x$eofs)

  # One panel row per variable, each with its own fill scale: loadings of
  # different variables are not comparable on a shared color scale
  plot_one <- function(v) {
    ggplot2::ggplot() +
      stars::geom_stars(data = x$eofs[v]) +
      overlay_layer +
      do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
      scico::scale_fill_scico(palette = 'vik', midpoint = 0, na.value = NA) +
      ggplot2::coord_sf() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right") +
      ggplot2::labs(fill = if (length(vars) == 1) "Loading" else v)
  }

  if (length(vars) == 1) {
    return(plot_one(vars))
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package needed for multivariate EOF plots. Showing first variable only.")
    return(plot_one(vars[1]))
  }

  patchwork::wrap_plots(purrr::map(vars, plot_one), ncol = 1)
}
```

- [ ] **Step 4: Guard `.plot_amplitudes_internal` and scale combined heights**

In `.plot_amplitudes_internal`, at the top of the `scale == "raw"` branch add:

```r
    if (length(x$eofs) > 1) {
      cli::cli_abort(
        "Raw amplitude scaling mixes units across variables and is not supported for multivariate patterns. Use scale = 'standardized' or 'variance'.",
        class = "tidyeof_multivariate_unsupported"
      )
    }
```

In `plot.patterns`'s combined branch, give EOF rows more vertical space per variable — replace `heights = layout$heights` in the `wrap_plots` call with:

```r
heights = c(layout$heights[1] * length(x$eofs), layout$heights[2])
```

- [ ] **Step 5: Guard the teleconnection entry points**

In `R/teleconnections.R`, add `check_single_attribute(dat)` as the first line of both `get_correlation` and `get_fdr` (before the `amplitudes` default). Document in their roxygen: "For multivariate analyses, correlate one variable at a time (e.g. `get_correlation(dat["tmean"], pat)`)."

- [ ] **Step 6: Run the tests, then the full suite**

Run: `Rscript -e 'devtools::test(filter = "multivariate")'` — Expected: PASS
Run: `Rscript -e 'devtools::test()'` — Expected: all pass

- [ ] **Step 7: Commit**

```bash
git add R/plot_patterns.R R/teleconnections.R tests/testthat/test-multivariate.R
git commit -m "multivariate plotting and teleconnection guards"
```

---

### Task 11: Documentation, NEWS, and final verification

**Files:**
- Modify: roxygen blocks in `R/patterns_core.R`, `R/get_climatology.R`, `R/project_patterns.R`, `R/reconstruct_field.R`, `R/common_patterns.R`
- Modify: `NEWS.md`
- Regenerate: `man/`, `NAMESPACE` via `devtools::document()`

- [ ] **Step 1: Update roxygen docs**

`patterns()` `@param dat`, replace with:

```
#' @param dat A `stars` object containing spatial and temporal dimensions.
#'   Multiple attributes (e.g. temperature and precipitation on the same grid)
#'   are analyzed jointly as combined EOFs: each variable contributes a block
#'   of the space dimension, modes share one amplitude time series, and
#'   `scale = TRUE` is required so variables with different units contribute
#'   comparably.
```

and add to `patterns()` `@details`:

```
#' For multivariate input, per-pixel standardization is undefined where the
#' climatological standard deviation is ~0 (e.g. precipitation in arid cells);
#' such cells are dropped automatically like NA cells. Strongly skewed
#' variables such as precipitation often benefit from a sqrt or log transform
#' before analysis.
```

`common_patterns()`: add to its description: "Each source must currently contain a single variable; combining `common_patterns()` with multivariate (multi-attribute) input is not yet supported."

`project_patterns()`: document the attribute rule: "For multivariate patterns, `newdata` must contain the same variables as the training data (any order); univariate patterns accept any single-attribute object regardless of name."

`reconstruct()`: note the return is a multi-attribute stars object for multivariate patterns, with each variable's climatology and units restored.

`tune_cca()`/`tune_eof()` `@param metrics`: add "For multivariate fields each metric also gets per-variable columns (e.g. `rmse_tmean`); the plain name is the pooled score (sd-normalized RMS for `rmse`, mean for correlations)."

- [ ] **Step 2: Add NEWS entry**

At the top of `NEWS.md`:

```markdown
# tidyeof (development version)

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
```

- [ ] **Step 3: Regenerate docs and run the full check**

Run: `Rscript -e 'devtools::document()'` — Expected: man/ pages regenerate without roxygen errors
Run: `Rscript -e 'devtools::test()'` — Expected: all tests pass
Run: `Rscript -e 'devtools::check()'` — Expected: 0 errors, 0 warnings (notes at parity with the pre-change baseline; run `git stash && Rscript -e 'devtools::check()'` first if a baseline is needed, then `git stash pop`)

- [ ] **Step 4: Commit**

```bash
git add R/ man/ NAMESPACE NEWS.md
git commit -m "document multivariate EOF/CCA support"
```

---

## Out of scope (per spec)

- Variables on different grids (list-of-stars interface)
- Block normalization (`scale = FALSE` multivariate)
- Multivariate sources in `common_patterns()` (the `check_single_attribute` gate there is intentional)
- Vignette updates (working tree has uncommitted vignette edits; follow up separately)
