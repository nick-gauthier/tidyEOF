prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()

# Helper: load test data and create mock second source
load_test_data <- function() {
  prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
    readRDS()

  # Mock second source: scaled version + small noise
  set.seed(42)
  mock2 <- prism %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 0.2), "\u00b0C"))

  list(prism = prism, mock2 = mock2)
}

# --- Construction tests ---

test_that("common_patterns returns correct S3 class", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  expect_s3_class(cpat, "common_patterns")
  expect_equal(cpat$sources, c("a", "b"))
  expect_equal(cpat$k, 3)
  expect_equal(length(cpat$source_patterns), 2)
})

test_that("source-specific patterns are extractable with $ and [[", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  # $ accessor
  pat_a <- cpat$a
  expect_s3_class(pat_a, "patterns")
  expect_equal(pat_a$k, 3)

  # [[ accessor
  pat_b <- cpat[["b"]]
  expect_s3_class(pat_b, "patterns")
  expect_equal(pat_b$k, 3)

  # Non-source names fall through to standard list access
  expect_equal(cpat$k, 3)
  expect_equal(cpat[["k"]], 3)
})

# --- Shared structure tests ---

test_that("sources share EOFs and proj_matrix", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  expect_identical(cpat$a$eofs, cpat$b$eofs)
  expect_identical(cpat$a$proj_matrix, cpat$b$proj_matrix)
  expect_identical(cpat$a$valid_pixels, cpat$b$valid_pixels)
  expect_identical(cpat$a$eigenvalues, cpat$b$eigenvalues)
})

test_that("sources have different climatologies", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  # Climatologies should differ since data differs
  expect_false(identical(cpat$a$climatology, cpat$b$climatology))
})

test_that("sources have different amplitudes", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  # Amplitudes should differ (different anomalies projected onto same EOFs)
  expect_false(identical(cpat$a$amplitudes, cpat$b$amplitudes))

  # But time columns should match original sources
  a_times <- stars::st_get_dimension_values(d$prism, "time")
  b_times <- stars::st_get_dimension_values(d$mock2, "time")
  expect_equal(cpat$a$amplitudes$time, a_times)
  expect_equal(cpat$b$amplitudes$time, b_times)
})

# --- Integration with couple/predict ---

test_that("common patterns work with couple and predict", {
  d <- load_test_data()

  # Split into train/test
  all_times <- stars::st_get_dimension_values(d$prism, "time")
  n <- length(all_times)
  split_idx <- floor(n * 0.7)
  train_times <- all_times[1:split_idx]
  test_times <- all_times[(split_idx + 1):n]

  train_a <- dplyr::filter(d$prism, time %in% train_times)
  train_b <- dplyr::filter(d$mock2, time %in% train_times)
  test_a <- dplyr::filter(d$prism, time %in% test_times)

  cpat <- common_patterns(
    list(a = train_a, b = train_b),
    k = 3, weight = FALSE
  )

  # Use source a patterns as both predictor and response (self-downscaling test)
  fine_pat <- patterns(train_a, k = 3, weight = FALSE)
  coupled <- couple(cpat$a, fine_pat, k = 2)

  # Predict with same source
  pred <- predict(coupled, test_a)
  expect_s3_class(pred, "stars")

  # Predict without reconstruction
  pred_amps <- predict(coupled, test_a, reconstruct = FALSE)
  expect_s3_class(pred_amps, "tbl_df")
  expect_true("time" %in% names(pred_amps))
})

test_that("cross-source prediction works with predictor_patterns", {
  d <- load_test_data()

  all_times <- stars::st_get_dimension_values(d$prism, "time")
  n <- length(all_times)
  split_idx <- floor(n * 0.7)
  train_times <- all_times[1:split_idx]
  test_times <- all_times[(split_idx + 1):n]

  train_a <- dplyr::filter(d$prism, time %in% train_times)
  train_b <- dplyr::filter(d$mock2, time %in% train_times)
  test_b <- dplyr::filter(d$mock2, time %in% test_times)

  cpat <- common_patterns(
    list(a = train_a, b = train_b),
    k = 3, weight = FALSE
  )

  fine_pat <- patterns(train_a, k = 3, weight = FALSE)
  coupled <- couple(cpat$a, fine_pat, k = 2)

  # Cross-source prediction: use source b data with source b patterns

  pred <- predict(coupled, test_b, predictor_patterns = cpat$b)
  expect_s3_class(pred, "stars")
})

test_that("predictor_patterns validation catches incompatible dimensions", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  fine_pat <- patterns(d$prism, k = 3, weight = FALSE)
  coupled <- couple(cpat$a, fine_pat, k = 2)

  # Create incompatible patterns with different k
  bad_pat <- patterns(d$prism, k = 5, weight = FALSE)

  expect_error(
    predict(coupled, d$prism, predictor_patterns = bad_pat),
    "predictor_patterns"
  )
})

# --- CV with common_with ---

test_that("prep_cv_folds works with common_with", {
  d <- load_test_data()

  cv <- prep_cv_folds(
    d$prism, d$prism,
    common_with = list(mock = d$mock2),
    kfolds = 2,
    max_k_pred = 3,
    max_k_resp = 3,
    weight = FALSE
  )

  expect_s3_class(cv, "cv_folds")
  expect_equal(cv$common_with_sources, "mock")

  # Folds should contain valid patterns
  fold1 <- cv$folds[[1]]
  expect_s3_class(fold1$train_pred_patterns, "patterns")
  expect_s3_class(fold1$train_resp_patterns, "patterns")
  expect_equal(fold1$train_pred_patterns$k, 3)
})

test_that("CV with common_with runs through tune_cca", {
  d <- load_test_data()

  cv <- prep_cv_folds(
    d$prism, d$prism,
    common_with = list(mock = d$mock2),
    kfolds = 2,
    max_k_pred = 3,
    max_k_resp = 3,
    weight = FALSE
  )

  results <- tune_cca(cv,
                      k_pred = 2:3,
                      k_resp = 2:3,
                      metrics = "rmse")

  expect_s3_class(results, "tbl_df")
  expect_true(nrow(results) > 0)
  expect_true("rmse" %in% names(results))
  expect_true(all(results$rmse > 0))
})

# --- Valid pixels intersection ---

test_that("valid pixels intersection handles different NA patterns", {
  d <- load_test_data()

  # Introduce NAs in different locations for each source
  prism_na <- d$prism
  mock_na <- d$mock2

  # Set some pixels to NA in source a only
  arr_a <- prism_na[[1]]
  arr_a[1, 1, ] <- NA
  prism_na[[1]] <- arr_a

  # Set different pixels to NA in source b only
  arr_b <- mock_na[[1]]
  arr_b[2, 2, ] <- NA
  mock_na[[1]] <- arr_b

  cpat <- common_patterns(
    list(a = prism_na, b = mock_na),
    k = 3, weight = FALSE
  )

  # Both sources should have the same valid_pixels (intersection)
  expect_identical(cpat$a$valid_pixels, cpat$b$valid_pixels)

  # The intersection should exclude both NA locations
  # Valid pixels from a alone would include pixel at (2,2)
  # Valid pixels from b alone would include pixel at (1,1)
  # Intersection excludes both
  vp <- cpat$a$valid_pixels
  expect_true(length(vp) > 0)
})

# --- Single source ---

test_that("single source common_patterns is equivalent to patterns", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(only = d$prism),
    k = 3, weight = FALSE
  )

  expect_s3_class(cpat, "common_patterns")
  expect_equal(cpat$sources, "only")

  pat <- cpat$only
  expect_s3_class(pat, "patterns")
  expect_equal(pat$k, 3)

  # Direct patterns() call for comparison (scale = TRUE to match common_patterns default)
  direct <- patterns(d$prism, k = 3, weight = FALSE, scale = TRUE)

  # Eigenvalues should be very close (same PCA)
  expect_equal(
    pat$eigenvalues$eigenvalues[1:3],
    direct$eigenvalues$eigenvalues[1:3],
    tolerance = 1e-6
  )

  # Amplitudes should correlate very highly (up to sign)
  for (i in 1:3) {
    pc <- paste0("PC", i)
    r <- cor(pat$amplitudes[[pc]], direct$amplitudes[[pc]])
    expect_true(abs(r) > 0.999)
  }
})

# --- Input validation ---

test_that("common_patterns rejects unnamed list", {
  d <- load_test_data()

  expect_error(
    common_patterns(list(d$prism, d$mock2), k = 3),
    "named list"
  )
})

test_that("common_patterns rejects empty list", {
  expect_error(
    common_patterns(list(), k = 3),
    "non-empty"
  )
})

test_that("common_patterns rejects non-stars elements", {
  expect_error(
    common_patterns(list(a = data.frame(x = 1:10)), k = 3),
    "stars"
  )
})

test_that("common_patterns rejects mismatched spatial grids", {
  d <- load_test_data()

  # Create a dataset with different spatial coordinates (shifted x offset)
  bad_grid <- stars::st_set_dimensions(d$prism, "x", offset = -100)

  expect_error(
    common_patterns(list(a = d$prism, b = bad_grid), k = 3),
    "grid_mismatch|coordinates differ"
  )
})

# --- Print method ---

test_that("print.common_patterns works without error", {
  d <- load_test_data()

  cpat <- common_patterns(
    list(a = d$prism, b = d$mock2),
    k = 3, weight = FALSE
  )

  expect_no_error(print(cpat))
})

# --- Different time ranges ---

test_that("sources with different time ranges work", {
  d <- load_test_data()

  all_times <- stars::st_get_dimension_values(d$prism, "time")
  n <- length(all_times)

  # Source a: first 80% of times
  src_a <- dplyr::filter(d$prism, time %in% all_times[1:floor(n * 0.8)])
  # Source b: last 80% of times (overlapping middle)
  src_b <- dplyr::filter(d$mock2, time %in% all_times[ceiling(n * 0.2):n])

  cpat <- common_patterns(
    list(a = src_a, b = src_b),
    k = 3, weight = FALSE
  )

  expect_s3_class(cpat, "common_patterns")

  # Each source should have its own number of time steps
  n_a <- nrow(cpat$a$amplitudes)
  n_b <- nrow(cpat$b$amplitudes)
  expect_equal(n_a, floor(n * 0.8))
  expect_equal(n_b, n - ceiling(n * 0.2) + 1)
})

test_that("reconstruct names a common_patterns source by its own variable, not the reference", {
  era <- prism                       # attribute "tmean"
  phyda <- setNames(prism, "temp")   # same grid/times, attribute "temp"

  cpat <- common_patterns(list(era = era, phyda = phyda),
                          k = 3, scale = TRUE, weight = FALSE)

  rec <- reconstruct(cpat$phyda)
  expect_named(rec, "temp")
  expect_equal(units(rec[["temp"]]), units(prism[["tmean"]]))

  # The reference source still reconstructs correctly too
  rec_era <- reconstruct(cpat$era)
  expect_named(rec_era, "tmean")
})
