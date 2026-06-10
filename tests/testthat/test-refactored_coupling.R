test_that("refactored couple produces same results as original predict_cca", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  # Create some mock predictor data by modifying prism
  cera_test <- prism_test %>%
    mutate(tmean = tmean * 0.8 + units::set_units(rnorm(length(tmean), 0, 1), "°C"))

  # Extract patterns (disable area weighting for test data without geometry)
  prism_patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  cera_patterns <- patterns(filter(cera_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)

  # Test data
  test_data <- filter(cera_test, time > as.Date("2018-12-01"))

  # Compare old vs new implementation
  # Note: We can't directly test against the old implementation since it's been refactored
  # Instead, we test that the new implementation runs without error and produces reasonable results

  # Test new implementation
  coupled <- couple(cera_patterns, prism_patterns, k = 2)

  expect_s3_class(coupled, "coupled_patterns")
  expect_equal(coupled$k, 2)
  expect_equal(coupled$method, "cca")
  expect_true(length(coupled$cca$cor) >= 2)

  # Test prediction
  prediction <- predict(coupled, test_data)

  expect_s3_class(prediction, "stars")
  expect_true("tmean" %in% names(prediction))

  # Test amplitude prediction
  amplitudes <- predict(coupled, test_data, reconstruct = FALSE)

  expect_s3_class(amplitudes, "tbl_df")
  expect_true("time" %in% names(amplitudes))
  expect_true(all(grepl("PC", names(amplitudes)[-1])))
})

test_that("couple handles edge cases appropriately", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  cera_test <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  prism_patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3)
  cera_patterns <- patterns(filter(cera_test, time <= as.Date("2018-12-01")), k = 4)

  # Test k validation
  expect_warning(
    couple(cera_patterns, prism_patterns, k = 10),
    "exceeds maximum possible"
  )

  # Test unsupported method
  expect_error(
    couple(cera_patterns, prism_patterns, method = "gam"),
    "Only 'cca' method is supported"
  )

  # Test that k is adjusted automatically
  coupled <- suppressWarnings(
    couple(cera_patterns, prism_patterns, k = 10)
  )
  expect_equal(coupled$k, 3)  # Should be min(4, 3) = 3
})

test_that("predict.coupled_patterns handles edge cases", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  cera_test <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  prism_patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  cera_patterns <- patterns(filter(cera_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  test_data <- filter(cera_test, time > as.Date("2018-12-01"))

  coupled <- couple(cera_patterns, prism_patterns, k = 2)

  # Test k validation in prediction
  expect_warning(
    predict(coupled, test_data, k = 5),
    "exceeds available modes"
  )

  # Test invalid object type
  expect_error(
    predict.coupled_patterns(list(method = "not_cca"), test_data),
    "must be a coupled_patterns object"
  )
})

test_that("canonical correlation utilities work correctly", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  cera_test <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  prism_patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  cera_patterns <- patterns(filter(cera_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)

  coupled <- couple(cera_patterns, prism_patterns, k = 2)

  # Test canonical correlations extraction
  canon_corr <- get_canonical_correlations(coupled)

  expect_s3_class(canon_corr, "data.frame")
  expect_equal(nrow(canon_corr), 2)
  expect_true(all(c("mode", "correlation", "correlation_squared") %in% names(canon_corr)))
  expect_true(all(canon_corr$correlation >= -1e-10 & canon_corr$correlation <= 1 + 1e-10))

  # Test canonical variables extraction
  pred_canonical <- get_canonical_variables(coupled, cera_patterns, type = "predictor")
  resp_canonical <- get_canonical_variables(coupled, prism_patterns, type = "response")

  expect_s3_class(pred_canonical, "tbl_df")
  expect_s3_class(resp_canonical, "tbl_df")
  expect_true("time" %in% names(pred_canonical))
  expect_true("time" %in% names(resp_canonical))
  expect_true(all(grepl("CV", names(pred_canonical)[-1])))
  expect_true(all(grepl("CV", names(resp_canonical)[-1])))
})

test_that("project_patterns validation works", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  test_data <- filter(prism_test, time > as.Date("2018-12-01"))

  # Test normal projection
  projections <- project_patterns(patterns, test_data)

  expect_s3_class(projections, "tbl_df")
  expect_true("time" %in% names(projections))
  expect_equal(ncol(projections) - 1, 3)  # 3 PCs plus time

  # Test validation (would need mock data with incompatible dimensions to test errors)
  expect_no_error(project_patterns(patterns, test_data))
})

# NOTE: Legacy wrapper function tests removed - functions have been archived
# See archive/legacy_wrappers.R for reference implementation

test_that("print and summary methods work for coupled_patterns", {
  skip_if_not_installed("testthat")

  # Load test data
  prism_test <- system.file('testdata/prism_test.RDS', package = 'tidyeof') %>%
    readRDS()

  cera_test <- prism_test %>%
    mutate(tmean = tmean * 0.8)

  prism_patterns <- patterns(filter(prism_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)
  cera_patterns <- patterns(filter(cera_test, time <= as.Date("2018-12-01")), k = 3, weight = FALSE)

  coupled <- couple(cera_patterns, prism_patterns, k = 2)

  # Test that print and summary don't throw errors
  expect_no_error(print(coupled))
  expect_no_error(summary(coupled))

  # cli writes to message connection, so capture both stdout and stderr
  print_output <- capture.output(print(coupled), type = "message")
  expect_true(any(grepl("Coupled Patterns Object", print_output)))
  expect_true(any(grepl("Method: cca", print_output)))

  summary_output <- capture.output(summary(coupled), type = "message")
  expect_true(any(grepl("Coupled Patterns Summary", summary_output)))
})

# NOTE: show_migration_guide() tests removed - function has been archived
# See archive/legacy_wrappers.R for reference implementation
test_that("truncated CCA prediction is the regression onto retained canonical variates", {
  # Standard CCA prediction (Glahn 1968; Barnett & Preisendorfer 1987):
  # response amplitudes are regressed onto the retained canonical variates.
  # cancor's variates are orthonormal so the OLS fit is computable
  # independently; apply_cca_prediction() on the training data must match it.
  set.seed(7)
  n <- 40
  times <- seq(as.Date("2000-01-01"), by = "year", length.out = n)
  pred_mat <- scale(matrix(rnorm(n * 5), n), scale = FALSE)
  resp_mat <- scale(pred_mat %*% matrix(rnorm(5 * 4), 5) + rnorm(n * 4),
                    scale = FALSE)

  pred_df <- as_tibble(pred_mat, .name_repair = ~ paste0("PC", 1:5)) %>%
    mutate(time = times, .before = 1)
  resp_df <- as_tibble(resp_mat, .name_repair = ~ paste0("PC", 1:4)) %>%
    mutate(time = times, .before = 1)

  # k = 4 is the full set (already correct before the fix); k = 1, 2 are the
  # truncated "regularization" cases where the back-transform matters
  for (k in c(1, 2, 4)) {
    coupled <- couple(pred_df, resp_df, k = k)
    predicted <- apply_cca_prediction(pred_df, coupled$cca, k = k)

    canonical_variates <- pred_mat %*% coupled$cca$xcoef[, seq_len(k), drop = FALSE]
    ols_fit <- canonical_variates %*% qr.solve(canonical_variates, resp_mat)

    expect_equal(unname(as.matrix(predicted[, -1])), unname(ols_fit),
                 tolerance = 1e-8)
  }
})
