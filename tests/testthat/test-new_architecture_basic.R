test_that("basic architecture functions work", {
  # Test that the new functions exist and have correct signatures
  expect_true(exists("couple"))
  expect_true(exists("predict.coupled_patterns"))
  expect_true(exists("project_patterns"))

  # Test that coupled_patterns class methods exist
  expect_true(exists("print.coupled_patterns"))
  expect_true(exists("summary.coupled_patterns"))

  # NOTE: Legacy wrappers (predict_cca, predict_gam, etc.) have been archived
  # See archive/legacy_wrappers.R for reference implementation
})

test_that("couple validation works", {
  # Test input validation
  expect_error(
    couple(NULL, NULL),
    "predictor_patterns must be a patterns object"
  )

  expect_error(
    couple(list(), list()),
    "predictor_patterns must be a patterns object"
  )

  # Test method validation
  expect_error(
    couple(data.frame(time = 1:10), data.frame(time = 1:10), method = "invalid"),
    "Only 'cca' method is supported"
  )
})

test_that("canonical correlation utilities exist", {
  expect_true(exists("get_canonical_correlations"))
  expect_true(exists("get_canonical_variables"))
  expect_true(exists("apply_cca_prediction"))
})

# NOTE: Legacy wrapper deprecation tests removed - functions have been archived
# See archive/legacy_wrappers.R for reference implementation

test_that("helper functions exist", {
  expect_true(exists("has_geometry_dimension"))
  expect_true(exists("extract_amplitudes_matrix"))
  expect_true(exists("validate_patterns_compatibility"))
})

test_that("simple matrix-based CCA works", {
  # Create simple test matrices
  set.seed(123)
  n <- 20
  predictor_matrix <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  response_matrix <- matrix(rnorm(n * 2), nrow = n, ncol = 2)

  # Test that cancor works (this is what our implementation uses)
  cca_result <- cancor(predictor_matrix, response_matrix)

  expect_true(is.list(cca_result))
  expect_true("cor" %in% names(cca_result))
  expect_true("xcoef" %in% names(cca_result))
  expect_true("ycoef" %in% names(cca_result))
  expect_equal(length(cca_result$cor), 2)  # min(3, 2) = 2
})

test_that("amplitude matrix extraction works", {
  # Test the helper function with simple data frames
  test_df1 <- data.frame(time = 1:5, PC1 = 1:5, PC2 = 6:10)
  test_df2 <- list(amplitudes = test_df1)
  class(test_df2) <- "patterns"

  # Test with data frame
  result1 <- extract_amplitudes_matrix(test_df1)
  expect_true(is.matrix(result1))
  expect_equal(ncol(result1), 2)  # PC1 and PC2
  expect_equal(nrow(result1), 5)  # 5 time steps

  # Test with patterns object
  result2 <- extract_amplitudes_matrix(test_df2)
  expect_true(is.matrix(result2))
  expect_equal(ncol(result2), 2)
  expect_equal(nrow(result2), 5)
})

test_that("CCA centering is honored in prediction and canonical variables", {
  set.seed(42)
  n <- 30
  predictor_matrix <- matrix(rnorm(n * 3), nrow = n, ncol = 3) + 2
  response_matrix <- matrix(rnorm(n * 2), nrow = n, ncol = 2) - 1
  colnames(predictor_matrix) <- paste0("PC", 1:3)
  colnames(response_matrix) <- paste0("PC", 1:2)

  cca_result <- cancor(predictor_matrix, response_matrix, xcenter = TRUE, ycenter = TRUE)

  new_amplitudes <- data.frame(time = seq_len(n), predictor_matrix, check.names = FALSE)

  predicted <- apply_cca_prediction(new_amplitudes, cca_result, k = 2)

  pred_centered <- sweep(predictor_matrix, 2, cca_result$xcenter, "-")
  canonical_predictors <- pred_centered %*% cca_result$xcoef[, 1:2, drop = FALSE]
  canonical_responses <- canonical_predictors %*% diag(cca_result$cor[1:2], nrow = 2)
  ycoef_inv <- if (nrow(cca_result$ycoef) == ncol(cca_result$ycoef)) {
    solve(cca_result$ycoef)
  } else {
    MASS::ginv(cca_result$ycoef)
  }
  expected_resp <- canonical_responses %*% ycoef_inv[1:2, , drop = FALSE]
  expected_resp <- sweep(expected_resp, 2, cca_result$ycenter, "+")
  colnames(expected_resp) <- paste0("PC", 1:ncol(expected_resp))

  expect_equal(as.matrix(predicted[, -1]), expected_resp, tolerance = 1e-8)

  coupled <- list(cca = cca_result, k = 2)
  canonical_vars <- get_canonical_variables(coupled, data = new_amplitudes, type = "predictor", k = 2)
  expected_cv <- pred_centered %*% cca_result$xcoef[, 1:2, drop = FALSE]
  colnames(expected_cv) <- paste0("CV", 1:2)

  expect_equal(as.matrix(canonical_vars[, -1]), expected_cv, tolerance = 1e-8)
})

test_that("CCA truncation uses leading rows of the full ycoef inverse", {
  set.seed(99)
  n <- 25
  predictor_matrix <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  response_matrix <- matrix(rnorm(n * 2), nrow = n, ncol = 2)
  colnames(predictor_matrix) <- paste0("PC", 1:3)
  colnames(response_matrix) <- paste0("PC", 1:2)

  cca_result <- cancor(predictor_matrix, response_matrix, xcenter = TRUE, ycenter = TRUE)
  k <- 1

  new_amplitudes <- data.frame(time = seq_len(n), predictor_matrix, check.names = FALSE)
  predicted <- apply_cca_prediction(new_amplitudes, cca_result, k = k)

  pred_centered <- sweep(predictor_matrix, 2, cca_result$xcenter, "-")
  canonical_predictors <- pred_centered %*% cca_result$xcoef[, 1:k, drop = FALSE]
  canonical_responses <- canonical_predictors %*% diag(cca_result$cor[1:k], nrow = k)

  # Regression of response PCs on the retained canonical variates =
  # leading k rows of the full ycoef (pseudo)inverse
  expected_resp <- canonical_responses %*%
    MASS::ginv(cca_result$ycoef)[1:k, , drop = FALSE]
  expected_resp <- sweep(expected_resp, 2, cca_result$ycenter, "+")

  expect_equal(unname(as.matrix(predicted[, -1])), unname(expected_resp), tolerance = 1e-8)
})
