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
