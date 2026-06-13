prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()

# Regression for the PC-name padding bug: names0() zero-pads by total count, so
# a fit at k >= 10 stores PC01.. while names0(5) gives PC1... project_patterns()
# used to regenerate names0(patterns$k), diverging from the truncated object's
# stored (training) names. The name-based centering reindex in
# apply_cca_prediction() then silently NA'd the predictor amplitudes, yielding
# all-NA predictions only for small k after truncation from a larger fit.

test_that("project_patterns columns match the object's stored amplitude names", {
  pp <- patterns(prism, k = 12)[1:5]
  proj <- project_patterns(pp, prism)
  expect_equal(setdiff(names(proj), "time"),
               setdiff(names(pp$amplitudes), "time"))
})

test_that("predict on a truncated CCA predictor is finite (not all-NA)", {
  pp <- patterns(prism, k = 12)
  rp <- patterns(prism, k = 12)
  coupled <- couple(pp[1:5], rp[1:5], k = 5, validate = FALSE)
  pred <- predict(coupled, prism)
  expect_true(all(is.finite(pred[[1]])))
})

test_that("predict on a truncated PCR predictor is finite (not all-NA)", {
  pp <- patterns(prism, k = 12)
  rp <- patterns(prism, k = 12)
  coupled <- couple(pp[1:5], rp[1:5], method = "pcr", validate = FALSE)
  pred <- predict(coupled, prism)
  expect_true(all(is.finite(pred[[1]])))
})

test_that("tune_cca yields finite scores for k_pred below max_k_pred", {
  cv <- suppressMessages(
    prep_cv_folds(prism, prism, kfolds = 3, max_k_pred = 10, max_k_resp = 10)
  )
  res <- suppressMessages(
    tune_cca(cv, k_pred = 4, k_resp = 6, k_cca = 3, metrics = "rmse")
  )
  expect_true(all(is.finite(res$rmse)))
})

test_that("centering that cannot be name-aligned errors instead of NA-ing", {
  # Defense in depth: a mismatch between the projected amplitude names and the
  # trained centering names must fail loudly, never silently corrupt to NA.
  new_amplitudes <- tibble::tibble(
    time = 1:3,
    PC1 = 1:3, PC2 = 1:3, PC3 = 1:3, PC4 = 1:3, PC5 = 1:3
  )
  cca_result <- list(
    xcenter = setNames(rep(0, 5), paste0("PC0", 1:5)),  # PC01..PC05
    xcoef = diag(5),
    cor = rep(0.5, 5),
    ycoef = diag(5),
    ycenter = FALSE
  )
  expect_error(
    tidyeof:::apply_cca_prediction(new_amplitudes, cca_result, k = 5),
    class = "tidyeof_centering_mismatch"
  )
})
