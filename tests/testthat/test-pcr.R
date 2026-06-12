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

test_that("predict warns when k is passed to a PCR object", {
  pred <- patterns(prism, k = 4, weight = FALSE)
  resp <- patterns(prism, k = 3, weight = FALSE)
  cpl <- couple(pred, resp, method = "pcr")
  expect_warning(predict(cpl, prism, k = 2, reconstruct = FALSE),
                 class = "tidyeof_k_ignored")
  # default k = NULL does not warn
  expect_no_warning(predict(cpl, prism, reconstruct = FALSE))
})

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
