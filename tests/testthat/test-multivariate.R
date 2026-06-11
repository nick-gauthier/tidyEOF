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
  expect_type(m, "double")
  expect_false(inherits(m, "units"))
})

test_that("multivariate anomalies equal per-variable anomalies (annual and monthly)", {
  for (use_monthly in c(FALSE, TRUE)) {
    anom_mv <- get_anomalies(prism_mv, scale = TRUE, monthly = use_monthly)
    anom_t <- get_anomalies(prism_mv["tmean"], scale = TRUE, monthly = use_monthly)
    anom_p <- get_anomalies(prism_mv["ppt"], scale = TRUE, monthly = use_monthly)
    expect_named(anom_mv, c("tmean", "ppt"))
    expect_equal(units::drop_units(anom_mv)[["tmean"]],
                 units::drop_units(anom_t)[[1]], tolerance = 1e-12)
    expect_equal(units::drop_units(anom_mv)[["ppt"]],
                 units::drop_units(anom_p)[[1]], tolerance = 1e-12)
  }
})

test_that("multivariate climatology round-trips through restore_climatology", {
  for (use_monthly in c(FALSE, TRUE)) {
    clim <- get_climatology(prism_mv, monthly = use_monthly)
    expect_named(clim$mean, c("tmean", "ppt"))
    anom <- get_anomalies(prism_mv, clim, scale = TRUE, monthly = use_monthly)
    restored <- restore_climatology(anom, clim, scale = TRUE, monthly = use_monthly)
    expect_equal(units::drop_units(restored[["tmean"]]),
                 units::drop_units(prism_mv[["tmean"]]), tolerance = 1e-8)
    expect_equal(units::drop_units(restored[["ppt"]]),
                 units::drop_units(prism_mv[["ppt"]]), tolerance = 1e-8)
  }
})

test_that("get_anomalies rejects climatology with mismatched attribute count", {
  clim <- get_climatology(prism_mv["tmean"])
  expect_error(get_anomalies(prism_mv, clim), class = "tidyeof_attribute_mismatch")

  anom <- get_anomalies(prism_mv, scale = TRUE, monthly = TRUE)
  clim_uni_monthly <- get_climatology(prism_mv["tmean"], monthly = TRUE)
  expect_error(restore_climatology(anom, clim_uni_monthly, scale = TRUE, monthly = TRUE),
               class = "tidyeof_attribute_mismatch")
})
