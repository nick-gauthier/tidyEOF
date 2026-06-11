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
