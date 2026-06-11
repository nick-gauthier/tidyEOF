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
  expect_true(all(is.finite(units::drop_units(rec[["tmean"]]))))
  expect_true(all(is.finite(units::drop_units(rec[["ppt"]]))))
})

test_that("eof_loading_matrix stacks variable blocks", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  m <- tidyeof:::eof_loading_matrix(pat)
  expect_equal(dim(m), c(2 * 2601, 3))
  expect_equal(m[1:2601, ], matrix(pat$eofs[["tmean"]], nrow = 2601, ncol = 3))
})

test_that("multivariate reconstruction preserves per-variable NA masks", {
  z <- prism_mv
  arr <- units::drop_units(z[["ppt"]])
  arr[2, 1, ] <- NA          # mask one ppt cell across all times
  arr[3, 5, ] <- NA          # and a second ppt cell
  z$ppt <- units::set_units(arr, "mm")

  pat <- patterns(z, k = 30, scale = TRUE, weight = FALSE)
  rec <- reconstruct(pat)

  # tmean has no masked cells; ppt's two masked cells must come back NA
  expect_false(any(is.na(units::drop_units(rec[["tmean"]]))))
  rec_ppt <- units::drop_units(rec[["ppt"]])
  expect_true(all(is.na(rec_ppt[2, 1, ])))
  expect_true(all(is.na(rec_ppt[3, 5, ])))
  # all other ppt cells reconstructed finite
  other <- rec_ppt
  other[2, 1, ] <- 0; other[3, 5, ] <- 0
  expect_true(all(is.finite(other)))
})

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

test_that("compute_spatial_metrics rejects mismatched attribute sets", {
  pat <- patterns(prism_mv, k = 3, scale = TRUE)
  rec <- reconstruct(pat)
  obs_bad <- setNames(prism_mv, c("tmean", "precip"))
  expect_error(tidyeof:::compute_spatial_metrics(rec, obs_bad),
               class = "tidyeof_attribute_mismatch")
})

test_that("multivariate metric subset returns only requested metric columns", {
  pat <- patterns(prism_mv, k = 4, scale = TRUE)
  rec <- reconstruct(pat)
  m <- tidyeof:::compute_spatial_metrics(rec, prism_mv, metrics = "rmse")
  expect_named(m, c("rmse", "rmse_tmean", "rmse_ppt"))
  expect_false(any(grepl("cor_", names(m))))
})
