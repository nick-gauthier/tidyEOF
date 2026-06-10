# Tests for geometry-dimension (sf-based) stars support
library(testthat)
library(stars)
library(sf)

# Helper: build a small geometry-based stars object (geometry × time)
make_geometry_stars <- function(n_geom = 20, n_time = 36, seed = 42) {
  set.seed(seed)
  polys <- st_sfc(lapply(seq_len(n_geom), function(i)
    st_polygon(list(matrix(c(i, 0, i + 1, 0, i + 1, 1, i, 1, i, 0),
                           ncol = 2, byrow = TRUE)))))
  polys <- st_set_crs(polys, 4326)
  times <- seq(as.Date("2000-01-01"), length.out = n_time, by = "month")
  arr <- array(rnorm(n_geom * n_time, mean = 20, sd = 5),
               dim = c(n_geom, n_time))
  st_as_stars(list(values = arr),
              dimensions = st_dimensions(geometry = polys, time = times))
}

# ---- get_climatology (annual) ------------------------------------------------

test_that("get_climatology annual works for geometry stars", {
  dat <- make_geometry_stars()
  clim <- get_climatology(dat)

  expect_type(clim, "list")
  expect_named(clim, c("mean", "sd"))
  expect_s3_class(clim$mean, "stars")
  expect_s3_class(clim$sd, "stars")

  # Should have geometry dimension but no time
  expect_true("geometry" %in% names(st_dimensions(clim$mean)))
  expect_false("time" %in% names(st_dimensions(clim$mean)))

  # Values should match manual calculation
  manual_mean <- rowMeans(dat[[1]])
  expect_equal(as.numeric(clim$mean[[1]]), manual_mean, tolerance = 1e-10)
})

# ---- get_climatology (monthly) -----------------------------------------------

test_that("get_climatology monthly works for geometry stars", {
  dat <- make_geometry_stars()
  clim <- get_climatology(dat, monthly = TRUE)

  expect_type(clim, "list")
  expect_s3_class(clim$mean, "stars")
  expect_s3_class(clim$sd, "stars")

  # Should have geometry + month dimensions
  dims <- st_dimensions(clim$mean)
  expect_true("geometry" %in% names(dims))
  expect_true("month" %in% names(dims))
  expect_false("time" %in% names(dims))

  # 12 months
  expect_equal(dim(clim$mean)[["month"]], 12)

  # Verify a known month: January values
  times <- seq(as.Date("2000-01-01"), length.out = 36, by = "month")
  jan_idx <- which(format(times, "%B") == "January")
  manual_jan_mean <- rowMeans(dat[[1]][, jan_idx])
  expect_equal(as.numeric(clim$mean[[1]][, 1]), manual_jan_mean, tolerance = 1e-10)
})

# ---- get_anomalies + restore_climatology roundtrip ---------------------------

test_that("get_anomalies and restore_climatology roundtrip (annual) for geometry stars", {
  dat <- make_geometry_stars()

  for (sc in c(FALSE, TRUE)) {
    clim <- get_climatology(dat)
    anom <- get_anomalies(dat, scale = sc)
    restored <- restore_climatology(anom, clim, scale = sc)
    expect_equal(dat[[1]], restored[[1]], tolerance = 1e-10,
                 label = paste("annual, scale =", sc))
  }
})

test_that("get_anomalies and restore_climatology roundtrip (monthly) for geometry stars", {
  dat <- make_geometry_stars()

  for (sc in c(FALSE, TRUE)) {
    clim <- get_climatology(dat, monthly = TRUE)
    anom <- get_anomalies(dat, monthly = TRUE, scale = sc)
    restored <- restore_climatology(anom, clim, monthly = TRUE, scale = sc)
    expect_equal(dat[[1]], restored[[1]], tolerance = 1e-10,
                 label = paste("monthly, scale =", sc))
  }
})

test_that("get_anomalies monthly preserves dimensions for geometry stars", {
  dat <- make_geometry_stars()
  anom <- get_anomalies(dat, monthly = TRUE)

  expect_equal(dim(anom), dim(dat))
  expect_equal(st_get_dimension_values(anom, "time"),
               st_get_dimension_values(dat, "time"))
  expect_true("geometry" %in% names(st_dimensions(anom)))
})

# ---- patterns() on geometry stars --------------------------------------------

test_that("patterns() works on geometry stars (annual)", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 3)

  expect_s3_class(pat, "patterns")
  expect_equal(pat$k, 3)
  expect_s3_class(pat$eofs, "stars")
  expect_true("geometry" %in% names(st_dimensions(pat$eofs)))
  expect_true("PC" %in% names(st_dimensions(pat$eofs)))
})

test_that("patterns() works on geometry stars (monthly)", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 3, monthly = TRUE)

  expect_s3_class(pat, "patterns")
  expect_true(pat$monthly)
})

# ---- reconstruct() on geometry stars -----------------------------------------

test_that("reconstruct() roundtrip on geometry stars", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 5)
  recon <- reconstruct(pat)

  # Reconstructed should have same dimensions as original
  expect_equal(dim(recon), dim(dat))
  expect_true("geometry" %in% names(st_dimensions(recon)))
  expect_true("time" %in% names(st_dimensions(recon)))
})

test_that("reconstruct() roundtrip on geometry stars (monthly)", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 5, monthly = TRUE)
  recon <- reconstruct(pat)

  expect_equal(dim(recon), dim(dat))
})

# ---- get_correlation() on geometry stars -------------------------------------

test_that("get_correlation() works on geometry stars", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 3)
  cors <- get_correlation(dat, pat)

  expect_s3_class(cors, "stars")
  # Should have geometry + PC dimensions
  dims <- st_dimensions(cors)
  expect_true("geometry" %in% names(dims))
  expect_true("PC" %in% names(dims))
  expect_equal(dim(cors)[["PC"]], 3)

  # Correlation values should be between -1 and 1
  vals <- cors[[1]]
  expect_true(all(abs(vals[!is.na(vals)]) <= 1))
})

# ---- get_fdr() errors for geometry stars -------------------------------------

test_that("get_fdr() errors clearly for geometry stars", {
  dat <- make_geometry_stars()
  pat <- patterns(dat, k = 3)
  expect_error(get_fdr(dat, pat), "raster data")
})

# ---- project_patterns() on geometry stars ------------------------------------

test_that("project_patterns() works on geometry stars", {
  dat <- make_geometry_stars(n_time = 48)
  # Use first 36 months for training, last 12 for projection
  train <- dat[, , 1:36]
  test_dat <- dat[, , 37:48]

  pat <- patterns(train, k = 3)
  proj <- project_patterns(pat, test_dat)

  expect_s3_class(proj, "tbl_df")
  expect_true("time" %in% names(proj))
  expect_equal(nrow(proj), 12)
  expect_equal(ncol(proj), 4)  # time + 3 PCs
})

# ---- area weights guard for zero-area (point) geometries ---------------------

make_point_stars <- function(n_geom = 20, n_time = 24, seed = 7) {
  set.seed(seed)
  pts <- st_sfc(lapply(seq_len(n_geom), function(i) st_point(c(i, i %% 5))))
  pts <- st_set_crs(pts, 4326)
  times <- seq(as.Date("2000-01-01"), length.out = n_time, by = "month")
  arr <- array(rnorm(n_geom * n_time, mean = 20, sd = 5), dim = c(n_geom, n_time))
  st_as_stars(list(values = arr),
              dimensions = st_dimensions(geometry = pts, time = times))
}

test_that("area_weights errors clearly for zero-area (point) geometries", {
  expect_error(area_weights(make_point_stars()), class = "tidyeof_zero_area")
})

test_that("patterns guides users to weight = FALSE for point geometries", {
  dat <- make_point_stars()
  expect_error(patterns(dat, k = 3), class = "tidyeof_zero_area")
  expect_no_error(patterns(dat, k = 3, weight = FALSE))
})
