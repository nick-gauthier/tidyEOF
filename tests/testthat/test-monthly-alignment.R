# Tests for calendar-indexed monthly climatology.
#
# The month dimension of a monthly climatology is always 1:12 in calendar
# order, so applying a climatology to data is direct indexing by month number.
# These tests pin down the behaviors that positional month matching broke:
# new data starting in a different calendar month than the training data,
# single-year inputs, and partial (non-whole-year) records.

library(testthat)
library(stars)

# Data with a pure seasonal cycle: value = calendar month number, no noise.
# Anomalies w.r.t. a correctly aligned monthly climatology are exactly zero.
mk_raster <- function(start, n_months) {
  times <- seq(as.Date(start), by = "month", length.out = n_months)
  vals <- as.integer(format(times, "%m"))
  arr <- array(rep(vals, each = 9), c(3, 3, n_months))
  names(dim(arr)) <- c("x", "y", "time")
  st_as_stars(arr) |>
    setNames("t2m") |>
    st_set_dimensions(1, values = 1:3, name = "x") |>
    st_set_dimensions(2, values = 1:3, name = "y") |>
    st_set_dimensions(3, values = times, name = "time")
}

mk_geom <- function(start, n_months) {
  times <- seq(as.Date(start), by = "month", length.out = n_months)
  vals <- as.integer(format(times, "%m"))
  geom <- sf::st_sfc(lapply(1:5, function(i)
    sf::st_polygon(list(rbind(c(i, 0), c(i + 1, 0), c(i + 1, 1), c(i, 1), c(i, 0))))))
  m <- matrix(rep(vals, each = 5), nrow = 5, ncol = n_months)
  st_as_stars(list(t2m = m),
              dimensions = st_dimensions(geometry = geom, time = times))
}

test_that("monthly climatology has a calendar month dimension (1:12) regardless of start month", {
  clim <- get_climatology(mk_raster("2000-07-01", 24), monthly = TRUE)
  expect_equal(as.integer(st_get_dimension_values(clim$mean, "month")), 1:12)
  for (mm in 1:12) {
    expect_true(all(clim$mean$t2m[, , mm] == mm))
  }
})

test_that("monthly anomalies are correct when newdata starts in a different month (raster)", {
  clim <- get_climatology(mk_raster("2000-01-01", 24), monthly = TRUE)
  anom <- get_anomalies(mk_raster("2003-07-01", 24), clim = clim, monthly = TRUE)
  expect_true(all(abs(anom$t2m) < 1e-12))
})

test_that("monthly anomalies are correct when newdata starts in a different month (geometry)", {
  clim <- get_climatology(mk_geom("2000-01-01", 24), monthly = TRUE)
  anom <- get_anomalies(mk_geom("2003-07-01", 24), clim = clim, monthly = TRUE)
  expect_true(all(abs(anom$t2m) < 1e-12))
})

test_that("a single year of new data works with a monthly climatology", {
  clim <- get_climatology(mk_raster("2000-01-01", 24), monthly = TRUE)
  anom <- get_anomalies(mk_raster("2003-01-01", 12), clim = clim, monthly = TRUE)
  expect_true(all(abs(anom$t2m) < 1e-12))
})

test_that("partial years compute a climatology and round-trip cleanly", {
  dat <- mk_raster("2000-01-01", 30)  # 2.5 years
  set.seed(42)
  dat[["t2m"]] <- dat[["t2m"]] + rnorm(length(dat[["t2m"]]), sd = 0.1)

  clim <- suppressMessages(get_climatology(dat, monthly = TRUE))

  for (scale in c(FALSE, TRUE)) {
    anom <- get_anomalies(dat, clim = clim, scale = scale, monthly = TRUE)
    restored <- restore_climatology(anom, clim, scale = scale, monthly = TRUE)
    expect_equal(dat$t2m, restored$t2m, tolerance = 1e-8)
  }
})

test_that("anomalies abort clearly when newdata has a month missing from the climatology", {
  # Climatology from January-June only
  clim <- suppressMessages(get_climatology(mk_raster("2000-01-01", 6), monthly = TRUE))
  expect_error(
    get_anomalies(mk_raster("2003-07-01", 3), clim = clim, monthly = TRUE),
    class = "tidyeof_missing_month"
  )
})

test_that("monthly anomalies preserve units for unscaled data", {
  dat <- mk_raster("2000-01-01", 24)
  dat[["t2m"]] <- units::set_units(dat[["t2m"]], "degC")
  anom <- get_anomalies(dat, monthly = TRUE)
  expect_true(inherits(anom$t2m, "units"))
})

test_that("monthly CV works with folds that are not whole years", {
  times <- seq(as.Date("2000-01-01"), by = "month", length.out = 36)
  set.seed(1)
  arr <- array(rnorm(3 * 3 * 36), c(3, 3, 36))
  names(dim(arr)) <- c("x", "y", "time")
  dat <- st_as_stars(arr) |>
    setNames("t2m") |>
    st_set_dimensions(1, values = 1:3, name = "x") |>
    st_set_dimensions(2, values = 1:3, name = "y") |>
    st_set_dimensions(3, values = times, name = "time")

  # 36 months / 5 folds -> 8,7,7,7,7-month folds; training sets of 28-29
  # months. Under the complete-years requirement this aborted mid-run.
  results <- suppressMessages(
    tune_eof(dat, k = 1:2, kfolds = 5, monthly = TRUE, n_reps = 2)
  )
  expect_s3_class(results, "tbl_df")
  expect_equal(nrow(results), 10)  # 2 k values x 5 folds
  expect_true(all(is.finite(results$rmse)))
})
