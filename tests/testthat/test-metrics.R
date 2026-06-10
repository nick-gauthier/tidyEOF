# Spatial correlation must measure agreement of the anomaly patterns, not the
# shared climatology. Two fields that sit on top of the same large, static
# spatial climatology look ~identical to a naive across-space correlation even
# when their anomalies are unrelated -- so the metric has to remove each cell's
# temporal mean first.

make_field <- function(clim_vec, anom, n_time) {
  clim <- matrix(rep(clim_vec, each = n_time), n_time, length(clim_vec))
  clim + anom
}

test_that("spatial correlation ignores shared climatology (unrelated anomalies -> ~0)", {
  set.seed(1)
  n_time <- 12; n_space <- 80
  clim_vec <- rnorm(n_space, sd = 20)          # large static climatology
  obs  <- make_field(clim_vec, matrix(rnorm(n_time * n_space), n_time), n_time)
  pred <- make_field(clim_vec, matrix(rnorm(n_time * n_space), n_time), n_time)

  expect_lt(abs(calc_cor_spatial(pred, obs)), 0.2)
})

test_that("spatial correlation is high when anomaly patterns agree", {
  set.seed(2)
  n_time <- 12; n_space <- 80
  clim_vec <- rnorm(n_space, sd = 20)
  anom <- matrix(rnorm(n_time * n_space), n_time)
  obs  <- make_field(clim_vec, anom, n_time)
  pred <- make_field(clim_vec, anom * 0.9 + matrix(rnorm(n_time * n_space, sd = 0.1), n_time), n_time)

  expect_gt(calc_cor_spatial(pred, obs), 0.8)
})

test_that("temporal correlation already ignores per-cell climatology", {
  # cor over time per cell is invariant to additive constants, so adding a
  # climatology must not change it -- a guard against regressions
  set.seed(3)
  n_time <- 20; n_space <- 30
  anom_o <- matrix(rnorm(n_time * n_space), n_time)
  anom_p <- anom_o * 0.8 + matrix(rnorm(n_time * n_space, sd = 0.3), n_time)
  clim_vec <- rnorm(n_space, sd = 15)

  bare <- calc_cor_temporal(anom_p, anom_o)
  with_clim <- calc_cor_temporal(make_field(clim_vec, anom_p, n_time),
                                 make_field(clim_vec, anom_o, n_time))
  expect_equal(bare, with_clim, tolerance = 1e-10)
})
