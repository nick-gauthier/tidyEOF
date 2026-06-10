# tune_eof() must measure genuine out-of-sample skill so that adding noise
# EOFs eventually hurts. It does this with a speckled holdout (Bro et al.
# 2008): in each held-out fold a random scatter of grid cells is hidden, the
# mode amplitudes are estimated from the visible cells, and the hidden cells
# are predicted. RMSE on hidden cells then has a true minimum near the real
# rank, instead of decreasing monotonically with k.

# Planted-rank field: exactly `ktrue` orthogonal spatial modes plus white
# noise. A correct tuner should select k at or near ktrue, not max(k).
make_planted_stars <- function(nt = 60, nx = 12, ny = 12, ktrue = 2,
                               sds = c(10, 6), noise = 2, seed = 11) {
  set.seed(seed)
  ns <- nx * ny
  basis <- qr.Q(qr(matrix(rnorm(ns * ktrue), ns, ktrue)))
  scores <- matrix(rnorm(nt * ktrue), nt, ktrue) %*% diag(sds, ktrue)
  mat <- scores %*% t(basis) + matrix(rnorm(nt * ns, sd = noise), nt, ns)
  arr <- array(t(mat), dim = c(nx, ny, nt))
  stars::st_as_stars(arr) %>%
    stars::st_set_dimensions(1, values = seq_len(nx), names = "x") %>%
    stars::st_set_dimensions(2, values = seq_len(ny), names = "y") %>%
    stars::st_set_dimensions(3,
      values = seq(as.Date("2000-01-01"), by = "year", length.out = nt),
      names = "time"
    ) %>%
    sf::st_set_crs("EPSG:32633")
}

test_that("tune_eof selects the planted rank, not the maximum k", {
  dat <- make_planted_stars(ktrue = 2, sds = c(10, 6))
  res <- tune_eof(dat, k = 1:8, kfolds = 4, weight = FALSE)
  best <- summarize_eof_cv(res, metric = "rmse")
  best_k <- attr(best, "best_k")

  expect_false(best_k == max(res$k))
  expect_lte(abs(best_k - 2), 1)
})

test_that("held-out RMSE rises again once k exceeds the planted rank", {
  dat <- make_planted_stars(ktrue = 2, sds = c(10, 6))
  res <- tune_eof(dat, k = 1:8, kfolds = 4, weight = FALSE)
  by_k <- res %>%
    dplyr::group_by(k) %>%
    dplyr::summarize(rmse = mean(rmse), .groups = "drop")

  rmse_at_truth <- by_k$rmse[by_k$k == 2]
  rmse_at_max <- by_k$rmse[by_k$k == 8]
  expect_lt(rmse_at_truth, rmse_at_max)
})

test_that("tune_eof is deterministic across repeated calls", {
  dat <- make_planted_stars()
  r1 <- tune_eof(dat, k = 1:4, kfolds = 3, weight = FALSE)
  r2 <- tune_eof(dat, k = 1:4, kfolds = 3, weight = FALSE)
  expect_equal(r1, r2)
})

test_that("hidden-cell masks are shared across k for a fair comparison", {
  # If masks differed by k, the k=2 and k=3 RMSEs would be scored on
  # different cells. With shared masks and a rank-2 truth, the residual at
  # k = 3 (one redundant mode) stays very close to k = 2.
  dat <- make_planted_stars(ktrue = 2, sds = c(12, 7), noise = 1.5)
  res <- tune_eof(dat, k = 2:3, kfolds = 4, weight = FALSE)
  by_k <- tapply(res$rmse, res$k, mean)
  expect_lt(abs(by_k[["3"]] - by_k[["2"]]) / by_k[["2"]], 0.1)
})
