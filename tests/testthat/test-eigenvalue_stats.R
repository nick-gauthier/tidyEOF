# Eigenvalue statistics: percent variance, North et al. (1982) error bars,
# and Rule N significance. These must give the same answers whether the full
# eigenvalue spectrum is available (base prcomp) or only the leading k modes
# (IRLBA), since prcomp_irlba reports the true total variance via $totalvar.

# Synthetic field with variance spread across modes: 4 planted signal modes
# (score sds 8, 5, 3, 2) on top of white noise, so percent-variance errors
# from a truncated-spectrum denominator are large and unmistakable.
make_spread_stars <- function(nt = 30, nx = 12, ny = 12, seed = 42) {
  set.seed(seed)
  ns <- nx * ny
  basis <- qr.Q(qr(matrix(rnorm(ns * 4), ns, 4)))
  scores <- matrix(rnorm(nt * 4), nt, 4) %*% diag(c(8, 5, 3, 2))
  mat <- scores %*% t(basis) + matrix(rnorm(nt * ns, sd = 0.5), nt, ns)
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

test_that("IRLBA eigenvalue table matches full-spectrum PCA", {
  skip_if_not_installed("irlba")
  dat <- make_spread_stars()

  pat_full <- patterns(dat, k = 4)
  suppressMessages(pat_irlba <- patterns(dat, k = 4, irlba_threshold = 1))

  for (col in c("percent", "cumulative", "low", "hi")) {
    expect_equal(pat_irlba$eigenvalues[[col]][1:4],
                 pat_full$eigenvalues[[col]][1:4],
                 tolerance = 1e-4)
  }
  # percent of the retained modes must reflect total variance, not sum to 100
  expect_lt(sum(pat_irlba$eigenvalues$percent[1:4]), 99)
})

test_that("North error bars use the number of time steps", {
  # more time steps (30) than grid cells (16), so any min(n, p)-based
  # shortcut for N gives the wrong answer
  dat <- make_spread_stars(nt = 30, nx = 4, ny = 4)
  pat <- patterns(dat, k = 3)
  expect_equal(unique(pat$eigenvalues$error), sqrt(2 / 30))

  skip_if_not_installed("irlba")
  dat2 <- make_spread_stars()
  suppressMessages(pat_irlba <- patterns(dat2, k = 4, irlba_threshold = 1))
  expect_equal(unique(pat_irlba$eigenvalues$error), sqrt(2 / 30))
})

test_that("eigen_test accepts total_var for truncated spectra", {
  set.seed(99)
  X <- scale(matrix(rnorm(40 * 100), 40, 100), scale = FALSE)
  lambdas <- prcomp(X, center = FALSE)$sdev^2

  for (k in 1:3) {
    expect_equal(
      eigen_test(lambdas[1:4], k = k, M = 100, n = 40, total_var = sum(lambdas)),
      eigen_test(lambdas, k = k, M = 100, n = 40)
    )
  }
})

test_that("rule_n_cutoff works for IRLBA-truncated eigenvalue tables", {
  skip_if_not_installed("irlba")
  dat <- make_spread_stars()

  pat_full <- patterns(dat, k = 4)
  suppressMessages(pat_irlba <- patterns(dat, k = 4, irlba_threshold = 1))

  expect_no_error(cutoff_irlba <- rule_n_cutoff(pat_irlba))
  cutoff_full <- rule_n_cutoff(pat_full)
  expect_equal(cutoff_irlba, min(cutoff_full, 4L))
})

test_that("print.patterns reports percent of total variance for IRLBA patterns", {
  skip_if_not_installed("irlba")
  dat <- make_spread_stars()

  pat_full <- patterns(dat, k = 4)
  suppressMessages(pat_irlba <- patterns(dat, k = 4, irlba_threshold = 1))

  out <- paste(cli::cli_fmt(print(pat_irlba)), collapse = " ")
  expected <- round(pat_full$eigenvalues$percent[1], 1)
  expect_match(out, as.character(expected), fixed = TRUE)
})
