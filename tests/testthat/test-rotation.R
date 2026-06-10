# Varimax rotation follows the standard REOF convention (Hannachi et al.
# 2007): the rotation criterion is computed on sqrt(eigenvalue)-scaled
# loadings, and the stored spatial patterns are those rotated loadings
# (unit-norm, with the variance carried by the amplitudes). The displayed
# patterns must therefore be exactly the varimax solution.

make_rotation_stars <- function(nt = 40, nx = 10, ny = 10, seed = 1) {
  set.seed(seed)
  ns <- nx * ny
  basis <- qr.Q(qr(matrix(rnorm(ns * 4), ns, 4)))
  scores <- matrix(rnorm(nt * 4), nt, 4) %*% diag(c(9, 6, 4, 2))
  mat <- scores %*% t(basis) + matrix(rnorm(nt * ns, sd = 0.3), nt, ns)
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

test_that("rotated EOFs are the varimax solution (standard REOF convention)", {
  dat <- make_rotation_stars()
  k <- 4
  # weight = FALSE so the stored patterns live in plain (unweighted) space,
  # directly comparable to a reference computed from the anomaly matrix
  pat <- patterns(dat, k = k, rotate = TRUE, weight = FALSE)

  # Independent reference: varimax on Kaiser-scaled loadings
  anom <- get_anomalies(dat, get_climatology(dat))
  X <- flatten_time_space(anom[1])$matrix
  pca <- prcomp(X, center = FALSE)
  L <- pca$rotation[, 1:k] %*% diag(pca$sdev[1:k])
  LR <- unclass(varimax(L)$loadings)
  LR <- LR[, order(colSums(LR^2), decreasing = TRUE)]
  ref <- sweep(LR, 2, sqrt(colSums(LR^2)), `/`)  # unit-norm patterns

  eof_mat <- matrix(pat$eofs[[1]], nrow = prod(dim(dat)[1:2]), ncol = k)
  for (j in seq_len(k)) {
    s <- sign(sum(eof_mat[, j] * ref[, j]))  # sign convention is arbitrary
    expect_equal(eof_mat[, j], s * ref[, j], tolerance = 1e-6)
  }
})

test_that("rotated amplitudes are uncorrelated with sd = sqrt(rotated eigenvalue)", {
  pat <- patterns(make_rotation_stars(), k = 4, rotate = TRUE, weight = FALSE)
  amps <- as.matrix(pat$amplitudes[, -1])

  expect_equal(unname(apply(amps, 2, sd)),
               pat$eigenvalues$std.dev[1:4], tolerance = 1e-6)

  cors <- cor(amps)
  expect_equal(unname(cors), diag(4), tolerance = 1e-6)
})

test_that("project_patterns reproduces rotated amplitudes", {
  dat <- make_rotation_stars()
  for (w in c(FALSE, TRUE)) {
    pat <- patterns(dat, k = 4, rotate = TRUE, weight = w)
    reproj <- project_patterns(pat, dat)
    expect_equal(as.matrix(reproj[, -1]), as.matrix(pat$amplitudes[, -1]),
                 tolerance = 1e-8)
  }
})

# Sign-orientation must never multiply a mode by zero. A perfectly balanced
# pattern (loadings summing to exactly zero) previously got sign(0) = 0, which
# zeroed the EOF, its amplitudes, and the projection matrix.
test_that("balanced EOFs (zero-sum loadings) get sign +1, not 0", {
  arr <- array(0, dim = c(2, 2, 2))
  arr[, , 1] <- 1                          # PC1: sums positive
  arr[, , 2] <- matrix(c(1, -1, 1, -1), 2) # PC2: sums to exactly zero
  eofs <- stars::st_as_stars(arr) %>%
    stars::st_set_dimensions(1, names = "x") %>%
    stars::st_set_dimensions(2, names = "y") %>%
    stars::st_set_dimensions(3, values = c("PC1", "PC2"), names = "PC")

  signs <- compute_eof_signs(eofs)
  expect_setequal(unique(signs), c(1))   # both modes oriented +1
  expect_false(any(signs == 0))
})
