prism <- system.file("testdata/prism_test.RDS", package = "tidyeof") %>%
  readRDS()

# North et al. (1982) sampling-error bars and the multiplet/degeneracy
# coloring are diagnostics on the true eigenspectrum, used to choose k before
# rotating. Varimax produces non-eigenvalues, so those overlays are dropped on
# rotated patterns; the unrotated scree plot keeps them.

has_geom <- function(p, geom) {
  any(vapply(p$layers, function(l) inherits(l$geom, geom), logical(1)))
}

test_that("unrotated screeplot draws North error bars", {
  pat <- patterns(prism, k = 4, rotate = FALSE)
  p <- screeplot(pat)
  expect_true(has_geom(p, "GeomLinerange"))
})

test_that("rotated screeplot omits North error bars and multiplet coloring", {
  pat <- patterns(prism, k = 4, rotate = TRUE)
  p <- screeplot(pat)
  expect_false(has_geom(p, "GeomLinerange"))
  expect_true(has_geom(p, "GeomPoint"))
})

test_that("rotated screeplot warns and skips rule_n cutoff", {
  pat <- patterns(prism, k = 4, rotate = TRUE)
  expect_warning(screeplot(pat, rule_n = TRUE), class = "tidyeof_rule_n_rotated")
})
