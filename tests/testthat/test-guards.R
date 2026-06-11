# Input guards and alignment fixes:
# - rotation requires k > 1 (previously silently skipped, stored rotate = TRUE)
# - multi-attribute stars inputs abort (previously silently used first attribute)
# - [.patterns only allows leading contiguous subsets (eigenvalue stats and
#   projection are undefined for non-contiguous subsets)
# - amplitude extraction and teleconnection maps align by time value, not
#   storage order

library(testthat)
library(stars)

prism <- readRDS(system.file("testdata/prism_test.RDS", package = "tidyeof"))

test_that("patterns aborts for rotate = TRUE with k = 1", {
  expect_error(
    patterns(prism, k = 1, rotate = TRUE),
    class = "tidyeof_invalid_option"
  )
})

test_that("patterns aborts for multi-attribute input", {
  two <- c(setNames(prism, "a"), setNames(prism, "b"))
  expect_error(patterns(two, k = 2), class = "tidyeof_multiple_attributes")
})

test_that("common_patterns aborts for multi-attribute input", {
  two <- c(setNames(prism, "a"), setNames(prism, "b"))
  expect_error(
    common_patterns(list(src = two), k = 2),
    class = "tidyeof_multiple_attributes"
  )
})

test_that("project_patterns aborts for multi-attribute newdata", {
  pat <- patterns(prism, k = 2)
  two <- c(setNames(prism, "a"), setNames(prism, "b"))
  expect_error(project_patterns(pat, two), class = "tidyeof_multiple_attributes")
})

test_that("[.patterns allows leading contiguous subsets only", {
  pat <- patterns(prism, k = 3)

  truncated <- pat[1:2]
  expect_s3_class(truncated, "patterns")
  expect_equal(truncated$k, 2)

  # Character subsetting resolves against the actual amplitude names
  expect_equal(pat[c("PC1", "PC2")]$k, 2)

  expect_error(pat[c(1, 3)], class = "tidyeof_invalid_subset")
  expect_error(pat[2:3], class = "tidyeof_invalid_subset")
  expect_error(pat[2], class = "tidyeof_invalid_subset")
  expect_error(pat[1:4], class = "tidyeof_invalid_subset")
})

test_that("extract_amplitudes_matrix with a times filter returns time-sorted rows", {
  amps <- tibble::tibble(
    time = as.Date(c("2002-01-01", "2000-01-01", "2001-01-01")),
    PC1 = c(3, 1, 2)
  )
  m <- extract_amplitudes_matrix(amps, times = amps$time)
  expect_equal(as.vector(m), c(1, 2, 3))
})

test_that("get_correlation aligns by time value when the time dimension is unsorted", {
  times_shuffled <- as.Date(c(
    "2000-05-01", "2000-01-01", "2000-03-01", "2000-02-01", "2000-04-01"
  ))
  vals <- as.numeric(times_shuffled)
  arr <- array(rep(vals, each = 4), c(2, 2, 5))
  names(dim(arr)) <- c("x", "y", "time")
  dat <- st_as_stars(arr) |>
    setNames("v") |>
    st_set_dimensions(1, values = 1:2, name = "x") |>
    st_set_dimensions(2, values = 1:2, name = "y") |>
    st_set_dimensions(3, values = times_shuffled, name = "time")

  # Amplitudes equal the pixel series, so a correctly aligned correlation is 1
  amps <- tibble::tibble(
    time = sort(times_shuffled),
    PC1 = as.numeric(sort(times_shuffled))
  )

  cors <- get_correlation(dat, patterns = NULL, amplitudes = amps)
  expect_true(all(abs(cors[["PC1"]] - 1) < 1e-12))
})
