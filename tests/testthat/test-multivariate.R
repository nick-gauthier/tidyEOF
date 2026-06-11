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
})
