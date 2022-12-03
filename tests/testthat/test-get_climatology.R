prism <- system.file('testdata/prism_test.RDS', package = 'tidyEOF') %>%
  readRDS()
clim <- get_climatology(prism)
clim_nounits <- get_climatology(units::drop_units(prism)) # does drop units also drop dimnmes? it doesn't on just prism but it does after get climatology?
clim_mon <- get_climatology(prism, monthly = TRUE)
clim_mon_nounits <- get_climatology(units::drop_units(prism), monthly = TRUE)

test_that('get_climatology() returns correct format', {
## add tests for multiple attributes?

  expect_true(inherits(prism, 'stars'))
  expect_true(inherits(clim, 'stars'))
  expect_true(inherits(clim_nounits, 'stars'))
  expect_true(inherits(clim_mon, 'stars'))
  expect_true(inherits(clim_mon_nounits, 'stars'))

  expect_equal(names(prism), 'tmean')
  expect_equal(names(clim), 'tmean')
  expect_equal(names(clim_nounits), 'tmean')
  expect_equal(names(clim_mon), 'tmean')
  expect_equal(names(clim_mon_nounits), 'tmean')

  expect_true(inherits(prism[[1]], 'units'))
  expect_true(inherits(clim[[1]], 'units'))
  expect_true(inherits(clim_mon[[1]], 'units'))
  expect_false(inherits(clim_nounits[[1]], 'units'))
  expect_false(inherits(clim_mon_nounits[[1]], 'units'))

  expect_equal(length(dim(prism)), 3)
  expect_equal(length(dim(clim)), 3)
  expect_equal(length(dim(clim_nounits)), 3)
  expect_equal(length(dim(clim_mon)), 4)
  expect_equal(length(dim(clim_mon_nounits)), 4)

  expect_equal(st_get_dimension_values(clim_mon, 'month'), month.name)
  expect_equal(st_get_dimension_values(clim_mon_nounits, 'month'), month.name)
})

test_that("get_climatology() works", {
  # random test that february averages are calculate correctly, could do more here
  expect_equal(st_apply(prism[,,,c(2, 14, 26)], 1:2, mean, rename = FALSE), abind::adrop(clim_mon_nounits[,,,2,1]))
  expect_equal(st_apply(prism[,,,c(2, 14, 26)], 1:2, sd, rename = FALSE), abind::adrop(clim_mon_nounits[,,,2,2]))

  })
