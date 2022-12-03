prism <- system.file('testdata/prism_test.RDS', package = 'tidyEOF') %>%
  readRDS()
prism_coarse <- prism %>% st_warp(cellsize = 0.5, method = 'average', use_gdal = TRUE, no_data_value = -99999) %>%
  setNames(names(prism)) %>%
  mutate(across(everything(), ~units::set_units(.x, units(prism[[1]]), mode = 'standard'))) %>%
  st_set_dimensions('band', values = st_get_dimension_values(prism, 'time'), names = 'time')

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
