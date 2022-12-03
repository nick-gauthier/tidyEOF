prism <- system.file('testdata/prism_test.RDS', package = 'tidyEOF') %>%
  readRDS()

test_that("get_pcs() runs", {
  expect_no_error(get_pcs(prism))
})
