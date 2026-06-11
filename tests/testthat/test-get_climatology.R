# Tests for get_climatology function
library(testthat)
library(stars)
library(units)
library(lubridate)

test_that("get_climatology handles basic annual calculations correctly", {
  # Create test data with known values
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- seq(0, 1, length.out = 5)
  y <- seq(0, 1, length.out = 5)

  set.seed(123)
  test_array <- array(rnorm(5*5*36), c(5, 5, 36))

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') |>
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate climatology
  result <- get_climatology(test_data, monthly = FALSE)

  # Basic structure tests — now returns a list
  expect_type(result, "list")
  expect_named(result, c("mean", "sd"))
  expect_s3_class(result$mean, "stars")
  expect_s3_class(result$sd, "stars")

  # Calculate expected values
  manual_mean <- apply(test_array, c(1,2), mean)
  names(dim(manual_mean)) <- c('x', 'y')
  manual_sd <- apply(test_array, c(1,2), sd)
  names(dim(manual_sd)) <- c('x', 'y')

  # Compare values
  expect_equal(result$mean$t2m, manual_mean, tolerance = 1e-7)
  expect_equal(result$sd$t2m, manual_sd, tolerance = 1e-7)
})

test_that("get_climatology handles monthly calculations correctly", {
  # Create test data with known monthly pattern
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- seq(0, 1, length.out = 3)
  y <- seq(0, 1, length.out = 3)

  # Create data where each month has a constant value
  month_values <- 1:12
  test_array <- array(0, c(3, 3, 36))
  for(i in 1:36) {
    test_array[,,i] <- month_values[(i-1) %% 12 + 1]
  }

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') |>
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate monthly climatology
  result <- get_climatology(test_data, monthly = TRUE)

  # Test structure
  expect_type(result, "list")
  expect_named(result, c("mean", "sd"))
  expect_s3_class(result$mean, "stars")

  # Test values — check mean and sd for each month
  for(i in 1:12) {
    # All grid cells should have value i for mean
    ref <- matrix(c(i, i, i), ncol = 3)
    names(dim(ref)) <- c('x', 'y')
    mean_vals <- result$mean$t2m[,,i]
    expect_equal(unique(mean_vals), ref)
    # All grid cells should have sd of 0 (constant values)
    ref <- matrix(c(0, 0, 0), ncol = 3)
    names(dim(ref)) <- c('x', 'y')
    sd_vals <- result$sd$t2m[,,i]
    expect_equal(unique(sd_vals), ref)
  }
})

test_that("get_climatology preserves units correctly", {
  # Create test data with units
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- seq(0, 1, length.out = 3)
  y <- seq(0, 1, length.out = 3)

  test_array <- array(rnorm(3*3*36, mean = 20, sd = 5), c(3, 3, 36))

  test_data <- st_as_stars(test_array) %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time") %>%
    mutate(A1 = set_units(A1, "degC"))

  # Calculate climatologies
  annual_result <- get_climatology(test_data, monthly = FALSE)
  monthly_result <- get_climatology(test_data, monthly = TRUE)

  # Check units are preserved on both mean and sd
  expect_true(inherits(annual_result$mean[[1]], "units"))
  expect_true(inherits(annual_result$sd[[1]], "units"))
  expect_true(inherits(monthly_result$mean[[1]], "units"))

  # Compare unit strings
  expect_equal(as.character(units(annual_result$mean[[1]])), "°C")
  expect_equal(as.character(units(monthly_result$mean[[1]])), "°C")
})

test_that("get_climatology handles invalid inputs appropriately", {
  # Test non-stars input
  expect_error(get_climatology(matrix(1:12, 3, 4)),
               "Input must be a stars object")
})

test_that("get_anomalies handles annual calculations correctly", {
  # Create test data
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- seq(0, 1, length.out = 5)
  y <- seq(0, 1, length.out = 5)

  set.seed(123)
  test_array <- array(rnorm(5*5*36, mean = 10), c(5, 5, 36))

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') |>
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate anomalies
  result <- get_anomalies(test_data)

  # Tests
  expect_s3_class(result, "stars")
  expect_equal(dim(result), dim(test_data))

  # Check calculations
  manual_mean <- apply(test_array, c(1,2), mean)
  expect_equal(
    mean(result$t2m),
    0,
    tolerance = 1e-10
  )
})

test_that("get_anomalies handles monthly calculations correctly", {
  # Create test data with known monthly pattern
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- 1:3
  y <- 1:3

  # Create data where each month has a constant value
  month_pattern <- rep(1:12, 3)  # 3 years
  test_array <- array(rep(month_pattern, each = 9), c(3, 3, 36))

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') |>
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate monthly anomalies
  result <- get_anomalies(test_data, monthly = TRUE)

  # Each month's anomalies should be zero since values are constant per month
  expect_true(all(abs(result$t2m) < 1e-10))
})

test_that("get_anomalies monthly calculation works", {
  # Create test data
  times <- seq(ymd("2000-01-01"), ymd("2002-12-31"), by = "month")
  x <- 1:3
  y <- 1:3

  # Create data with known monthly pattern (value equals month number)
  test_array <- array(0, c(3, 3, 36))
  for(i in 1:36) {
    test_array[,,i] <- (i-1) %% 12 + 1
  }

  test_data <- st_as_stars(test_array) %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate anomalies
  result <- get_anomalies(test_data, monthly = TRUE)

  # Check dimensions
  expect_equal(dim(result), c(x = 3, y = 3, time = 36))
  expect_equal(st_get_dimension_values(result, "time"), times)

  # For each month, anomalies should be zero (each month's values equal its mean)
  for(m in 1:12) {
    month_vals <- as.vector(result[,,,which(month(times) == m)]$A1)
    expect_equal(month_vals, rep(0, length(month_vals)))
  }
})


test_that("get_anomalies and restore_climatology are inverses", {
  # Test with dates
  times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
  x <- 1:3
  y <- 1:3

  set.seed(123)
  test_array <- array(rnorm(3*3*36, mean = 20, sd = 5), c(3, 3, 36))
  names(dim(test_array)) <- c("x", "y", "time")
  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Test both monthly and non-monthly cases
  for (monthly in c(TRUE, FALSE)) {
    for (scale in c(TRUE, FALSE)) {
      # Get anomalies and restore
      anom <- get_anomalies(test_data, monthly = monthly, scale = scale)
      restored <- restore_climatology(anom, get_climatology(test_data, monthly = monthly),
                                      monthly = monthly, scale = scale)

      # Compare original and restored values
      expect_equal(test_data$t2m, restored$t2m, tolerance = 1e-10)
    }
  }
})

test_that("functions handle incomplete years gracefully", {
  # Create data with incomplete year
  times <- seq(as.Date("2000-01-01"), as.Date("2002-05-31"), by = "month")  # 29 months
  x <- 1:3
  y <- 1:3

  test_array <- array(1:261, c(3, 3, 29))  # 29 months of data

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Should message about unequal month sample sizes but still compute
  expect_message(clim <- get_climatology(test_data, monthly = TRUE),
                 class = "tidyeof_unbalanced_months")

  # Results should still have correct structure
  expect_type(clim, "list")
  expect_s3_class(clim$mean, "stars")

  # Anomalies no longer require complete years and round-trip exactly
  anom <- get_anomalies(test_data, clim = clim, monthly = TRUE)
  restored <- restore_climatology(anom, clim, monthly = TRUE)
  expect_equal(test_data$t2m, restored$t2m, tolerance = 1e-10)
})

test_that("functions handle non-January start dates correctly", {
  # Create test data starting in July
  times <- seq(ymd("2000-07-01"), ymd("2002-06-30"), by = "month")
  x <- 1:3
  y <- 1:3

  # Create data where each month has a unique value
  test_array <- array(rep(1:12, each = 9), c(3, 3, 24))  # 2 complete years
  names(dim(test_array)) <- c("x", "y", "time")
  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Calculate climatology
  clim <- get_climatology(test_data, monthly = TRUE)

  # Month dimension is always calendar-ordered 1:12, regardless of start month
  expect_equal(as.integer(st_get_dimension_values(clim$mean, "month")), 1:12)

  # Climatology of calendar month m matches the data values for that month
  for(m in 1:12) {
    month_data <- test_array[,,which(lubridate::month(times) == m)]
    expect_equal(unique(as.vector(clim$mean$t2m[,,m])), unique(as.vector(month_data)))
  }

  # Calculate anomalies (zero, since each month is constant across years)
  anom <- get_anomalies(test_data, monthly = TRUE)
  expect_true(all(abs(anom$t2m) < 1e-10))

  # Test restore_climatology
  restored <- restore_climatology(anom, clim, monthly = TRUE)
  expect_equal(test_data$t2m, restored$t2m, tolerance = 1e-10)
})

test_that("functions handle partial years correctly with non-standard start", {
  # Create test data with partial years starting in October
  times <- seq(ymd("2000-10-01"), ymd("2002-03-31"), by = "month")  # 18 months
  x <- 1:3
  y <- 1:3

  test_array <- array(1:162, c(3, 3, 18))

  test_data <- st_as_stars(test_array) %>%
    setNames('t2m') %>%
    st_set_dimensions(1, values = x, name = "x") %>%
    st_set_dimensions(2, values = y, name = "y") %>%
    st_set_dimensions(3, values = times, name = "time")

  # Should message about unequal month sample sizes but still compute
  expect_message(clim <- get_climatology(test_data, monthly = TRUE),
                 class = "tidyeof_unbalanced_months")

  # Month dimension is calendar-ordered; every month is present in this record
  expect_equal(as.integer(st_get_dimension_values(clim$mean, "month")), 1:12)
  expect_true(all(!is.na(clim$mean$t2m)))
})

test_that("cycle of operations preserves values for different start months", {
  # Test multiple different starting months
  start_months <- c(1, 4, 7, 10)  # Test quarterly starts

  for(start_month in start_months) {
    # Create 2 years of data starting in different months
    times <- seq(ymd(paste0("2000-", start_month, "-01")),
                 length.out = 24, by = "month")
    x <- 1:3
    y <- 1:3

    test_array <- array(rnorm(3*3*24, mean = 20, sd = 5), c(3, 3, 24))
    names(dim(test_array)) <- c("x", "y", "time")

    test_data <- st_as_stars(test_array) %>%
      setNames('t2m') %>%
      st_set_dimensions(1, values = x, name = "x") %>%
      st_set_dimensions(2, values = y, name = "y") %>%
      st_set_dimensions(3, values = times, name = "time")

    # Full cycle of operations
    clim <- get_climatology(test_data, monthly = TRUE)
    anom <- get_anomalies(test_data, monthly = TRUE)
    restored <- restore_climatology(anom, clim, monthly = TRUE)

    # Month dimension is calendar-ordered 1:12 regardless of start month
    expect_equal(as.integer(st_get_dimension_values(clim$mean, "month")), 1:12)

    # Check value preservation
    expect_equal(test_data$t2m, restored$t2m, tolerance = 1e-10)
  }
})
