# Calculate climatological mean and standard deviation for spatial data

Computes climatological statistics (mean and standard deviation) for a
spatial field, either annually or monthly. Preserves spatial dimensions
and units from the input data.

## Usage

``` r
get_climatology(dat, monthly = FALSE)
```

## Arguments

- dat:

  A stars object containing a spatial field with dimensions (x, y, time)

- monthly:

  Logical. If TRUE, computes monthly climatology. If FALSE (default),
  computes statistics over the entire period.

## Value

A list with two stars objects:

- mean:

  Climatological mean with original spatial dimensions and units

- sd:

  Climatological standard deviation with same structure

## Examples

``` r
# Create sample data
library(stars)
#> Loading required package: abind
#> Loading required package: sf
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE
times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
x <- seq(0, 1, length.out = 10)
y <- seq(0, 1, length.out = 10)
dat <- stars::st_as_stars(array(rnorm(10*10*36), c(10, 10, 36))) %>%
  st_set_dimensions(1, values = x, name = "x") %>%
  st_set_dimensions(2, values = y, name = "y") %>%
  st_set_dimensions(3, values = times, name = "time")
#> Error in stars::st_as_stars(array(rnorm(10 * 10 * 36), c(10, 10, 36))) %>%     st_set_dimensions(1, values = x, name = "x") %>% st_set_dimensions(2,     values = y, name = "y") %>% st_set_dimensions(3, values = times,     name = "time"): could not find function "%>%"

# Calculate annual climatology
clim <- get_climatology(dat)
#> Error: object 'dat' not found

# Calculate monthly climatology
monthly_clim <- get_climatology(dat, monthly = TRUE)
#> Error: object 'dat' not found
```
