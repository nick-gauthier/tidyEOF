
# devtools::install_github('khufkens/ecmwfr')
library(ecmwfr)

# set a key to the keychain
wf_set_key(user = "ngauthier91@gmail.com",
           key = "d0b10ca0286891178d1f4f672e2f96eb",
           service = "webapi")

# get the dates you want
dates <- paste0(paste0(rep(1901:2010, each = 12), paste0(c(paste0(0, 1:9), 10:12),'01')), collapse = '/')

prec_request <- list(
  area    = "70/-130/30/-60",
  class   = "ep",
  dataset = "cera20c",
  date    = dates,
  expver  = "1",
  grid    = "2.0/2.0",
  levtype = "sfc",
  number  = "0/1/2/3/4/5/6/7/8/9",
  param   = "228.128",
  stream  = "edmo",
  type    = "fc",
  format  = "netcdf",
  target  = "CERA-20C_precipitation.nc"
)

wf_request(request = prec_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data")

snow_request <- list(
  area    = "70/-130/30/-60",
  class   = "ep",
  dataset = "cera20c",
  date    = paste0(paste0(1901:2010,'03','01'), collapse = '/'),
  expver  = "1",
  grid    = "1.0/1.0",
  levtype = "sfc",
  number  = "0/1/2/3/4/5/6/7/8/9",
  param   = "33.128/141.128",
  stream  = "edmo",
  format  = "netcdf",
  target  = "CERA-20C_snow.nc"
)

wf_request(request = snow_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data")


era_dates <- paste0(paste0(1982:2017,'03','01'), collapse = '/')

era_snow_request <- list(
  area    = "70/-130/30/-60",
  class   = "ei",
  dataset = "interim",
  date    = era_dates,
  expver  = "1",
  grid    = "1.0/1.0",
  levtype = "sfc",
  param   = "33.128/141.128",
  stream  = "moda",
  type    = "an",
  format  = "netcdf",
  target  = "era-interim_snow.nc"
)

wf_request(request = era_snow_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data")
