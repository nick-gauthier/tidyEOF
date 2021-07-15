library(ecmwfr)

# set a key to the keychain
wf_set_key(user = "ngauthier91@gmail.com",
           key = "ebde449d7da12df582f01a20f233353b",
           service = "webapi")

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

dates <- paste0(sort(c(paste0(rep(1981:2009, each = 3) ,paste0(c('10','11', '12'), '01')),
                       paste0(rep(1982:2010, each = 3) ,paste0(c('01','02', '03'), '01')))), collapse = '/')
geop_request <- list(
  area    = "90/0/-90/360",
  class   = "ep",
  dataset = "cera20c",
  date    = dates,
  expver  = "1",
  grid    = "1.0/1.0",
  levelist= "500",
  levtype = "pl",
  number  = "0/1/2/3/4/5/6/7/8/9",
  param   = "129.128",
  stream  = "edmo",
  type    = 'an',
  format  = "netcdf",
  target  = "CERA-20C_geop.nc"
)

sst_request <- list(
  area    = "90/0/-90/360",
  class   = "ep",
  dataset = "cera20c",
  date    = dates,
  expver  = "1",
  grid    = "1.0/1.0",
  levtype = "sfc",
  number  = "0/1/2/3/4/5/6/7/8/9",
  param   = "34.128",
  stream  = "edmo",
  type    = 'an',
  format  = "netcdf",
  target  = "CERA-20C_sst.nc"
)

wf_request(request = snow_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data-raw")
wf_request(request = geop_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data-raw")
wf_request(request = sst_request, user = "ngauthier91@gmail.com",
           transfer = TRUE, path = "data-raw")
