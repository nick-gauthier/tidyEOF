## code to prepare `DATASET` dataset goes here

library(raster) # processing raster data
library(tidyverse) # data manipulation and visualization
library(stars)

# download state boundary files for cropping and plotting
state_names <- c('arizona', 'new mexico', 'colorado',
                 'california', 'utah', 'nevada',
                 'oregon', 'washington', 'idaho',
                 'wyoming', 'montana')

states_wus <- maps::map('state', regions = state_names,
                        fill = TRUE, plot = FALSE) %>%
  st_as_sf()

world <- maps::map('world',wrap=c(0,360),fill = TRUE, plot = FALSE) %>%
  st_as_sf()


## Geographic Data

#Define a study area to constrain all computations.
bbox <- extent(c(-125, -102, 31, 49))


preprocess <- function(x, var, bbox, flip = FALSE, regrid = FALSE, daily = FALSE) {
  maps <- brick(x, varname = var)

  indices <- getZ(maps) %>% # find the index for April 1st
    str_detect('-03-')%>%
    which()

  raster::subset(maps, indices) %>%
    {if(flip) rotate(.) else .} %>%
    crop(bbox) %>%
    {if(daily) mean(.) else .} %>% # if a file is one year of daily values
    {if(regrid) aggregate(., fact = 2) else .} # resample to lower res to speed up analysis
}

prism <- list.files('data-raw/UA-SWE', pattern = 'SWE_Depth', full.names = TRUE) %>%
  map(preprocess, var = 'SWE', bbox = bbox, regrid = TRUE, daily = TRUE) %>%
  brick() %>%
  mask(., all(near(., 0)), maskvalue = 1) %>% # mask pixels that never receive snow
  mask(states_wus) %>% # crop to state boundaries
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 1982:2017, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm))

## cera
cera <- map(1:10, ~ (brick('data-raw/CERA-20C_snow.nc', varname = 'sd', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  crop(bbox) %>%
  mask(., all(near(., 0)), maskvalue = 1) %>%
  mask(states_wus) %>%
  `*`(1000) %>% # convert to mm
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 1901:2010, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm))

# the CERA data are in mm SWE for the land fraction of the grid cell. Need to correct for this
land_frac <- read_ncdf('data-raw/_grib2netcdf-webmars-public-svc-blue-003-6fe5cac1a363ec1525f54343b6cc9fd8-tyqPrq.nc',
                       var = 'lsm') %>%
  adrop(4) %>%
  slice('number', 1) %>%
  .[st_bbox(cera)]

cera <- cera * land_frac

# cesm
cesm_h2osno <- preprocess('data-raw/b.e11.BLMTRC5CN.f19_g16.001.clm2.h0.H2OSNO.085001-184912.nc',
                          var = 'H2OSNO',
                          flip = TRUE,
                          bbox = bbox)

cesm_h2osno_ext <- preprocess('data-raw/b.e11.BLMTRC5CN.f19_g16.001.clm2.h0.H2OSNO.185001-200512.nc',
                          var = 'H2OSNO',
                          flip = TRUE,
                          bbox = bbox)

cesm <- c(cesm_h2osno, cesm_h2osno_ext) %>%
  brick() %>%
  raster::resample(as(cera, 'Raster')) %>%
  mask(mean(as(cera, 'Raster'))) %>%
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(SWE = units::set_units(if_else(SWE < 0, 0, SWE), mm)) # fix edge effects from resampling

# fix edge effects from resampling
#cesm[cesm < 0] <- 0

ccsm_lm <- preprocess('data-raw/snw_LImon_CCSM4_past1000_r1i1p1_085001-185012.nc',
                   var = 'snw',
                   flip = TRUE,
                   bbox = bbox)

ccsm_ext <- preprocess('data-raw/snw_LImon_CCSM4_historical_r1i2p1_185001-200512.nc',
                   var = 'snw',
                   flip = TRUE,
                   bbox = bbox)[[-1]] # the datasets overlap in year 1850

ccsm <- c(ccsm_lm, ccsm_ext) %>%
  brick() %>%
  raster::resample(as(cera, 'Raster')) %>%
  mask(mean(as(cera, 'Raster'))) %>%
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(SWE = units::set_units(if_else(SWE < 0, 0, SWE), mm)) # fix edge effects from resampling

# fix edge effects from resampling
#ccsm[ccsm < 0] <- 0

#plot(mean(prism));plot(mean(cera));plot(mean(cesm));plot(mean(ccsm))

###### climate data for teleconnection analysis

get_prism_time <- function(x) {
  st_get_dimension_values(x, 'band') %>%
  str_split('_') %>%
  map_chr(~.[[5]]) %>%
  parse_date(format = '%Y%m')
}

# precipitation
ppt <- list.files('data-raw/PRISM_ppt_stable_4kmM3_198101_201904_bil', full.names = TRUE, pattern = '.bil$') %>%
  map(~raster(.) %>% crop(bbox)) %>%
  brick %>%
  aggregate(fact = 2) %>%
  crop(states_wus) %>%
  st_as_stars() %>%
  setNames('ppt') %>%
  st_set_dimensions(., 'band', values = get_prism_time(.), names = 'time')

ppt_jfm <- ppt %>%
  filter(format(time, "%m") %in% c('01', '02', '03')) %>%
  aggregate(by = 'years', sum, na.rm = TRUE) %>%
  st_set_dimensions(., 'time', values = 1981:2019)

ppt_ond <- ppt %>%
  filter(format(time, "%m") %in% c('10', '11', '12')) %>%
  aggregate(by = 'years', sum, na.rm = TRUE) %>%
  st_set_dimensions(., 'time', values = (1981:2018) + 1) # + 1 for water year

# temperature
tmean <- list.files('data-raw/PRISM_tmean_stable_4kmM3_198101_201904_bil', full.names = TRUE, pattern = '.bil$') %>%
  map(~raster(.) %>% crop(bbox)) %>%
  brick %>%
  aggregate(fact = 2) %>%
  crop(states_wus) %>%
  st_as_stars() %>%
  setNames('tmean') %>%
  st_set_dimensions(., 'band', values = get_prism_time(.), names = 'time')

tmean_jfm <- tmean %>%
  filter(format(time, "%m") %in% c('01', '02', '03')) %>%
  aggregate(by = 'years', mean, na.rm = TRUE) %>%
  st_set_dimensions('time', values = 1981:2019)

tmean_ond <- tmean %>%
  filter(format(time, "%m") %in% c('10', '11', '12')) %>%
  aggregate(by = 'years', mean, na.rm = TRUE) %>%
  st_set_dimensions('time', values = (1981:2018) + 1) # + 1 for water year

# combine
prism_clim <- c(ppt_jfm[,-1] + ppt_ond,
                (tmean_jfm[,-1] + tmean_ond) / 2) %>%
  mutate(ppt = units::set_units(ppt, mm),
         tmean = units::set_units(tmean, '°C'))

# sst
sst_jfm <- read_ncdf('data-raw/sst.mnmean.nc', sub = 'sst') %>%
  st_crop(st_bbox(c(xmin = -1, xmax = 359, ymin = -75, ymax = 75), crs = 4326)) %>%
  filter(between(time, as.POSIXct('1981-09-01'), as.POSIXct('2017-03-01')),
         format(time, '%m') %in% c('01', '02', '03')) %>%
  aggregate(by = 'years', mean, na.rm = TRUE) %>%
  st_set_dimensions('time', values = 1982:2017)

sst_ond <- read_ncdf('data-raw/sst.mnmean.nc', sub = 'sst') %>%
  st_crop(st_bbox(c(xmin = -1, xmax = 359, ymin = -75, ymax = 75), crs = 4326)) %>%
  filter(between(time, as.POSIXct('1981-09-01'), as.POSIXct('2017-03-01')),
         format(time, '%m') %in% c('10','11','12')) %>%
  aggregate(by = 'years', mean, na.rm = TRUE) %>%
  st_set_dimensions('time', values = (1981:2016) + 1) # + 1 for water year

sst <- ((sst_jfm + sst_ond) / 2) %>%
  mutate(sst = units::set_units(sst, '°C'))

geop_jfm <- brick('data-raw/adaptor.mars.internal-1584737536.594777-6445-11-8d0fcc69-bd7d-40e2-a423-cf86c741ac79.nc')[[-(1:3)]] %>%
   stackApply(rep(1:36, each = 3), 'mean') %>%
  aggregate(fact = 8) %>%
  st_as_stars() %>%
  st_set_dimensions('band', values = 1982:2017, names = 'time')

geop_ond <- brick('data-raw/adaptor.mars.internal-1584737536.594777-6445-11-8d0fcc69-bd7d-40e2-a423-cf86c741ac79.nc')[[-(1:3)]] %>%
  stackApply(rep(1:36, each = 3), 'mean') %>%
  #   `-`(., mean(.)) %>% # is this necessary?
  aggregate(fact = 4, method = 'bilinear') %>%
  setNames(1982:2017) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = parse_number(layer)) %>%
  select(-layer)

#mountains <- read_sf('../data/ne_10m_geography_regions_polys.shp') %>%
#  filter(name %in% c('SIERRA NEVADA', 'CASCADE RANGE', 'ROCKY MOUNTAINS')) %>%
#  dplyr::select(name, geometry)

#mountains_mask <- mountains %>%
#  mutate(rasters = split(., 1:3) %>% map(~mask(prism[[1]], .))) %>%
#  st_drop_geometry() %>%
#  mutate(rasters = map(rasters, as.data.frame, na.rm = TRUE, xy = TRUE)) %>%
#  unnest(rasters) %>%
#  dplyr::select(-X1982)

#Turn the raster bricks into data frames and join.
# prism_dat <- prism %>%
#   as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
#   mutate(year = round(parse_number(layer))) %>%
#   rename(SWE = value) %>%
#   dplyr::select(-layer)
#
# cera_dat <- cera %>%
#   as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
#   mutate(year = round(parse_number(layer))) %>%
#   rename(SWE = value) %>%
#   dplyr::select(-layer)
#
# cesm_dat <- cesm %>%
#   as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
#   mutate(year = str_sub(layer, 2) %>% parse_number() %>% round) %>%
#   rename(SWE = value) %>%
#   select(-layer)
#
# ccsm_dat <- ccsm %>%
#   as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
#   mutate(year = str_sub(layer, 2) %>% parse_number() %>% round) %>%
#   rename(SWE = value) %>%
#   select(-layer)

#usethis::use_data(prism_dat, cera_dat, cesm_dat, ccsm_dat, tmean_dat, ppt_dat, geop_jfm, sst_dat, areas_prism, areas_cera, states_wus, world, internal = TRUE, overwrite = TRUE)

usethis::use_data(prism, cera, cesm, ccsm, prism_clim, geop_jfm, sst, states_wus, world, internal = TRUE, overwrite = TRUE)


## noaa
noaa <- read_ncdf('~/Downloads/weasd.mon.mean.nc', var = 'weasd') %>%
  filter(lubridate::month(time) == 3,
         dplyr::between(lubridate::year(time), 1901, 2010)) %>% # alternatively format(myDate,"%m")
  as('Raster') %>%
  raster::rotate() %>%
  raster::crop(bbox) %>%
  raster::mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_crop(states_wus)

## Bilinear interpolation


# raster::resample(cera[[99]] / mean(cera[[82:110]]), prism) %>% plot
# raster::resample((cera[[99]] / mean(cera[[82:110]])), prism) %>% plot
#
#
# filter(year >= 1982) %>%
#     group_by(x,y) %>%
#   mutate(SWE = (SWE - mean(SWE))/sd(SWE)) %>%
#   ungroup() %>%
#   filter(year == 1999) %>%
#   ggplot(aes(x, y)) +
#   geom_raster(aes(fill = SWE)) +
#   coord_quickmap() +
#   scale_fill_viridis_c(limits = c(-5,5))
#
# cera[cera < 0.001] <- 0
#
# delta_mult <- raster::resample(cera / mean(cera[[82:110]]) , prism) * mean(prism)
# delta_mult2 <- raster::resample(cera / mean(cera) , prism) * mean(prism)
#
# delta_add <- raster::resample(cera - mean(cera[[82:110]]) , prism) + mean(prism)
