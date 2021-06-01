## code to prepare `DATASET` dataset goes here

library(raster) # processing raster data
library(tidyverse) # data manipulation and visualization
library(stars)
sf_use_s2(TRUE) # use s2 spherical geometry

# download state boundary files for cropping and plotting
states_wus <- rnaturalearth::ne_states(country = 'United States of America', returnclass = 'sf') %>%
  filter(name %in% c('Arizona', 'New Mexico', 'Colorado',
                     'California', 'Utah', 'Nevada',
                     'Oregon', 'Washington', 'Idaho',
                     'Wyoming', 'Montana'))

d <- rnaturalearth::ne_countries(returnclass = "sf", scale = 'small') %>%
  st_set_precision(1e7) %>%
  st_as_s2() %>%
  st_as_sf()
x1 <- st_crop(d, st_bbox(c(xmin = 0, xmax = 180, ymin = -90, ymax = 90)))
x2 <- st_crop(d, st_bbox(c(xmin = -180, xmax = -0.001, ymin = -90, ymax = 90)))

x <- rbind(x1, x2)
## mollwiede on lon=180
prj <- "+proj=moll +lon_0=180 +ellps=WGS84 +datum=WGS84 +no_defs"
## fudge above means we get a seam in Antarctica and Russia ....
world <- sf::st_transform(x[1], prj)

## Geographic Data

#Define a study area to constrain all computations.
bbox <- extent(st_bbox(states_wus))# extent(c(-125, -102, 31, 49))


preprocess <- function(x, var, bbox, flip = FALSE, regrid = FALSE, daily = FALSE, month = '-03-') {
  maps <- brick(x, varname = var)

  indices <- getZ(maps) %>% # find the index for the time steps in question
    str_detect(month)%>%
    which()

  raster::subset(maps, indices) %>%
    {if(flip) rotate(.) else .} %>%
    crop(bbox, snap = 'out') %>%
    {if(daily) mean(.) else .} %>% # if a file is one year of daily values
    {if(regrid) aggregate(., fact = 2) else .} # resample to lower res to speed up analysis
}

prism <- list.files('data-raw/UA-SWE', pattern = 'SWE_Depth', full.names = TRUE) %>%
  map(preprocess, var = 'SWE', bbox = bbox, regrid = TRUE, daily = TRUE) %>%
  brick() %>%
  mask(., all(near(., 0)), maskvalue = 1) %>% # mask pixels that never receive snow
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 1982:2017, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm)) %>%
  st_set_crs(4326) %>% # the documentation for UA SWE says its EPSG 4326, but the file says EPSG 4269 -- not a huge difference at this scale so just transform without reprojecting
  st_crop(states_wus) %>%
  .[,2:273,1:212,] # remove blank cells at edges of domain

## cera
# the CERA data are in mm SWE for the land fraction of the grid cell. Need to correct for this
land_frac <- read_ncdf('data-raw/_grib2netcdf-webmars-public-svc-blue-003-6fe5cac1a363ec1525f54343b6cc9fd8-tyqPrq.nc',
                       var = 'lsm') %>%
  adrop(4) %>%
  slice('number', 1)

cera_raw <- map(1:10, ~ (brick('data-raw/CERA-20C_snow.nc', varname = 'sd', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  `*`(1000) %>% # convert from m water equivalent to mm
  st_as_stars() %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 1901:2010, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm)) %>%
  `*`(land_frac)

# mask out cells that never receive snow
cera_mask <- cera_raw %>%
  units::drop_units() %>%
  near(0) %>%
  st_apply(1:2, all) %>%
  transmute(all = as.numeric(na_if(!all, 0)))

cera <- (cera_raw * cera_mask) %>%
  st_crop(states_wus) %>%
  .[,2:23,2:18]

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
  #mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_warp(slice(cera, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm)) %>%
  st_crop(st_as_sf(cera[,,,1]))# make sure na cells in CERA are na here too

cesm_raw <- c(cesm_h2osno, cesm_h2osno_ext) %>%
  brick() %>%
  #mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
 # st_warp(slice(cera, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm)) #%>%
#  st_crop(st_as_sf(cera[,,,1]))

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
  #mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_warp(slice(cera, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
  setNames('SWE') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(SWE = units::set_units(SWE, mm)) %>%
  st_crop(st_as_sf(cera[,,,1]))# make sure na cells in CERA are na here too

ccsm_lm_prec <- 'data-raw/pr_Amon_CCSM4_past1000_r1i1p1_085001-185012.nc' %>%
  preprocess(var = 'pr', bbox = bbox, flip = TRUE,
             month = paste0('-', c('01', '02', '03', 10:12), '-')) %>%
  stackApply(rep(1:2002, each = 3), sum) %>%
  .[[-c(1, 2002)]] %>%
  stackApply(rep(1:1000, each = 2), sum)

# note the datasets overlap in year 1850
ccsm_ext_prec <- 'data-raw/pr_Amon_CCSM4_historical_r1i2p1_185001-200512.nc' %>%
  preprocess(var = 'pr', bbox = bbox, flip = TRUE,
             month = paste0('-', c('01', '02', '03', 10:12), '-')) %>%
  stackApply(rep(1:312, each = 3), sum) %>%
  .[[-c(1, 312)]] %>%
  stackApply(rep(1:155, each = 2), sum)

ccsm_prec <- c(ccsm_lm_prec, ccsm_ext_prec) %>%
  brick() %>%
  #mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_warp(slice(cera, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
  setNames('pr') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(pr = units::set_units(pr, mm)) %>%
  st_crop(st_as_sf(cera[,,,1])) # make sure na cells in CERA are na here too

ccsm_lm_temp <- 'data-raw/tas_Amon_CCSM4_past1000_r1i1p1_085001-185012.nc' %>%
  preprocess(var = 'tas', bbox = bbox, flip = TRUE,
             month = paste0('-', c('01', '02', '03', 10:12), '-')) %>%
  stackApply(rep(1:2002, each = 3), mean) %>%
  .[[-c(1, 2002)]] %>%
  stackApply(rep(1:1000, each = 2), mean)

# note the datasets overlap in year 1850
ccsm_ext_temp <- 'data-raw/tas_Amon_CCSM4_historical_r1i2p1_185001-200512.nc' %>%
  preprocess(var = 'tas', bbox = bbox, flip = TRUE,
             month = paste0('-', c('01', '02', '03', 10:12), '-')) %>%
  stackApply(rep(1:312, each = 3), mean) %>%
  .[[-c(1, 312)]] %>%
  stackApply(rep(1:155, each = 2), mean)

ccsm_temp <- c(ccsm_lm_temp, ccsm_ext_temp) %>%
  brick() %>%
  #mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_warp(slice(cera, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
  setNames('tas') %>%
  st_set_dimensions('band', values = 850:2005, names = 'time') %>%
  mutate(tas = units::set_units(tas, '°C')) %>%
  st_crop(st_as_sf(cera[,,,1]))# make sure na cells in CERA are na here too

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
# why does ppt have no nas but tmean does?
prism_clim <- c(ppt_jfm[,-1] + ppt_ond,
                (tmean_jfm[,-1] + tmean_ond) / 2) %>%
  mutate(ppt = units::set_units(ppt, mm),
         tmean = units::set_units(tmean, '°C')) %>%
  st_warp(prism) %>%
  st_crop(states_wus)

# CERA-20C geopotential
geop <- map(1:10, ~ (brick('data-raw/CERA-20C_geop.nc', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  stackApply(rep(1:29, each = 6), 'mean') %>%
 # shift(180) %>%
  #rotate() %>%
  #shift(180) %>%
  st_as_stars() %>%
  st_set_dimensions('band', values = 1982:2010, names = 'time') #%>%
 # st_set_crs(4326)

# CERA-20C sst
land_mask <- read_ncdf('~/Downloads/_grib2netcdf-webmars-public-svc-blue-000-6fe5cac1a363ec1525f54343b6cc9fd8-6_792e.nc',
                       var = 'lsm') %>%
  adrop(4) %>%
  slice('number', 1) %>%
  mutate(lsm = if_else(lsm < 0.5, 1, NA_real_))

sst <- map(1:10, ~ (brick('data-raw/CERA-20C_sst.nc', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  stackApply(rep(1:29, each = 6), 'mean') %>%
  #shift(180) %>%
  #rotate() %>%
  #shift(180) %>%
  st_as_stars() %>%
  st_set_dimensions('band', values = 1982:2010, names = 'time') %>%
  `*`(land_mask)
 # st_set_crs(4326)

sst_cera_jfm <- map(1:10, ~ (brick('data-raw/CERA-20C_sst_jfm.nc', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  stackApply(rep(1:29, each = 3), 'mean') %>%
  #shift(180) %>%
  #rotate() %>%
  #shift(180) %>%
  st_as_stars() %>%
  st_set_dimensions('band', values = 1982:2010, names = 'time') %>%
  `*`(land_mask)

geop_cera_jfm <- map(1:10, ~ (brick('data-raw/CERA-20C_geop_jfm.nc', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  stackApply(rep(1:29, each = 3), 'mean') %>%
  # shift(180) %>%
  #rotate() %>%
  #shift(180) %>%
  st_as_stars() %>%
  st_set_dimensions('band', values = 1982:2010, names = 'time') #%>%
# st_set_crs(4326)

#mountains <- read_sf('../data/ne_10m_geography_regions_polys.shp') %>%
#  filter(name %in% c('SIERRA NEVADA', 'CASCADE RANGE', 'ROCKY MOUNTAINS')) %>%
#  dplyr::select(name, geometry)

#mountains_mask <- mountains %>%
#  mutate(rasters = split(., 1:3) %>% map(~mask(prism[[1]], .))) %>%
#  st_drop_geometry() %>%
#  mutate(rasters = map(rasters, as.data.frame, na.rm = TRUE, xy = TRUE)) %>%
#  unnest(rasters) %>%
#  dplyr::select(-X1982)

usethis::use_data(prism, cera, cera_raw, cesm, ccsm, ccsm_prec, ccsm_temp, prism_clim, geop, sst, noaa, states_wus, world, internal = TRUE, overwrite = TRUE)

## noaa
noaa <- read_ncdf('~/Downloads/weasd.mon.mean.nc', var = 'weasd') %>%
  filter(lubridate::month(time) == 3) %>%
      #   dplyr::between(lubridate::year(time), 1901, 2010)) %>% # alternatively format(myDate,"%m")
  as('Raster') %>%
  raster::rotate() %>%
  raster::crop(bbox) %>%
  raster::mask(., all(near(., 0)), maskvalue = 1) %>%
  st_as_stars() %>%
  st_crop(states_wus) %>%
  setNames('SWE') %>%
  mutate(SWE = units::set_units(SWE, mm)) %>%
  st_set_dimensions('band', values = 1836:2015, names = 'time')

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
