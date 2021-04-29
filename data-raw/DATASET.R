## code to prepare `DATASET` dataset goes here

library(raster) # processing raster data
library(tidyverse) # data manipulation and visualization
library(sf)

# download state boundary files for cropping and plotting
state_names <- c('arizona', 'new mexico', 'colorado',
                 'california', 'utah', 'nevada',
                 'oregon', 'washington', 'idaho',
                 'wyoming', 'montana')

states_wus <- maps::map('state', regions = state_names,
                        fill = TRUE, plot = FALSE) %>%
  st_as_sf()

world <- maps::map('world',wrap=c(0,360),fill = TRUE, plot = FALSE)


## Geographic Data

#Define a study area to constrain all computations.
bbox <- extent(c(-125, -102, 31, 49))


#Import the snow observation data from https://nsidc.org/data/nsidc-0719.^[What's up with this warning message? long_name=CRS definition
#spatial_ref=GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.01745329251994328,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]]
#GeoTransform=-125.0208 0.04166662697178698 0 49.9375 0 -0.04166662697178698].

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

prism <- list.files('../data', pattern = 'SWE_Depth', full.names = TRUE) %>%
  map(preprocess, var = 'SWE', bbox = bbox, regrid = TRUE, daily = TRUE) %>%
  brick() %>%
  mask(., all(near(., 0)), maskvalue = 1) %>% # mask pixels that never receive snow
  mask(states_wus) %>% # crop to state boundaries
  setNames(1982:2017)

cera <- map(1:10, ~ (brick('../data/CERA-20C_snow.nc', varname = 'sd', level = .))) %>%
  # this averages over the full ensemble
  reduce(`+`) %>%
  `/`(10) %>%
  crop(bbox) %>%
  mask(., all(near(., 0)), maskvalue = 1) %>%
  mask(states_wus) %>%
  `*`(1000) # convert to mm

cesm_h2osno <- preprocess('../data/b.e11.BLMTRC5CN.f19_g16.001.clm2.h0.H2OSNO.085001-184912.nc',
                          var = 'H2OSNO',
                          flip = TRUE,
                          bbox = bbox)

cesm_h2osno_ext <- preprocess('../data/b.e11.BLMTRC5CN.f19_g16.001.clm2.h0.H2OSNO.185001-200512.nc',
                          var = 'H2OSNO',
                          flip = TRUE,
                          bbox = bbox)

cesm <- c(cesm_h2osno, cesm_h2osno_ext) %>%
  brick() %>%
  raster::resample(cera) %>%
  mask(mean(cera))

# fix edge effects from resampling
cesm[cesm < 0] <- 0

ccsm_lm <- preprocess('../../snow-wna/data/snw_LImon_CCSM4_past1000_r1i1p1_085001-185012.nc',
                   var = 'snw',
                   flip = TRUE,
                   bbox = bbox)

ccsm_ext <- preprocess('~/Downloads/snw_LImon_CCSM4_historical_r1i2p1_185001-200512.nc',
                   var = 'snw',
                   flip = TRUE,
                   bbox = bbox)[[-1]] # the datasets overlap in year 1850

ccsm <- c(ccsm_lm, ccsm_ext) %>%
  brick() %>%
  raster::resample(cera) %>%
  mask(mean(cera))

# fix edge effects from resampling
ccsm[ccsm < 0] <- 0

#plot(mean(prism));plot(mean(cera));plot(mean(cesm));plot(mean(ccsm))

#climate data for teleconnection analysis

ppt <- list.files('data-raw/PRISM_ppt_stable_4kmM3_198101_201904_bil/', full.names = TRUE, pattern = '.bil$') %>%
  map(~raster(.) %>% crop(bbox)) %>%
  brick %>%
  aggregate(fact = 2) %>%
  crop(states_wus)

ppt_dat1 <- names(ppt) %>%
  str_split('_') %>%
  map_chr(~.[[5]]) %>%
  setNames(ppt, .) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(time = parse_number(layer)) %>%
  separate(time, into = c('year', 'month'), sep = -2, convert = TRUE) %>%
  select(-layer)

ppt_jfm <- ppt_dat1 %>%
  filter(month %in% c(1,2,3)) %>%
  group_by(x, y, year) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>% # diff from temp
  ungroup()

ppt_dat <- ppt_dat1 %>%
  filter(month %in% c(10,11,12)) %>%
  group_by(x, y, year) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>% # diff from temp
  ungroup() %>%
  mutate(year = year + 1) %>%
  inner_join(ppt_jfm, by = c("x", "y", "year")) %>%
  mutate(value = value.x + value.y, .keep = 'unused')

tmean <- list.files('data-raw/PRISM_tmean_stable_4kmM3_198101_201904_bil', full.names = TRUE, pattern = '.bil$') %>%
  map(~raster(.) %>% crop(bbox)) %>%
  brick %>%
  aggregate(fact = 2) %>%
  crop(states_wus)

tmean_dat1 <- names(tmean) %>%
  str_split('_') %>%
  map_chr(~.[[5]]) %>%
  setNames(tmean, .) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(time = parse_number(layer)) %>%
  separate(time, into = c('year', 'month'), sep = -2, convert = TRUE) %>%
  select(-layer)

tmean_jfm <- tmean_dat1 %>%
  filter(month %in% c(1,2,3)) %>%
  group_by(x, y, year) %>%
  summarise(value = mean(value, na.rm = TRUE)) %>% # diff from precip
  ungroup()

tmean_dat <- tmean_dat1 %>%
    filter(month %in% c(10,11,12)) %>%
    group_by(x, y, year) %>%
    summarise(value = mean(value, na.rm = TRUE)) %>% # diff from precip
    ungroup() %>%
    mutate(year = year + 1) %>%
    inner_join(tmean_jfm, by = c("x", "y", "year")) %>%
    mutate(value = (value.x + value.y) / 2, .keep = 'unused')

#sst_mean <- mean(brick('~/gdrive/Projects/snow-wna/data/sst.mon.ltm.1981-2010.nc', varname = 'sst')[[1:3]])

sst_brick <- brick('data-raw/sst.mnmean.nc', varname = 'sst')

time_sst <- getZ(sst_brick)
id_jfm <- which(time_sst >= as.Date('1982-01-01') & time_sst <= as.Date('2017-03-01') & str_detect(time_sst, c('-01-', '-02-', '-03-')))

sst_jfm <- subset(sst_brick, id_jfm) %>%
  stackApply(rep(1:36, each = 3), 'mean') %>%
 # `-`(sst_mean) %>% # is this necessary?
  crop(extent(c(-1, 359, -75, 75))) %>%
  disaggregate(fact = 2, method = 'bilinear') %>%
  setNames(1982:2017) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = parse_number(layer)) %>%
  select(-layer)

id_ond <- which(time_sst >= as.Date('1981-10-01') & time_sst <= as.Date('2017-03-01') & str_detect(time_sst, c('-10-', '-11-', '-12-')))

sst_ond <- subset(sst_brick, id_ond) %>%
  stackApply(rep(1:36, each = 3), 'mean') %>%
  # `-`(sst_mean) %>% # is this necessary?
  crop(extent(c(-1, 359, -75, 75))) %>%
  disaggregate(fact = 2, method = 'bilinear') %>%
  setNames(1982:2017) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = parse_number(layer)) %>%
  select(-layer)

sst_dat <- inner_join(sst_jfm, sst_ond, by = c("x", "y", "year")) %>%
  mutate(value = (value.x + value.y) / 2, .keep = 'unused')

geop_jfm <- brick('data-raw/adaptor.mars.internal-1584737536.594777-6445-11-8d0fcc69-bd7d-40e2-a423-cf86c741ac79.nc')[[-(1:3)]] %>%
   stackApply(rep(1:36, each = 3), 'mean') %>%
 #   `-`(., mean(.)) %>% # is this necessary?
  aggregate(fact = 4, method = 'bilinear') %>%
  setNames(1982:2017) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = parse_number(layer)) %>%
  select(-layer)

#geop_ond <- brick('data-raw/adaptor.mars.internal-1584737536.594777-6445-11-8d0fcc69-bd7d-40e2-a423-cf86c741ac79.nc')[[-(1:3)]] %>%
  stackApply(rep(1:36, each = 3), 'mean') %>%
  #   `-`(., mean(.)) %>% # is this necessary?
  aggregate(fact = 4, method = 'bilinear') %>%
  setNames(1982:2017) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = parse_number(layer)) %>%
  select(-layer)

# supporting data
landfrac <- raster('../data/_grib2netcdf-webmars-public-svc-blue-003-6fe5cac1a363ec1525f54343b6cc9fd8-tyqPrq.nc')

areas_prism <- area(prism) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(area = layer)

areas_cera <- (landfrac *  raster::area(cera)) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(area = layer)

mountains <- read_sf('../data/ne_10m_geography_regions_polys.shp') %>%
  filter(name %in% c('SIERRA NEVADA', 'CASCADE RANGE', 'ROCKY MOUNTAINS')) %>%
  dplyr::select(name, geometry)

mountains_mask <- mountains %>%
  mutate(rasters = split(., 1:3) %>% map(~mask(prism[[1]], .))) %>%
  st_drop_geometry() %>%
  mutate(rasters = map(rasters, as.data.frame, na.rm = TRUE, xy = TRUE)) %>%
  unnest(rasters) %>%
  dplyr::select(-X1982)

#Turn the raster bricks into data frames and join.
prism_dat <- prism %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = round(parse_number(layer))) %>%
  rename(SWE = value) %>%
  dplyr::select(-layer)

cera_dat <- cera %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = round(parse_number(layer))) %>%
  rename(SWE = value) %>%
  dplyr::select(-layer)

cesm_dat <- cesm %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = str_sub(layer, 2) %>% parse_number() %>% round) %>%
  rename(SWE = value) %>%
  select(-layer)

ccsm_dat <- ccsm %>%
  as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
  mutate(year = str_sub(layer, 2) %>% parse_number() %>% round) %>%
  rename(SWE = value) %>%
  select(-layer)

usethis::use_data(prism_dat, cera_dat, cesm_dat, ccsm_dat, tmean_dat, ppt_dat, geop_jfm, sst_dat, areas_prism, areas_cera, states_wus, world, internal = TRUE, overwrite = TRUE)


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
