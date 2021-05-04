#' Calculate teleconnections between a patterns object and another layer
#'
#'The function will check if there are any overlapping time steps between the
#'raster field and the PC amplitudes.
#' @param dat
#' @param patterns
#'
#' @return
#' @export
#'
#' @examples
get_correlation <- function(dat, patterns) {
  times_amps <- patterns$amplitudes$time
  times_dat <- st_get_dimension_values(dat, 'time')
  times_cor <- intersect(times_amps, times_dat)

  if(length(times_cor) <= 0) stop ('Need at least two time steps in common') # add to tests
  if(!(identical(times_cor, times_amps) & identical(times_cor, times_dat))) warning('Using the time period ', first(times_cor), ' to ', last(times_cor))

  amps <- filter(patterns$amplitudes, time %in% times_cor) %>%
    select(-time)

  suppressWarnings( # suppress warnings that sd is zero
  filter(dat, time %in% times_cor) %>%
  st_apply(c('x', 'y'), function(x) cor(x, amps), .fname = 'PC') %>%
    st_set_dimensions(., 'PC', values = paste0('PC', st_get_dimension_values(., 'PC')))
  )
}

# need to make sure time series are on same scale
#' @export
get_fdr <- function(dat, patterns, fdr = 0.1) { # could combine with above
  amps <- patterns$amplitudes %>%
    spread(PC, amplitude, sep ='')

  dat %>%
    group_by(x,y) %>%
    nest %>%
    mutate(corrs = map(data,  ~inner_join(., amps, by = 'year') %>%
                         summarise(PC1 = cor.test(value, PC1)$p.value,
                                   PC2 = cor.test(value, PC2)$p.value,
                                   PC3 = cor.test(value, PC3)$p.value,
                                   PC4 = cor.test(value, PC4)$p.value))) %>%
    select(-data) %>%
    unnest(corrs) %>%
    gather(PC, value, PC1:PC4) %>%
    group_by(PC) %>%
    mutate(value = p.adjust(value, method = 'fdr')) %>%
    filter(value < fdr) %>%
    ungroup()
}