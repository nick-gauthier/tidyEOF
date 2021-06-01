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
    st_set_dimensions(., 'PC', values = paste0('PC', st_get_dimension_values(., 'PC'))) %>%
    aperm(c(2,3,1))
  )
}

#' @export
get_fdr <- function(dat, patterns, fdr = 0.1) {
  times_amps <- patterns$amplitudes$time
  times_dat <- st_get_dimension_values(dat, 'time')
  times_cor <- intersect(times_amps, times_dat)

  if(length(times_cor) <= 0) stop ('Need at least two time steps in common') # add to tests
  if(!(identical(times_cor, times_amps) & identical(times_cor, times_dat))) warning('Using the time period ', first(times_cor), ' to ', last(times_cor))

  amps <- filter(patterns$amplitudes, time %in% times_cor) %>%
    select(-time)

  suppressWarnings( # suppress warnings that sd is zero
   fdr_rast <- filter(dat, time %in% times_cor) %>%
      st_apply(c('x', 'y'), fdr_fun, amps = amps, .fname = 'PC') %>%
     aperm(c(2,3,1)) %>%
      st_apply('PC', adjust) %>%
    setNames('FDR')
  )

fdr_rast %>%
 st_get_dimension_values('PC') %>%
seq_along() %>%
map(~slice(fdr_rast, 'PC', .x) %>%
      st_contour(contour_lines = TRUE, breaks = fdr) %>%
      transmute(PC = paste0('PC', .x))) %>%
do.call(rbind, .)
}

fdr_fun <- function(x, amps) {
  if(!any(is.na(x))) apply(amps, 2, function(y) cor.test(x, y)$p.value) else rep(NA, ncol(amps))
}

adjust <- function(x) {
  p.adjust(x, method = 'fdr', n = sum(!is.na(x)))
}

