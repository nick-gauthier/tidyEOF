#' Get annual or monthly climatologies, or rescale a grid using climatologies
#'
#' @param dat Input data.
#' @param monthly Calculate mean and standard deviation for each month instead
#' of the entire period? Defaults to FALSE.
#'
#' Based on the `climatology()` and `localScaling()` functions from transformeR.
#' These functions apply calculations per gridcell.
#'
#' @return
#' @export
#'
#' @examples
get_climatology <- function(dat, monthly = FALSE) {

  if(monthly) {

    mon_mean <- aggregate(dat, by_months, mean) %>%
      aperm(c(2,3,1)) %>% # aggregate puts time dimension first
      st_set_dimensions('geometry', names = 'month') # bug(?) in stars names new time dimension "geometery"

    mon_sd <- aggregate(dat, by_months, sd) %>%
      aperm(c(2,3,1)) %>%
      st_set_dimensions('geometry', names = 'month')

    c(mon_mean, mon_sd) %>%
      setNames(c('mean', 'sd'))
  } else {
    # unit <- units(dat[[1]])
    c(st_apply(dat, 1:2, mean, na.rm = TRUE),
      st_apply(dat, 1:2, sd, na.rm = TRUE)) #%>%
    # only works if there's one attribute
    # mutate(mean = units::set_units(mean, unit, mode = 'standard'),
    #       sd = units::set_units(sd, unit, mode = 'standard'))

  }
}

#' @export
get_anomalies <- function(dat, monthly = TRUE) {

}

# convenience function for monthly aggregation, based on example in aggregate.stars
by_months = function(x) {
  lubridate::month(x, label = TRUE, abbr = FALSE)
}



test <- read_ncdf('data-raw/sst.mnmean.nc')
get_climatology(test[,,,1:100], monthly = TRUE)['sd'] %>% plot
