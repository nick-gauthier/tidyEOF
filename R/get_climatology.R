#' Get annual or monthly climatologies, or rescale a grid using climatologies
#'
#' @param dat Input data.
#' @param fun Function to calculate climatologies. Defaults to mean, but can also be sd.
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
get_climatology <- function(dat, fun = mean, monthly = FALSE) {
  if(monthly) {
    new <- aggregate(dat, by_months, fun) %>%
      aperm(c(2,3,1)) %>% # aggregate puts time dimension first
      st_set_dimensions('geometry', names = 'month')
  } else {
    new <- st_apply(dat, 1:2, fun, na.rm = TRUE)
  }
  # restore units that were stripped by st_apply
  purrr::modify2(new, dat, ~units::set_units(.x, units(.y), mode = 'standard'))
}

#' @export
# get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {
#   if(is.null(clim)) clim <- get_climatology(dat, monthly = monthly)
#
#   if(monthly) { # what if the months don't start with january or are uneven? does this still work?
#     if(scale) {
#       purrr::map(1:12, ~(filter(dat, lubridate::month(time) == .x) - filter(clim['mean'], month == .x)) / filter(clim['sd'], month == .x))
#     } else {
#       purrr::map(1:12, ~filter(dat, lubridate::month(time) == .x) - filter(clim['mean'], month == .x))
#     }
#
#   } else {
#     out <- dat  -  clim['mean']  # center the field
#     if(scale) return(out / clim['sd']) else return(out)  ## scale the field (optional)
#   }
# }

# convenience function for monthly aggregation, based on example in aggregate.stars
by_months = function(x) {
  lubridate::month(x, label = TRUE, abbr = FALSE)
}







