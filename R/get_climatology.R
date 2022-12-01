#' Get annual or monthly means and standard deviations, or rescale a grid using climatology
#'
#' @param dat Input data.
#' @param monthly Calculate monthly climatology or anomalies instead of using
#' the entire period? Defaults to FALSE.
#' @param clim An optional climatology made with get_climatology().
#' Defaults to NULL, in which case get_climatology() is called internally.
#' @param scale Calculate standardized anomalies (x - mean) / sd? Defaults to
#' FALSE.
#'
#' Based on the `climatology()` and `scaleGrid()` functions from transformeR.
#' These functions apply calculations per grid cell.
#'
#' @return
#' @export
#'
#' @examples
get_climatology <- function(dat, monthly = FALSE) {
  if (monthly) {
    new_mn <- aggregate(dat, by_months, FUN = mean) %>%
      aperm(c(2, 3, 1)) %>% # aggregate puts time dimension first
      st_set_dimensions('geometry', names = 'month')

    new_sd <- aggregate(dat, by_months, FUN = sd) %>%
        aperm(c(2, 3, 1)) %>%
        st_set_dimensions('geometry', names = 'month')

  } else {
    new_mn <- st_apply(dat, 1:2, FUN = mean, na.rm = TRUE, rename = FALSE)
    new_sd <- st_apply(dat, 1:2, FUN = sd, na.rm = TRUE, rename = FALSE)
  }

  new <- c(mean = new_mn, stdev = new_sd, along = 'var')

  if (any(purrr::map_lgl(dat, inherits, 'units'))) {
    # do any of the attr. have units? if so, restore units
    new <-
      purrr::modify2(new, dat, ~ units::set_units(.x, units(.y), mode = 'standard')) %>%
      setNames(names(new)) %>%
      st_as_stars(dimensions = st_dimensions(new))
  }

  return(new)
}

#' @export
get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {

  if(is.null(clim)) {
    clim <- get_climatology(dat, monthly = monthly)
  }

  mn <- slice(clim, 'var', 1) # using slice instead of filter here because filter resets offset
  stdev <- slice(clim, 'var', 2)

  if(monthly) { # what if the months don't start with January or are uneven? does this still work?
    out <- sweep_months(dat, mn, '-')
    if(scale) {
      out <- sweep_months(out, stdev, '/')
    }

  } else {
    out <- dat - mn  # center the field
    if(scale) out <- out / stdev  ## scale the field (optional)
  }
  return(out)
}

# convenience function for monthly aggregation, based on example in aggregate.stars
by_months = function(x) {
  lubridate::month(x, label = TRUE, abbr = FALSE)
}

# convenience functions for sweeping monthly summary statistics.
sweep_months <- function(e1, e2, FUN) {
  FUN <- match.fun(FUN)
  purrr::map(1:12, ~ FUN(filter(e1, lubridate::month(time) == .x), adrop(filter(e2, month == month.name[.x])))) %>%
    do.call(c, .)
}
