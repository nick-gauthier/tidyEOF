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

  new <- c(mean = new_mn, stdev = new_sd, along = 'var') %>%
    mutate()
   #for some reason c.stars(,, along = 'var') drops dimnames from the underlying array, mutate restores them

  if (any(purrr::map_lgl(dat, inherits, 'units'))) {
    # do any of the attr. have units? if so, restore units
    new <- restore_units(new, dat)
  }

  return(new)
}

#' @export
get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {

  if(is.null(clim)) clim <- get_climatology(dat, monthly = monthly)

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

#' @export
restore_climatology <- function(anomalies, clim, scale = FALSE, monthly = FALSE) {
  target_mean <- slice(clim, 'var', 1) %>%
    units::drop_units() # is this necessary
  target_sd <- slice(clim, 'var', 2) %>%
    units::drop_units()

  if(monthly) {
    if(scale) anomalies <- sweep_months(anomalies, target_sd, '*')

    final <- sweep_months(anomalies, target_mean, '+')
  } else {
    if(scale) anomalies <- anomalies * target_sd

    final <- anomalies + target_mean
  }
# should restore units here?
  return(final)
}

# convenience function for monthly aggregation, based on example in aggregate.stars
by_months = function(x) {
  lubridate::month(x, label = TRUE, abbr = FALSE)
}

# convenience functions for sweeping monthly summary statistics.
# should check that time is posix?
sweep_months <- function(e1, e2, FUN) {
  FUN <- match.fun(FUN)
  purrr::map(1:12, ~ FUN(filter(e1, lubridate::month(time) == .x), abind::adrop(filter(e2, month == month.name[.x])))) %>%
    do.call(c, .) %>%
    slice(., 'time', order(time(.))) %>% # reshuffle so months/years in right order
  st_set_dimensions(., 'time', values = time(.)) # the lubridate command above results in interval times not dates
}

# convenience function to add back units . . . there has to be a better way!
restore_units <- function(new, ref) {
  old_units <- purrr::map(ref, purrr::possibly(units))
  apply_units <- function(x, y) purrr::modify_in(x, names(old_units)[y], ~units::set_units(.x, old_units[[y]], mode = 'standard'))

  seq_along(new) %>%
    purrr::reduce(apply_units, .init = new)
}