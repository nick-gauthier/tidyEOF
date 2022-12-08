#' Delta change methods
#'
#' @param pred
#' @param obs
#' @param newdata
#'
#' @return
#' @export
#'
#' @examples
delta_mul <- function(pred, obs, newdata = NULL, monthly = FALSE) {
  if(is.null(newdata)) newdata <- pred

  pred_clim <- get_climatology(pred, monthly = monthly) %>%
    slice('var', 1)
  obs_clim <- get_climatology(obs, monthly = monthly) %>%
    slice('var', 1)

  if(monthly) {
    sweep_months(newdata, pred_clim, '/') %>%
      delta_interp(obs, newdata) %>%
      sweep_months(obs_clim, '*')
  } else {
  (newdata / pred_clim) %>%
    delta_interp(obs, newdata) %>%
    `*`(obs_clim)
  }
}

#' @export
delta_add <- function(pred, obs, newdata = NULL, monthly = FALSE) {
  if(is.null(newdata)) newdata <- pred

  pred_clim <- get_climatology(pred, monthly = monthly) %>%
    slice('var', 1)
  obs_clim <- get_climatology(obs, monthly = monthly) %>%
    slice('var', 1)

  if(monthly) {
    sweep_months(newdata, pred_clim, '-') %>%
      delta_interp(obs, newdata) %>%
      sweep_months(obs_clim, '+')
  } else {
    (newdata - pred_clim) %>%
      delta_interp(obs, newdata) %>%
      `+`(obs_clim)
  }

}


#convenience function to interpolate to a higher grid and restore units and dimensions
delta_interp <- function(anom, obs, newdata) {
  # setup the new dimensions
  new_dims <- st_dimensions(obs)
  new_dims$time <- st_dimensions(newdata)$time # use time dimensions from newdata

  anom %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear', no_data_value = -99999) %>%
    setNames(names(newdata)) %>%
    # this only uses the first unit, should be able to support multiple units across attributes in the future
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    `st_dimensions<-`(new_dims) # restore original dimensions
}
