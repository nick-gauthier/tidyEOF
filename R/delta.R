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
  pred_clim <- get_climatology(pred, monthly = monthly)
  obs_clim <- get_climatology(obs, monthly = monthly)

  # setup the new dimensions
  new_dims <- st_dimensions(obs)
  new_dims$time <- st_dimensions(newdata)$time # use time dimensions from newdata

  (newdata / pred_clim) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear', no_data_value = -99999) %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    `st_dimensions<-`(new_dims) %>% # restore original dimensions
    `*`(obs_clim)
}

#' @export
delta_add <- function(pred, obs, newdata = NULL, monthly = FALSE) {
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred, monthly = monthly)
  obs_clim <- get_climatology(obs, monthly = monthly)

  # setup the new dimensions
  new_dims <- st_dimensions(obs)
  new_dims$time <- st_dimensions(newdata)$time # use time dimensions from newdata

  (newdata - pred_clim) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear', no_data_value = -99999) %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    `st_dimensions<-`(new_dims) %>% # restore original dimensions
    `+`(obs_clim)
}
