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

  (newdata / pred_clim) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    `*`(obs_clim)
}

#' @export
delta_add <- function(pred, obs, newdata = NULL, monthly = NULL) {
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred, monthly = monthly)
  obs_clim <- get_climatology(obs, monthly = monthly)

  (newdata - pred_clim) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    `+`(obs_clim)
}