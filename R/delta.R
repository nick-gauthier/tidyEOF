#' Delta change methods
#'
#' @param pred
#' @param obs
#' @param newdata
#' @param k A placeholder, must be NULL.
#'
#' @return
#' @export
#'
#' @examples
delta_mul <- function(pred, obs, newdata = NULL, k = NULL) { # k is just aplaceholder
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred)
  obs_clim <- get_climatology(obs)

  (newdata / pred_clim['mean']) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    `*`(obs_clim['mean'])
}

#' @export
delta_add <- function(pred, obs, newdata = NULL, k = NULL) { # k is just aplaceholder
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred)
  obs_clim <- get_climatology(obs)

  # simplify units here!
  (units::drop_units(newdata) - pred_clim['mean']) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    units::drop_units() %>%
    `+`(obs_clim['mean']) %>%
    # make varname generic!
    mutate(SWE = units::set_units(if_else(SWE < 0, 0, SWE), mm))
}