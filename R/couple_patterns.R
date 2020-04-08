#' Title
#'
#' @param patterns.x
#' @param patterns.y
#'
#' @return
#' @export
#'
#' @examples
couple_patterns <- function(patterns.x, patterns.y) {
  amps <- if('patterns' %in% class(patterns.x)) patterns.x$amplitudes else patterns.x
  predictors <- amps %>%
    mutate(PC = paste0('PC', PC)) %>%
    spread(PC, amplitude)

  predictands <- patterns.y$amplitudes %>%
    mutate(PC = paste0('PC', PC))

  gam_formula <- names(predictors) %>%
    map( ~ paste0('s(', ., ', k = 4)')) %>%
    paste(collapse = ' + ') %>%
    paste('amplitude ~ ', .) %>%
    as.formula()

  # should check how many overlapping years there are during model fitting, print result, and error if none
  model <- inner_join(predictands, predictors, by = 'year') %>%
    group_by(PC) %>%
    nest() %>%
    mutate(mod = map(data, ~ gam(
      gam_formula,
      data = .,
      method = 'REML',
      select = TRUE
    ))) %>%
    ungroup() %>%
    mutate(r2 = map_dbl(mod, ~ summary(.)$r.sq)) %>%
    mutate(
      data = map2(data, mod, add_predictions),
      data = map2(data, mod, add_residuals)
    )

  data <- model %>%
    dplyr::select(PC, data) %>%
    unnest(cols = c(data)) %>%
    dplyr::select(PC, amplitude, year, pred, resid)

  # would it be worth it to also return the reconstructed spatial field here or not?
  coupled_patterns <- list(model = dplyr::select(model, -data),
                           data = data)
  class(coupled_patterns) <- 'coupled_patterns'
  return(coupled_patterns)
}