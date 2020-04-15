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

  data_in <- prep_data(patterns.x, patterns.y)

  model <- fit_model(data_in)

  data_out <- model %>%
    dplyr::select(PC, data) %>%
    unnest(cols = c(data)) %>%
    dplyr::select(PC, amplitude, year, pred, resid)

  # would it be worth it to also return the reconstructed spatial field here or not?
  coupled_patterns <- list(model = dplyr::select(model, -data),
                           data = data_out)
  class(coupled_patterns) <- 'coupled_patterns'
  return(coupled_patterns)
}

prep_data <- function(patterns.x, patterns.y) {
  amps.x <- if('patterns' %in% class(patterns.x)) patterns.x$amplitudes else patterns.x
  amps.y <- if('patterns' %in% class(patterns.y)) patterns.y$amplitudes else patterns.y

    predictors <- amps.x %>%
    mutate(PC = paste0('PC', PC)) %>%
    spread(PC, amplitude)

  predictands <- amps.y %>%
    mutate(PC = paste0('PC', PC)) %>%
    dplyr::select(year, PC, amplitude)

  inner_join(predictands, predictors, by = 'year')
}

fit_model <- function(data_in) {
  gam_formula <- data_in %>%
    dplyr::select(-PC, -amplitude) %>% #, -year) %>%
    names() %>%
    map( ~ paste0("s(", ., ", bs = 'cr', k = 3)")) %>%
    paste(collapse = ' + ') %>%
    paste('amplitude ~ ', .) %>%
    as.formula()

  model <- data_in %>%
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
}