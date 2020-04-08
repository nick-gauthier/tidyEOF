#' Title
#'
#' @param coupled_patterns
#' @param new_amps
#'
#' @return
#' @export
#'
#' @examples
predict_patterns <- function(coupled_patterns, new_amps) {

  amps <- if('patterns' %in% class(new_amps)) new_amps$amplitudes else new_amps
  predictors <- amps %>%
    mutate(PC = paste0('PC', PC)) %>%
    spread(PC, amplitude)

  coupled_patterns$model$mod %>%
    map(~add_predictions(predictors, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(year, PC, amplitude)
}

reconstruct_field <- function(target_patterns, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  eofs <- target_patterns$eofs

  left_join(amplitudes, eofs, by = c('PC' = 'EOF')) %>%
    mutate(PC = as.numeric(PC)) %>%
    mutate(anomaly = weight * amplitude) %>%
    dplyr::select(-c(amplitude, value, weight)) %>%
    group_by(x, y, year) %>%
    summarize(anomaly = sum(anomaly)) %>%
    left_join(target_patterns$climatology, by = c("x", "y")) %>%
    mutate(SWE = anomaly + swe_mean,
           SWE = if_else(SWE < 0, 0, SWE))%>%
    dplyr::select(-c(swe_mean, swe_sd))
}