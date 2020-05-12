#' Title
#'
#' @param coupled_patterns
#' @param new_amps
#'
#' @return
#' @export
#'
#' @examples
predict_patterns <- function(coupled_patterns, predictors) {# new_amps) {

  preds <- if('patterns' %in% class(predictors)) predictors$amplitudes else predictors
  if('patterns' %in% class(predictors)) {
    preds <- preds %>%
      mutate(PC = paste0('PC', PC)) %>%
      spread(PC, amplitude)
  }

  models <- if('coupled_patterns' %in% class(coupled_patterns)) coupled_patterns$model else coupled_patterns
  models$mod %>%
    map(~add_predictions(preds, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(year, PC, amplitude)
}

reconstruct_field <- function(target_patterns, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  eofs <- target_patterns$eofs

  left_join(amplitudes, eofs, by = c('PC' = 'EOF')) %>%
    mutate(PC = as.numeric(PC)) %>%
    mutate(anomaly = weight * amplitude) %>%
    dplyr::select(-c(amplitude, weight)) %>%
    group_by(x, y, year) %>%
    summarize(anomaly = sum(anomaly)) %>%
    left_join(target_patterns$climatology, by = c("x", "y")) %>%
    mutate(SWE = anomaly + swe_mean,
           SWE = if_else(SWE < 0, 0, SWE))%>% # can't have negative swe
    dplyr::select(-c(swe_mean, swe_sd)) %>%
    ungroup()
}