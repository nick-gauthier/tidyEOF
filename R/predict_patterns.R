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
