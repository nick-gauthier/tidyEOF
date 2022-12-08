#' Title
#'
#' @param target_patterns
#' @param amplitudes
#'
#' @return
#' @export
#'
#' @examples
reconstruct_field <- function(target_patterns, amplitudes = NULL, nonneg = TRUE) {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  # is there a more robust way to do nonneg?
  # check (ncol(amplitudes) - 1) == number of PCs in eofs?
  # check margin 3 is time?

  target_mean <- slice(target_patterns$climatology, 'var', 1) %>%
    units::drop_units()
  target_sd <- slice(target_patterns$climatology, 'var', 2) %>%
    units::drop_units()

  anomalies <- amplitudes %>%
    rowwise() %>%
    mutate(PCs = list(c_across(-time)), .keep = 'unused') %>%
    ungroup() %>%
    deframe() %>%
    purrr::map(~sweep(target_patterns$eofs, MARGIN = 3, STATS = .x, FUN = "*")) %>%
    do.call('c', .) %>%
    stars::st_apply(1:2, sum) %>%
    merge(name = 'time') %>%
    stars::st_set_dimensions('time', values = amplitudes$time) %>%
    setNames(target_patterns$names)

  if(target_patterns$monthly) {
      if(target_patterns$scaled) anomalies <- sweep_months(anomalies, target_sd, '*')
      final <- sweep_months(anomalies, target_mean, '+')
    } else {
      if(target_patterns$scaled) anomalies <- anomalies * target_sd
      final <- anomalies + target_mean
    }

  if(target_patterns$weight) final <- final / lat_weights(final)
  if(nonneg) final <- mutate(final, across(everything(), ~if_else(.x < 0, 0, .x)))

  final %>%
    mutate(across(everything(), ~units::set_units(.x, target_patterns$units, mode = 'standard'))) # replace with modify2 for multiple variables?
}