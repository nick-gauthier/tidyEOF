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
  if(is(amplitudes, 'stars')) amplitudes <- project_patterns(target_patterns, amplitudes)
  # is there a more robust way to do nonneg?
  # check (ncol(amplitudes) - 1) == number of PCs in eofs?
  # check margin 3 is time?

  # this should work on new -raw- data and get the pc projection like cca

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

  final <- restore_climatology(anomalies,
                               clim = target_patterns$climatology,
                               scale = target_patterns$scaled,
                               monthly = target_patterns$monthly)

  if(target_patterns$weight) final <- final / lat_weights(final)
  if(nonneg) final <- mutate(final, across(everything(), ~if_else(.x < 0, 0, .x)))

  final %>%
    mutate(across(everything(), ~units::set_units(.x, target_patterns$units, mode = 'standard'))) # replace with modify2 for multiple variables?
}