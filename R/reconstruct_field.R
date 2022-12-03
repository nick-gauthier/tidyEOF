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
  setNames(target_patterns$names) %>%
  mutate(across(everything(), ~units::set_units(.x, target_patterns$units, mode = 'standard')))

if(target_patterns$monthly) {
  final <- anomalies %>%
    {if(target_patterns$scaled) sweep_months(., slice(target_patterns$climatology, 'var', 2), '*') else . } %>%
      sweep_months(slice(target_patterns$climatology, 'var', 1), '+')
} else {
  final <- anomalies %>%
    {if(target_patterns$scaled) . * slice(target_patterns$climatology, 'var', 2) else .} %>%
    `+`(slice(target_patterns$climatology, 'var', 1))
}
  final %>%
    units::drop_units() %>% # hacky doing this twice . . .
    {if(nonneg) mutate(., across(everything(), ~if_else(.x < 0, 0, .x))) else .} %>%
    mutate(across(everything(), ~units::set_units(.x, target_patterns$units, mode = 'standard'))) # replace with modify2 for multiple variables?
}
