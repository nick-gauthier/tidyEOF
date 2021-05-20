#' Title
#'
#' @param target_patterns
#' @param amplitudes
#'
#' @return
#' @export
#'
#' @examples
reconstruct_field <- function(target_patterns, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes

  # check (ncol(amplitudes) - 1) == number of PCs in eofs?
  # check margin 3 is time?

amplitudes %>%
  rowwise() %>%
  mutate(PCs = list(c_across(-time)), .keep = 'unused') %>%
  ungroup() %>%
  deframe() %>%
  map(~sweep(target_patterns$eofs, MARGIN = 3, STATS = .x, FUN = "*")) %>%
  do.call('c', .) %>%
  st_apply(1:2, sum) %>%
  merge(name = 'time') %>%
  st_set_dimensions('time', values = amplitudes$time) %>%
  setNames('SWE') %>%
  {. +  target_patterns$climatology['mean']} %>%
  # should make generic!
  mutate(SWE = if_else(SWE < 0, 0, SWE),
         SWE = units::set_units(SWE, mm)) # can't have negative swe
}




