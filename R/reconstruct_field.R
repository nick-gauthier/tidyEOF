#' Title
#'
#' @param target_patterns
#' @param amplitudes
#'
#' @return
#' @export
#'
#' @examples
reconstruct_field <- function(target_patterns, amplitudes = NULL, nonneg = TRUE, scale = FALSE) {
  if(is.null(amplitudes)) amplitudes <- target_patterns$amplitudes
  # is there a more robust way to do nonneg?
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
  setNames(target_patterns$names) %>%
  `+`(target_patterns$climatology['mean']) %>%
  {if(scale) . * target_patterns$climatology['sd'] else .} %>%
  {if(nonneg) mutate(., across(everything(), ~if_else(.x < 0, 0, .x))) else .} %>%
  mutate(across(everything(), ~units::set_units(.x, target_patterns$units, mode = 'standard')))
}



