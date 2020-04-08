#' Title
#'
#' @param patterns
#'
#' @return
#' @export
#'
#' @examples
reconstruct <- function(target_patterns, amplitudes = NULL) {
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