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
  eofs <- as_tibble(target_patterns$eofs)
  clim <- as_tibble(target_patterns$climatology)

  amplitudes %>%
    pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
    left_join(eofs, by = 'PC') %>%
    mutate(anomaly = weight * amplitude) %>%
    dplyr::select(-c(amplitude, weight)) %>%
    group_by(x, y, time) %>%
    summarize(anomaly = sum(anomaly, na.rm = TRUE), .groups = 'drop') %>%
    left_join(clim, by = c("x", "y")) %>%
    # should rewrite so SWE isn't hard coded
    mutate(SWE = anomaly + mean,
           SWE = if_else(SWE < 0, 0, SWE)) %>% # can't have negative swe
    dplyr::select(-c(mean, sd, anomaly)) %>%
    st_as_stars(dims = c('x','y','time')) %>%
    st_set_crs(st_crs(target_patterns$eofs)) %>%
    mutate(SWE = units::set_units(SWE, mm))
}

# loop through amplitudes, for each time row, multiply those numbers by the eofs, at the end join the results?
