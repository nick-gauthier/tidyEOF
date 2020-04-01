#' Get PCs
#'
#' @param pc_object
#' @param k
#'
#' @return
#' @export
#'
#' @examples
get_amps <- function(pc_object, k){
  broom::tidy(pc_object, 'samples') %>%
    dplyr::mutate(year = as.numeric(as.character(row))) %>%
    dplyr::select(-row) %>%
    dplyr::filter(PC <= k) %>%
    dplyr::mutate(PC = as.character(PC)) %>%
    dplyr::rename(amplitude = value)
}

plot_amps <- function(patterns) {
  patterns %>%
    dplyr::select(-patterns) %>%
    unnest(amplitudes) %>%
    mutate(PC = paste0('PC', PC)) %>%
    ggplot(aes(year, amplitude)) +
    geom_line() +
    facet_wrap(~PC) +
    theme_bw()
}
