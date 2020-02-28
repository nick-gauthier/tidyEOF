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
  tidy(pc_object, 'samples') %>%
    mutate(year = as.numeric(as.character(row))) %>%
    dplyr::select(-row) %>%
    filter(PC <= k) %>%
    mutate(PC = as.character(PC)) %>%
    rename(amplitude = value)
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
