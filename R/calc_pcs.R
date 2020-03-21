
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
calc_pcs <- function(dat, scale = TRUE){
  dat %>%
    dplyr::mutate(SWE = SWE * cos(y * pi / 180)) %>% # area weight
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    prcomp(scale. = scale)
}


