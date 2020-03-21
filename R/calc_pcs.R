
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
calc_pcs <- function(dat, scale = FALSE){
  dat %>%
    dplyr::group_by(x,y) %>%
    dplyr::mutate(SWE = SWE - mean(SWE)) %>% # anomalize before weighting
    {if(scale) dplyr::mutate(., SWE = SWE / sd(SWE)) else .} %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>% # area weight, sqrt means covariance matrix weighted by cos(lat)
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    prcomp()
}

