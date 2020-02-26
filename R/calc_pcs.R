
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
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    prcomp(scale. = scale)
}

area_weight <- function(x){
  names_x <- names(x)
  x %>%
    init('y') %>% # get a map of latitudes
    `*`(pi/180) %>% # convert to radians
    cos %>% # cosine
    sqrt %>%
    `*`(x) %>%
    `names<-`(names_x)
}
