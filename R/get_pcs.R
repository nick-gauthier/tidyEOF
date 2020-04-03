
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
get_pcs <- function(dat, scale = FALSE){
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

get_eigenvalues <- function(pca_object){
  n <- length(pca_object$sdev)

  pca_object %>%
    broom::tidy(matrix = 'pcs') %>%
    mutate(eigenvalues = std.dev ^ 2,
           error = sqrt(2 / n),
           low =  eigenvalues * (1 - error) * 100 / sum(eigenvalues),
           hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
           cumvar_line = hi + 0.02 * max(hi))
}

