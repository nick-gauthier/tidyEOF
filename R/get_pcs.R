
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

get_pcs <- function(dat, time_dim = 3, scale = FALSE, clim = NULL) {

  # get climatologies if not supplied
  # suppress warning about unequal dimensions
  if(is.null(clim)) clim <- suppressWarnings(get_climatology(dat))

  dat %>%
    `-`(clim['mean']) %>% # center the field
    {if(scale) . /  clim['sd'] else .} %>% # scale the field (optional)
    area_weight() %>% # weight by sqrt cosine latitude, in radians
    split(time_dim) %>% # split along the time dimension
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    t() %>% # transpose to space is columns and time rows
    prcomp(center = FALSE)
}


area_weight <- function(dat) {
  st_dim_to_attr(dat, which = 2) %>% # get the y coordinates
    `*`(pi / 180) %>% # convert to radians
    cos() %>% # cosine weighting
    sqrt() %>% # sqrt so the covariance matrix is weighted by cosine latitude
    `*`(dat)
}

#' @export
get_eigenvalues <- function(pca){
  n <- length(pca$sdev)

  pca %>%
    broom::tidy(matrix = 'pcs') %>%
    mutate(eigenvalues = std.dev ^ 2,
           error = sqrt(2 / n),
           low =  eigenvalues * (1 - error) * 100 / sum(eigenvalues),
           hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
           cumvar_line = hi + 0.02 * max(hi))
}

# older version for data frame, turn into method?
get_pcs2 <- function(dat, scale = FALSE){
  dat %>%
    dplyr::group_by(x,y) %>%
    dplyr::mutate(SWE = SWE - mean(SWE)) %>% # anomalize before weighting
    {if(scale) dplyr::mutate(., SWE = SWE / sd(SWE)) else .} %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>% # area weight, sqrt means covariance matrix weighted by cos(lat)
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    prcomp(center = FALSE)
}
