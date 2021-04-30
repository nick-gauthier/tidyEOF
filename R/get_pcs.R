
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

get_pcs <- function(dat, scale = FALSE, clim = NULL) {

  # get climatologies if not supplied
  # suppress warning about unequal dimensions
  if(is.null(clim)) clim <- suppressWarnings(get_climatology(dat))

  dat %>%
    units::drop_units() %>%
    `-`(clim['mean']) %>% # center the field
    {if(scale) . /  clim['sd'] else .} %>% # scale the field (optional)
    area_weight() %>% # weight by sqrt cosine latitude, in radians
    split('time') %>% # split along the time dimension
    setNames(st_get_dimension_values(dat, 'time')) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>% # transpose to space is columns and time rows
    prcomp(center = FALSE)
}


area_weight <- function(dat) {
  st_dim_to_attr(dat, which = 2) %>% # get the y coordinates
    `*`(pi / 180) %>% # convert to radians
    cos() %>% # cosine weighting
    sqrt() %>% # sqrt so the covariance matrix is weighted by cosine latitude
    `*`(dat, .)
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
