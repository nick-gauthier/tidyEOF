
#' Calculate principal components
#'
#' @param dat
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

get_pcs <- function(dat, scale = FALSE, clim = NULL, monthly = FALSE, weight = TRUE) {

  #if(weight) dat <- dat * lat_weights(dat) # weight by sqrt cosine latitude, in radians

  dat %>%
    get_anomalies(clim = clim, scale = scale, monthly = monthly) %>%
    {if(weight) . * lat_weights(.) else .} %>%
    split('time') %>% # split along the time dimension
    setNames(st_get_dimension_values(dat, 'time')) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>% # transpose to space is columns and time rows
    prcomp(center = FALSE)
}

#' @export
lat_weights <- function(dat) {
  lats <- st_dim_to_attr(dat, which = 2) # get the y coordinates -- this is brittle if not in x, y, time order
    # convert to radians then apply cosine weighting
    # sqrt so the covariance matrix is weighted by cosine latitude
  sqrt(cos(lats * pi / 180))
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
