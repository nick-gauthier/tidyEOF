#' Get EOFs and PCs from observations
#'
#' @param dat A `stars` object containing spatial and temporal dimensions.
#' @param k The number of PC/EOF modes to retain.
#'
#' @return A `patterns` object containing a tibble of PC amplitudes, a `stars`
#' object with EOF spatial patterns, and a `stars` object with the mean and
#' (optionally) standard deviation fields used.
#' @export
#'
#' @examples
#'
get_patterns <- function(dat, k = 4, scale = FALSE, rotate = FALSE, monthly = FALSE){

  climatology <- get_climatology(dat, monthly = monthly)
  pca <- get_pcs(dat, scale = scale, clim = climatology, monthly = monthly)
  eigenvalues <- get_eigenvalues(pca)
  eofs <- get_eofs(dat, pca, k, rotate)

  times <- st_get_dimension_values(dat, 3)

  amplitudes <- pca$x %>%
    .[,1:k, drop = FALSE] %>%
    scale() %>% # too strict? just divide by sqrt(eigenvalue)?
    {if(rotate & k > 1) . %*% eofs$rotation_matrix else .} %>%
    as_tibble(rownames = 'time') %>%
    mutate(time = as.numeric(time))
   # note that this object still has scale and center attributes

  patterns <- list(eofs = eofs$eofs,
       amplitudes = amplitudes,
       climatology = climatology,
       pca = pca,
       eigenvalues = eigenvalues,
       rotation = if(rotate & k > 1) eofs$rotation_matrix else NA,
       units = units(dat[[1]]), # only units of the 1st dataset
       names = names(dat),
       scaled = scale,
       monthly = monthly)

  class(patterns) <- 'patterns'
  return(patterns)
}
