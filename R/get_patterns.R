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
get_patterns <- function(dat, k = 4, scale = FALSE, rotate = FALSE){

  climatology <- get_climatology(dat)
  pca <- get_pcs(dat, scale = scale, clim = climatology)
  eigenvalues <- get_eigenvalues(pca)
  eofs <- get_eofs(dat, pca, k, rotate)

  times <- st_get_dimension_values(dat, 3)

  amplitudes <- pca$x %>%
    .[,1:k, drop = FALSE] %>%
    scale() %>%
    {if(rotate == TRUE) . %*% eofs$rotation_matrix else .} %>%
    as_tibble(rownames = 'time') %>%
    mutate(time = as.numeric(time))
   # note that this object still has scale and center attributes

# eofs_corr <- amplitudes %>%
#    mutate(amplitude = if_else(PC %in% c('1','2'), amplitude * -1, amplitude)) %>% # change signs so physically interpretable
#   full_join(dat, by = 'year') %>%
#  group_by(x, y, PC) %>%
#  summarise(correlation = cor(amplitude, SWE))

  patterns <- list(eofs = eofs$eofs,
                  # eofs_corr = eofs_corr,
       amplitudes = amplitudes,
       climatology = climatology,
       pca = pca,
       eigenvalues = eigenvalues)

  class(patterns) <- 'patterns'
  return(patterns)
}

#' @export
get_climatology <- function(dat) {

 # unit <- units(dat[[1]])
  c(st_apply(dat, 1:2, mean, na.rm = TRUE),
    st_apply(dat, 1:2, sd, na.rm = TRUE)) #%>%
    # only works if there's one attribute
  # mutate(mean = units::set_units(mean, unit, mode = 'standard'),
   #       sd = units::set_units(sd, unit, mode = 'standard'))
}


# combine these?
#get_anomalies <- function(dat, scale = FALSE) {
#  dat %>%
#    dplyr::group_by(x,y) %>%
#    dplyr::mutate(SWE = SWE - mean(SWE)) %>% # anomalize before weighting
#   {if(scale) dplyr::mutate(., SWE = SWE / sd(SWE)) else .} %>%
#    dplyr::ungroup()
#}
