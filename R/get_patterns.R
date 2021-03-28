#' Get EOFs and PCs from observations
#'
#' @param dat
#' @param k
#'
#' @return
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
    .[,1:k] %>%
    scale() %>%
    {if(rotate == TRUE) . %*% eofs$rotation_matrix else .} %>%
    as_tibble(rownames = 'time') %>%
    pivot_longer(-time, names_to = 'PC', values_to = 'amplitude')
    #mutate(time = as.numeric(times) )# need to tell if time is numeric or not
   # note that this object still has scale and center attributes

# eofs_corr <- amplitudes %>%
#    mutate(amplitude = if_else(PC %in% c('1','2'), amplitude * -1, amplitude)) %>% # change signs so physically interpetable
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
  c(st_apply(dat, 1:2, mean),
    st_apply(dat, 1:2, sd, na.rm = TRUE))
}

#older version for data frames
get_climatology2 <- function(dat) {
  dat %>%
  group_by(x, y) %>%
  summarise(swe_mean = mean(SWE), # make generic (i.e. var rather than swe)
            swe_sd = sd(SWE)) %>%
    ungroup()
}