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

get_patterns <- function(dat, k, mask = NULL, scale = FALSE, rotate = FALSE){

  if(!is.null(mask)) dat <- semi_join(dat, mask, by = c("x", "y"))

  climatology <- get_climatology(dat)
  pca <- get_pcs(dat, scale = scale)
  eigenvalues <- get_eigenvalues(pca)
  eofs <- get_eofs(dat, pca, k, eigenvalues, rotate)

  amplitudes <- pca$x %>%
    .[,1:k] %>%
    scale() %>%
    {if(rotate == TRUE) . %*% eofs$rotation_matrix else .} %>%
    as_tibble(rownames = 'year', .name_repair = ~as.character(1:k)) %>%
    gather(PC, amplitude, -year) %>%
    mutate(year = as.numeric(year)) # note that this object still has scale and center attributes

  eofs_corr <- amplitudes %>%
    mutate(amplitude = if_else(PC %in% c('1','2'), amplitude * -1, amplitude)) %>% # change signs so physically interpetable
    full_join(dat, by = 'year') %>%
    group_by(x, y, PC) %>%
    summarise(correlation = cor(amplitude, SWE)) %>%
    {if(!is.null(mask)) semi_join(., mask, by = c("x", "y")) else .}

  patterns <- list(eofs = eofs$eofs,
                   eofs_corr = eofs_corr,
       amplitudes = amplitudes,
       climatology = climatology,
       pca = pca,
       eigenvalues = eigenvalues)

  class(patterns) <- 'patterns'
  return(patterns)
}


get_climatology <- function(dat) {
  dat %>%
  group_by(x, y) %>%
  summarise(swe_mean = mean(SWE), # make generic (i.e. var rather than swe)
            swe_sd = sd(SWE))
}