
#' Title
#'
#' @param dat
#' @param pca
#' @param k
#' @param eigenvalues
#' @param rotate
#'
#' @return
#' @export
#'
#' @examples
get_eofs <- function(dat, pca, k, eigenvalues, rotate) {
  eofs <- pca %>%
    broom::tidy(matrix = 'variables') %>%
    filter(PC <= k) %>%
    left_join(eigenvalues[1:2], by = 'PC') %>%
    mutate(weight = value * std.dev, # scale by stdev (i.e. sqrt(eigenvalues)) for more robust rotation (hannachi et al 2007)
           EOF = as.character(PC),
           column = as.character(column)) %>%
    dplyr::select(-c(std.dev, PC))

  if(rotate == TRUE) {
    reofs <- eofs %>% # varimax rotation
      dplyr::select(-value) %>% # rotate on the weights
      pivot_wider(names_from = EOF, values_from = weight) %>%
      column_to_rownames(var = 'column') %>%
      as.matrix() %>%
      varimax()

    rotation_matrix <- reofs$rotmat # save the rotation matrix for the amplitudes

    eofs <- unclass(reofs$loadings) %>%
      as_tibble(rownames = 'column') %>%
      pivot_longer(-column, names_to = 'EOF', values_to = 'weight')
  } else {
    rotation_matrix <- NULL
  }

  eofs <- dat %>%
    spread(year, SWE) %>%
    mutate(column = as.character(1:n())) %>%
    dplyr::select(x, y, column) %>%
    full_join(eofs, by = 'column') %>%
    dplyr::select(-column)

  list(eofs = eofs,
       rotation_matrix = rotation_matrix)
}


