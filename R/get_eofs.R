
#' Title
#'
#' @param dat
#' @param pca
#' @param k
#' @param eigenvalues
#' @param rotate
#'
#'scale by stdev (i.e. sqrt(eigenvalues)) for more robust rotation (Hannachi et al 2007)
#'
#' @return
#' @export
#'
#' @examples
#'
get_eofs <-  function(dat, pca, k, rotate = FALSE) {

  eofs <- pca$rotation[, 1:k]

    if(rotate == TRUE) {
      reofs <- varimax(eofs %*% diag(pca$sdev, k, k))

      rotation_matrix <- reofs$rotmat # save the rotation matrix for the amplitudes

      eofs <- unclass(reofs$loadings) %>%
        `colnames<-`(paste0('PC', 1:k))
    } else {
      rotation_matrix <- NULL
    }

  eof_maps <- dat[,,,1] %>%
    st_coordinates() %>%
    dplyr::select(x, y) %>%
  bind_cols(as_tibble(eofs)) %>%
    st_as_stars() %>%
    merge(name = 'PC') %>%
    slice('PC', 1:k)

  list(eofs = eof_maps,
       rotation_matrix = rotation_matrix)
}

# old version for data frames
get_eofs2 <-  function(dat, pca, k, eigenvalues, rotate = FALSE) {

  eofs <- pca %>%
    broom::tidy(matrix = 'variables') %>%
    filter(PC <= k) %>%
    left_join(eigenvalues[1:2], by = 'PC') %>%
    mutate(EOF = as.character(PC),
           column = as.character(column),
           # scale by stdev (i.e. sqrt(eigenvalues)) for more robust rotation
           weight = value * std.dev) %>%
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

  list(eofs = eof_maps,
       rotation_matrix = rotation_matrix)
}

