
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

  pc_names <- names0(k, 'PC')
  eofs <- pca$rotation[, 1:k, drop = FALSE] %>% # drop = FALSE preserves PC names when there's only 1 PC
    `%*%`(diag(pca$sdev, k, k)) %>% # scale by sdev (sqrt(eigenvalues)) for more robust rotation
    `colnames<-`(pc_names)

    if(rotate == TRUE & k > 1) {
      reofs <- varimax(eofs)
      loadings <- unclass(reofs$loadings)
      rotation_matrix <- reofs$rotmat # save the rotation matrix to use later on the amplitudes


      colnames(rotation_matrix) <- pc_names
      eofs <- loadings
      colnames(eofs) <- pc_names

    } else {
      rotation_matrix <- NULL
    }

  # so you can't just merge by row becuase you may have lost values.
  # you could join by x and y but that'd be slow
  # but then you could use the old dimensions as is which is ideal
  ref <- dat[,,,1, drop = TRUE]

  y_dec <- st_get_dimension_values(ref, 'y') %>%
    is.unsorted(strictly = TRUE)

  eof_maps <- ref %>%
    as_tibble() %>%
    na.omit() %>% # can introduce issues if there are entire null rows/columns
    dplyr::select(x, y) %>%
  bind_cols(as_tibble(eofs)) %>%
    stars::st_as_stars(y_decreasing = y_dec) %>% # match ordering of y axis to original
    sf::st_set_crs(sf::st_crs(dat)) %>%
    mutate(dummy = 1) %>% # hacky way to get around 1 pc issue below
    merge(name = 'PC') %>% # the problem with this is that it doesn't work if there is only 1 pc!
    .[,,,1:k] %>%
    setNames('weight')
    #slice('PC', 1:k)

  list(eofs = eof_maps,
       rotation_matrix = rotation_matrix)
}
