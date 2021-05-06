
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

  eofs <- pca$rotation[, 1:k, drop = FALSE] # drop = FALSE preserves PC names when there's only 1 PC

    if(rotate == TRUE) {
      reofs <- varimax(eofs %*% diag(pca$sdev, k, k))

      rotation_matrix <- reofs$rotmat # save the rotation matrix for the amplitudes

      eofs <- unclass(reofs$loadings) %>%
        `colnames<-`(paste0('PC', 1:k))
    } else {
      rotation_matrix <- NULL
    }

  eof_maps <- dat[,,,1] %>%
    as_tibble() %>%
    na.omit() %>%
    dplyr::select(x, y) %>%
  bind_cols(as_tibble(eofs)) %>%
    st_as_stars() %>%
    st_set_crs(st_crs(dat)) %>%
    mutate(dummy = 1) %>% # hacky way to get around 1 pc issue bellow
    merge(name = 'PC') %>% # the problem with this is that it doesn't work if there is only 1 pc!
    .[,,,1:k]
    #slice('PC', 1:k)

  list(eofs = eof_maps,
       rotation_matrix = rotation_matrix)
}