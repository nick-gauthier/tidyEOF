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


get_patterns <- function(dat, k, mask = NULL, scale = FALSE, rotate = TRUE){

  if(!is.null(mask)) dat <- semi_join(dat, mask, by = c("x", "y"))

  climatology <- dat %>%
    group_by(x, y) %>%
    summarise(swe_mean = mean(SWE), # make generic (i.e. var rather than swe)
              swe_sd = sd(SWE))

  pca_object <- calc_pcs(dat, scale = scale)
  eigenvalues <- get_eigenvalues(dat)

  eofs <- get_eofs(dat, pca_object, k, eigenvalues, rotate)

  amplitudes <- pca_object$x %>%
    .[,1:k] %>%
    scale() %>% # scale each amplitude series by its sd
    {if(rotate == TRUE) . %*% eofs$rotation_matrix else .} %>%
    as_tibble(rownames = 'year', .name_repair = ~1:k) %>%
    gather(PC, amplitude, -year) %>%
    mutate(year = as.numeric(year))#%>%
    #group_nest(PC, .key = 'amplitudes')

  list(eofs = eofs$eofs,
       amplitudes = amplitudes,
       climatology = climatology,
       pca_object = pca_object,
       eigenvalues = eigenvalues)
}

get_eofs <- function(dat, pca_object, k, eigenvalues, rotate) {
  eofs <- pca_object %>%
    broom::tidy(matrix = 'variables') %>%
    filter(PC <= k) %>%
    left_join(eigenvalues[1:2], by = 'PC') %>%
    mutate(weight = value * std.dev,
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
    dplyr::select(-column) #%>%
    #group_nest(EOF, .key = 'patterns')

  list(eofs = eofs,
       rotation_matrix = rotation_matrix)
}


