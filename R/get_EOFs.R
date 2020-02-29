#' Get EOFs from the PC object
#'
#' @param dat
#' @param pca_object
#' @param eigenvalues
#' @param k
#'
#' @return
#' @export
#'
#' @examples
#'

get_patterns <- function(pc_object, data, eigs, k, mask = TRUE) {
amps <- get_amps(pc_object, k) %>%
  group_nest(PC, .key = 'amplitudes')

eofs <- get_eofs(data, pc_object, eigs, k) %>%
  {if(mask) semi_join(., states_mask) else .} %>%
  group_nest(EOF, .key = 'patterns')

left_join(amps, eofs, by = c('PC' = 'EOF'))
}

get_eofs <- function(dat, pca_object, eigenvalues, k){
  eofs <- pca_object %>%
    tidy(matrix = 'variables') %>%
    filter(PC <= k) %>%
    left_join(eigenvalues[1:2], by = 'PC') %>%
    mutate(weight = value * std.dev,
           EOF = as.character(PC)) %>%
    dplyr::select(-std.dev, -PC)

  dat %>%
    spread(year, SWE) %>%
    mutate(column = 1:n()) %>%
    dplyr::select(x, y, column) %>%
    full_join(eofs, by = 'column') %>%
    dplyr::select(-column)
}

plot_eof <- function(patterns, palette){
  eofs <- patterns %>%
    dplyr::select(-amplitudes) %>%
    unnest(patterns) %>%
    mutate(EOF = paste0('EOF', PC))

  ggplot(eofs) +
    geom_raster(aes(x, y, fill = weight)) +
    geom_sf(data = states_wus, fill = NA, color = 'black') +
    facet_wrap(~EOF) +
    scale_fill_scico(palette = palette, direction = 1, limits = c(-1, 1) * max(abs(eofs$weight))) +
    theme_void()+
    ggtitle('Observed March SWE EOFS')
}