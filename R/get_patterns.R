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


get_patterns <- function(dat, k, mask = NULL, scale = FALSE){

  if(!is.null(mask)) dat <- semi_join(dat, mask, by = c("x", "y"))

  climatology <- dat %>%
    group_by(x, y) %>%
    summarise(swe_mean = mean(SWE),
              swe_sd = sd(SWE))

  pca_object <- calc_pcs(dat, scale = scale)
  eigs <- get_eigs(dat)

  eofs <- pca_object %>%
    broom::tidy(matrix = 'variables') %>%
    filter(PC <= k) %>%
    left_join(eigs[1:2], by = 'PC') %>%
    mutate(weight = value * std.dev,
           EOF = as.character(PC),
           column = as.character(column)) %>%
    dplyr::select(-c(std.dev, PC))

  varim <- eofs %>% # varimax rotation
    dplyr::select(-value) %>%
    pivot_wider(names_from = EOF, values_from = weight) %>%
    column_to_rownames(var = 'column') %>%
    as.matrix %>%
    varimax
  rot_mat <- varim$rotmat # rotater

  reofs <- unclass(varim$loadings) %>%
    as_tibble(rownames = 'column') %>%
    pivot_longer(-column, names_to = 'EOF', values_to = 'weight')# %>%
   # right_join(eofs, by = c('column', 'EOF'))

  eofs <- dat %>%
    spread(year, SWE) %>%
    mutate(column = as.character(1:n())) %>%
    dplyr::select(x, y, column) %>%
    full_join(reofs, by = 'column') %>%
    dplyr::select(-column) %>%
    group_nest(EOF, .key = 'patterns')

  amps <- pca_object$x %>%
    .[,1:k] %>%
    scale() %>% # scale each amplitude series by its sd
    `%*%`(rot_mat) %>%
    as_tibble(rownames = 'year', .name_repair = ~1:k) %>%
    gather(PC, amplitude, -year) %>%
    mutate(year = as.numeric(year)) %>%
    group_nest(PC, .key = 'amplitudes')

  left_join(amps, eofs, by = c('PC' = 'EOF'))
}


plot_eof <- function(patterns, palette, normalized = FALSE){
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