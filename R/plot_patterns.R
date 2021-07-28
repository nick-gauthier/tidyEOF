#' Title
#'
#' @param patterns
#' @param palette
#' @param normalized
#'
#' @return
#' @export
#'
#' @examples
plot_eofs <- function(patterns, scaled = FALSE, rawdata = NULL){
  if(scaled){
    ggplot() +
      geom_stars(data = get_correlation(rawdata, patterns)) +
      facet_wrap(~PC) +
      scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1)) +
      coord_quickmap() +
      theme_void()
  } else {
    max_weight <- max(abs(patterns$eofs$weight))
    ggplot() +
      geom_stars(data = patterns$eofs) +
      facet_wrap(~ PC) +
      scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1) * max_weight) +
      coord_quickmap() +
      theme_void()
  }
}

#' @export
plot_amps <- function(patterns, scaled = TRUE, events = NULL) {
  if(!scaled) {
    eigs <- split(patterns$eofs) %>%
      as_tibble() %>%
      # get the sqrt of sum of squared loadings (works even if rotated)
      dplyr::summarise(across(starts_with('PC'), ~sqrt(sum(.x^2, na.rm = TRUE)))) %>%
      pivot_longer(everything(), names_to = 'PC', values_to = 'std.dev')

    amps <- patterns$amplitudes %>%
      pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude * std.dev)
  } else {
    amps <-  patterns$amplitudes %>%
      pivot_longer(-time, names_to = 'PC', values_to = 'amplitude')
  }

    ggplot(amps, aes(time, amplitude)) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = events, color = 'red', linetype = 2) +
    facet_wrap(~PC) +
    theme_bw()
}

#' @export
plot_scree <- function(dat, k, kmax = 10, scale = FALSE){
  get_pcs(dat, scale = scale) %>%
    get_eigenvalues() %>%
    dplyr::mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated))) %>%
    filter(PC <= kmax) %>%
    ggplot(aes(x = PC, y = percent * 100)) +
    geom_linerange(aes(x = PC, ymin = low, ymax = hi)) +
    geom_point(size = 2, aes(color = multiplet)) +
    geom_text(aes(x = PC, y = cumvar_line, label = paste0(round(cumulative * 100, 0), '%')), size = 2.5, vjust = 0) +
    labs(x = "Principal Component", y = "Normalized Eigenvalue") +
    geom_vline(xintercept = k + .5, linetype = 2, color = 'red', alpha = .7) +
    theme_bw() +
    guides(color = 'none') +
    scale_x_continuous(breaks = seq(0, 25, 5)) +
    scale_color_brewer(palette = 'Spectral')
}
#' @rdname plot_eofs
#' @export
