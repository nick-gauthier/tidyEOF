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
plot_eofs <- function(patterns, palette){
  max_weight <- max(abs(patterns$eofs$weight))
  ggplot(patterns$eofs) +
    geom_raster(aes(x, y, fill = weight)) +
    geom_sf(data = states_wus, fill = NA, color = 'black') +
    facet_wrap(~paste0('EOF', EOF)) +
    scale_fill_scico(palette = palette, direction = 1, limits = c(-1, 1) * max_weight) +
    theme_void()+
    ggtitle('Observed March SWE EOFS')
}

plot_amps <- function(patterns) {
  patterns$amplitudes %>%
    mutate(PC = paste0('PC', PC)) %>%
    ggplot(aes(year, amplitude)) +
    geom_line() +
    facet_wrap(~PC) +
    theme_bw()
}

plot_scree <- function(patterns, k, kmax = 20){
  patterns$eigenvalues %>%
    mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated))) %>%
    filter(PC <= kmax) %>%
    ggplot(aes(x = PC, y = percent * 100)) +
    geom_linerange(aes(x = PC, ymin = low, ymax = hi)) +
    geom_point(size = 2, aes(color = multiplet)) +
    geom_text(aes(x = PC, y = cumvar_line, label = paste0(round(cumulative * 100, 0), '%')), size = 2.5, vjust = 0) +
    labs(x = "Principal Component", y = "Normalized Eigenvalue") +
    geom_vline(xintercept = k + .5, linetype = 2, color = 'red', alpha = .7) +
    theme_bw() +
    guides(color = F) +
    scale_x_continuous(breaks = seq(0, 25, 5))
}

plot_climatology <- function(patterns) {
  a <- ggplot(patterns$climatology) +
    geom_raster(aes(x, y, fill = swe_mean)) +
    geom_sf(data = states_wus, fill = NA, color = 'black') +
    theme_void() +
    scale_fill_scico(palette = 'davos', direction = -1) +
    ggtitle('Mean March SWE, 1982-2017', 'PRISM/SNOTEL Observations')

  b <- ggplot(patterns$climatology) +
    geom_raster(aes(x, y, fill = swe_sd)) +
    geom_sf(data = states_wus, fill = NA, color = 'black') +
    theme_void() +
    scale_fill_scico(palette = 'davos', direction = -1)+
    ggtitle('Standard Deviation of March SWE, 1982-2017', 'PRISM/SNOTEL Observations')

  a + b
}