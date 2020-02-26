
#' Get eigenvalues from pc object
#'
#' @param pc_object
#' @param raster
#'
#' @return
#' @export
#'
#' @examples
get_eigs <- function(pc_object, raster){
  pc_object %>%
    broom::tidy(matrix = 'pcs') %>%
    mutate(eigenvalues = std.dev ^ 2,
           error = sqrt(2 / n_effective(raster)),
           low =  eigenvalues * (1 - error) * 100 / sum(eigenvalues),
           hi = eigenvalues * (1 + error) * 100 / sum(eigenvalues),
           cumvar_line = hi + 0.02 * max(hi))
}

n_effective <- function(x){
  n <- nlayers(x)
  x %>%
    area_weight %>%
    as.data.frame(na.rm = TRUE) %>%
    filter_all(any_vars(floor(.) != 0)) %>%
    t %>%
    as_tibble %>%
    gather(cell, value) %>%
    nest(data = c(value)) %>%
    mutate(rho = map_dbl(data, ~cor(.$value, lag(.$value), use = 'comp'))) %>%
    remove_missing() %>%
    mutate(effective_n = n * (1 - rho^2) / (1 + rho^2)) %>% # from Bretherton et al 1999
    summarise(mean(effective_n)) %>%
    pull
}

plot_scree <- function(eigs, k){
  eigs %>%
    mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated))) %>%
    filter(PC <= 25) %>%
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

