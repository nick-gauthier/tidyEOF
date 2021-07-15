#' Calculate teleconnections between a patterns object and another layer
#'
#' @param dat
#' @param patterns
#'
#' @return
#' @export
#'
#' @examples
get_corrs <- function(dat, patterns) { # should be able to take arbitrary number of amps rather than first 4
  amps <- patterns$amplitudes %>%
    spread(PC, amplitude, sep ='')

  dat %>%
    group_by(x, y) %>%
    nest %>%
    mutate(corrs = map(data,  ~inner_join(., amps, by = 'year') %>%
                         summarise(PC1 = cor(value, PC1),
                                   PC2 = cor(value, PC2),
                                   PC3 = cor(value, PC3),
                                   PC4 = cor(value, PC4)))) %>%
    select(-data) %>%
    unnest(corrs)
}

get_fdr <- function(dat, patterns, fdr = 0.1) { # could combine with above
  amps <- patterns$amplitudes %>%
    spread(PC, amplitude, sep ='')

  dat %>%
    group_by(x,y) %>%
    nest %>%
    mutate(corrs = map(data,  ~inner_join(., amps, by = 'year') %>%
                         summarise(PC1 = cor.test(value, PC1)$p.value,
                                   PC2 = cor.test(value, PC2)$p.value,
                                   PC3 = cor.test(value, PC3)$p.value,
                                   PC4 = cor.test(value, PC4)$p.value))) %>%
    select(-data) %>%
    unnest(corrs) %>%
    gather(PC, value, PC1:PC4) %>%
    group_by(PC) %>%
    mutate(value = p.adjust(value, method = 'fdr')) %>%
    filter(value < fdr) %>%
    ungroup()
}