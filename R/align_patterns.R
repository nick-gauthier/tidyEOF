
#' Title
#'
#' @param patterns
#'
#' @return
#' @export
#'
#' @examples
flip_patterns <- function(patterns) {
  # this function tries to make sure all the patterns have roughly the same major sign, which makes plotting easier
  # alternatively could try to align to original time series variable as well
  sums <- split(patterns$eofs) %>%
    as_tibble() %>%
    dplyr::summarise(across(starts_with('PC'), ~sign(sum(.x, na.rm = TRUE)))) %>%
    unlist()

    patterns$eofs <-  sweep(patterns$eofs, MARGIN = 3, STATS = sums, FUN = "*")

    patterns$amplitudes <- patterns$amplitudes %>%
      select(-time) %>%
      sweep(MARGIN = 2, STATS = sums, FUN = '*') %>%
      bind_cols(time = patterns$amplitudes$time, .)

    return(patterns)
}

# congruence <- function(x, y) {
#   # could check that both have the same dimensions
#   t1 <- as_tibble(x$eofs) %>%
#     tidyr::pivot_wider(names_from = PC, values_from = weight) %>%
#     dplyr::select(-x, -y) %>%
#     remove_missing() # not ideal but . . .
#
#   t2 <- as_tibble(y$eofs) %>%
#     pivot_wider(names_from = PC, values_from = weight) %>%
#     select(-x, -y) %>%
#     remove_missing()
#
#   psych::factor.congruence(t1, t2)
# }
#
# align <- function(x, y) {
#   t1 <- as_tibble(x$eofs) %>%
#     pivot_wider(names_from = PC, values_from = weight) %>%
#     select(-x, -y) %>%
#     remove_missing()
#
#   t2 <- as_tibble(y$eofs) %>%
#     pivot_wider(names_from = PC, values_from = weight) %>%
#     select(-x, -y) %>%
#     remove_missing()
#
#   vegan::procrustes(as.matrix(t1), as.matrix(t2), scale = FALSE)
#
#   #psych::factor.congruence(t1, t2)
# }

# congruence(ccsm_patterns_mon, era_patterns_mon) %>%
#   as_tibble(rownames = 'ref') %>%
#   pivot_longer(-ref) %>%
#   group_by(ref) %>%
#   arrange(-abs(value), .by_group = TRUE) %>%
#   mutate(fit = case_when(abs(value) >= 0.98 ~ 'Excellent',
#                          abs(value) >= 0.92 ~ 'Good',
#                          abs(value) >= 0.82 ~ 'Borderline',
#                          abs(value) >= 0.68 ~ 'Poor',
#                          TRUE ~ 'Terrible')) %>%
#   filter(abs(value) == max(abs(value)))
#
# align(ccsm_patterns_mon, mh_mon_patterns)$rotation
# congruence(ccsm_patterns_mon, mh_mon_patterns) %>% corrplot::corrplot()