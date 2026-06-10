#' Compute sign vector for consistent EOF orientation
#'
#' Determines the sign (+1 or -1) needed for each EOF component so that
#' the dominant loading is positive. Works for both raster and sf geometry
#' stars objects.
#'
#' @param eofs A stars object with EOF spatial patterns (must have a PC dimension)
#' @return Named numeric vector of +1/-1 values, one per PC
#' @keywords internal
compute_eof_signs <- function(eofs) {
  # Orient so the dominant loading is positive. A zero-sum (perfectly balanced)
  # pattern maps to +1 rather than sign(0) = 0, which would zero out the mode.
  orient <- function(x) if (sum(x, na.rm = TRUE) >= 0) 1 else -1
  if (has_geometry_dimension(eofs)) {
    eof_data <- eofs[[1]]
    sums <- apply(eof_data, 2, orient)
    names(sums) <- paste0("PC", seq_along(sums))
  } else {
    sums <- split(eofs) %>%
      as_tibble() %>%
      dplyr::summarise(across(starts_with('PC'), ~orient(.x))) %>%
      unlist()
  }
  sums
}

#' Apply sign flips to a patterns object
#'
#' Flips the sign of EOFs, amplitudes, and projection matrix according to the
#' supplied sign vector. This keeps all components synchronized.
#'
#' @param patterns A patterns object
#' @param signs Named numeric vector of +1/-1 values (from compute_eof_signs)
#' @return The patterns object with signs applied
#' @keywords internal
apply_sign_flips <- function(patterns, signs) {
  patterns$eofs <- sweep(patterns$eofs, MARGIN = length(dim(patterns$eofs)), STATS = signs, FUN = "*")

  patterns$amplitudes <- patterns$amplitudes %>%
    select(-time) %>%
    sweep(MARGIN = 2, STATS = signs, FUN = '*') %>%
    bind_cols(time = patterns$amplitudes$time, .)

  if (!is.null(patterns$proj_matrix)) {
    patterns$proj_matrix <- sweep(patterns$proj_matrix, 2, signs, `*`)
  }

  patterns
}

#' Flip EOF patterns to have consistent sign
#'
#' Ensures all EOF patterns have positive-dominant loadings by flipping the sign
#' of both the spatial pattern and corresponding amplitude time series when needed.
#' This makes plotting and interpretation more consistent across analyses.
#'
#' @param patterns A patterns object from patterns()
#' @return The patterns object with signs adjusted for consistency
#' @keywords internal
flip_patterns <- function(patterns) {
  signs <- compute_eof_signs(patterns$eofs)
  apply_sign_flips(patterns, signs)
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