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