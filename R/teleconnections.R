#' Calculate teleconnections between a patterns object and another layer
#'
#' Computes pixel-wise correlations between a spatiotemporal field and each
#' PC amplitude time series. Checks for overlapping time steps between the
#' raster field and the PC amplitudes.
#'
#' @param dat A stars object with a time dimension
#' @param patterns A patterns object from patterns()
#' @param amplitudes Optional amplitudes tibble (defaults to patterns$amplitudes)
#'
#' @return A stars object with correlation values for each PC
#' @export
get_correlation <- function(dat, patterns, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- patterns$amplitudes

  matched <- match_times(amplitudes, dat)
  k <- ncol(matched$amps)
  dat_filtered <- dplyr::filter(dat, time %in% matched$times)

  if (k == 1) {
    # k=1: st_apply drops the PC dimension, so compute directly
    suppressWarnings(
      result <- st_apply(dat_filtered, get_spatial_dimensions(dat_filtered),
                         function(x) cor(x, matched$amps[[1]]))
    )
    return(setNames(result, "PC1"))
  }

  suppressWarnings(
    result <- st_apply(dat_filtered, get_spatial_dimensions(dat_filtered),
                       function(x) cor(x, matched$amps), .fname = 'PC') %>%
      st_set_dimensions(., 'PC', values = paste0('PC', st_get_dimension_values(., 'PC')))
  )
  # Move PC (first dim) to last: works for both 2D (PC, geometry) and 3D (PC, x, y)
  n <- length(dim(result))
  aperm(result, c(2:n, 1))
}

#' Calculate FDR-corrected significance contours for teleconnections
#'
#' Computes pixel-wise correlations with FDR correction and returns significance
#' contour lines as sf polygons.
#'
#' @param dat A stars object with a time dimension
#' @param patterns A patterns object from patterns()
#' @param fdr False discovery rate threshold (default 0.1)
#' @param amplitudes Optional amplitudes tibble (defaults to patterns$amplitudes)
#'
#' @return An sf object with significance contour polygons for each PC
#' @export
get_fdr <- function(dat, patterns, fdr = 0.1, amplitudes = NULL) {
  if(is.null(amplitudes)) amplitudes <- patterns$amplitudes

  matched <- match_times(amplitudes, dat)
  k <- ncol(matched$amps)

  if (has_geometry_dimension(dat)) {
    cli::cli_abort(
      "FDR contour maps require raster data. Use {.fn get_correlation} for geometry-based stars.",
      class = "tidyeof_geometry_unsupported"
    )
  }

  dat_filtered <- dplyr::filter(dat, time %in% matched$times)

  if (k == 1) {
    # k=1: st_apply drops the PC dimension, so compute directly.
    # FDR correction across PCs is trivial with one PC.
    suppressWarnings(
      fdr_rast <- st_apply(dat_filtered, get_spatial_dimensions(dat_filtered),
                           function(x) {
                             if (!any(is.na(x))) cor.test(x, matched$amps[[1]])$p.value else NA
                           }) %>%
        setNames('FDR')
    )
    return(
      fdr_rast %>%
        st_contour(contour_lines = TRUE, breaks = fdr) %>%
        dplyr::transmute(PC = "PC1")
    )
  }

  suppressWarnings(
    fdr_rast <- st_apply(dat_filtered, get_spatial_dimensions(dat_filtered),
                         fdr_fun, amps = matched$amps, .fname = 'PC') %>%
      aperm(c(2,3,1)) %>%
      st_apply('PC', adjust) %>%
      setNames('FDR')
  )

  fdr_rast %>%
    st_get_dimension_values('PC') %>%
    seq_along() %>%
    purrr::map(~slice(fdr_rast, 'PC', .x) %>%
          st_contour(contour_lines = TRUE, breaks = fdr) %>%
          dplyr::transmute(PC = paste0('PC', .x))) %>%
    do.call(rbind, .)
}

#' Match time steps between amplitudes and a stars object
#'
#' @param amplitudes A tibble with a time column
#' @param dat A stars object with a time dimension
#' @return A list with `times` (common time steps) and `amps` (filtered amplitude matrix)
#' @keywords internal
match_times <- function(amplitudes, dat) {
  times_amps <- amplitudes$time
  times_dat <- st_get_dimension_values(dat, 'time')
  # Use subsetting to preserve Date/POSIXct class (intersect can strip it)
  times_cor <- sort(times_amps[times_amps %in% times_dat])

  if(length(times_cor) < 2) cli::cli_abort("Need at least two time steps in common.")
  if(!(identical(times_cor, times_amps) & identical(times_cor, times_dat))) {
    cli::cli_inform("Using the time period {times_cor[1]} to {times_cor[length(times_cor)]}.")
  }

  amps <- dplyr::filter(amplitudes, time %in% times_cor) %>%
    dplyr::arrange(time) %>%
    dplyr::select(-time)

  list(times = times_cor, amps = amps)
}

fdr_fun <- function(x, amps) {
  if(!any(is.na(x))) apply(amps, 2, function(y) cor.test(x, y)$p.value) else rep(NA, ncol(amps))
}

adjust <- function(x) {
  p.adjust(x, method = 'fdr', n = sum(!is.na(x)))
}
