#' Scree plot for EOF patterns
#'
#' @param x A patterns object from patterns()
#' @param k Optional number of components to highlight with vertical line
#' @param kmax Maximum number of components to show (default 10)
#' @param rule_n Logical, whether to show the modified Rule N significance
#'   cutoff as a dashed blue line (default FALSE)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot object
#' @export
#' @examples
#' \dontrun{
#' pat <- patterns(data, k = 5)
#' screeplot(pat)
#' screeplot(pat, k = 3, kmax = 8)
#' screeplot(pat, rule_n = TRUE)  # show significance cutoff
#' }
screeplot.patterns <- function(x, k = NULL, kmax = 10, rule_n = FALSE, ...) {
  # Compute Rule N cutoff if requested
  rule_n_k <- if (rule_n) rule_n_cutoff(x) else NULL

  x$eigenvalues %>%
    dplyr::mutate(separated = if_else(is.na(lag(low)), TRUE, hi < lag(low)),
           multiplet = as.factor(cumsum(separated)),
           cumvar_line = hi + 0.02 * max(hi)) %>%
    filter(PC <= kmax) %>%
    ggplot2::ggplot(aes(x = PC, y = percent)) +
    ggplot2::geom_linerange(aes(x = PC, ymin = low, ymax = hi)) +
    ggplot2::geom_point(size = 2, aes(color = multiplet)) +
    ggplot2::geom_text(aes(x = PC, y = cumvar_line, label = glue::glue("{round(cumulative, 0)}%")), size = 2.5, vjust = 0) +
    ggplot2::labs(x = "Principal Component", y = "Normalized Eigenvalue") +
    { if (!is.null(k)) ggplot2::geom_vline(xintercept = k + .5, linetype = 2, color = 'red', alpha = .7) } +
    { if (!is.null(rule_n_k) && rule_n_k > 0 && rule_n_k <= kmax) ggplot2::geom_vline(xintercept = rule_n_k + .5, linetype = 2, color = 'blue', alpha = .7) } +
    { if (!is.null(rule_n_k) && rule_n_k > kmax) ggplot2::annotate("text", x = kmax, y = Inf, label = glue::glue("All shown modes significant\n(Rule N cutoff: {rule_n_k})"), hjust = 1, vjust = 1.5, size = 3, color = 'blue') } +
    ggplot2::theme_bw() +
    ggplot2::guides(color = 'none') +
    ggplot2::scale_x_continuous(breaks = seq_len(kmax)) +
    ggplot2::scale_color_brewer(palette = 'Spectral')
}

#' Plot method for patterns objects
#'
#' @param x A patterns object
#' @param type Type of plot: "combined" (default), "eofs", or "amplitudes"
#' @param scaled For EOFs: show correlations (TRUE) or raw loadings (FALSE)
#' @param rawdata Optional raw data for correlation calculation when scaled = TRUE
#' @param scale For amplitudes: scaling method ("standardized", "variance", "raw")
#' @param scale_y For amplitudes: y-axis scaling ("fixed" or "free")
#' @param events For amplitudes: optional dates to mark with vertical lines
#' @param overlay Optional sf object to overlay on EOF maps (e.g., coastlines, boundaries)
#' @param overlay_color Color for overlay geometry (default "grey30")
#' @param overlay_fill Fill for overlay geometry (default NA for no fill)
#' @param ... Additional arguments (currently unused)
#' @return A ggplot2 object or patchwork object for combined plots
#' @export
plot.patterns <- function(x,
                          type = "combined",
                          scaled = FALSE,
                          rawdata = NULL,
                          scale = c("standardized", "variance", "raw"),
                          scale_y = c("fixed", "free"),
                          events = NULL,
                          overlay = NULL,
                          overlay_color = "grey30",
                          overlay_fill = NA,
                          ...) {

  type <- match.arg(type, c("combined", "eofs", "amplitudes"))
  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  if(type == "combined") {
    # Smart combined plot: EOFs above, PCs below with KISS layout logic
    if(!requireNamespace("patchwork", quietly = TRUE)) {
      warning("patchwork package needed for combined plots. Install with: install.packages('patchwork')\nShowing EOFs only.")
      return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill))
    } else {
      # Simple smart layout - covers 80% of use cases perfectly
      layout <- .simple_smart_layout(x$k)

      # Create EOF plot with smart layout
      p_eofs <- .plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, layout = layout$eof, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill)

      # Create amplitude plot with matching layout
      p_amps <- .plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events, layout = layout$pc)

      # Combine with patchwork using top-over-bottom layout
      return(
        patchwork::wrap_plots(p_eofs, p_amps, ncol = 1, heights = c(layout$heights[1] * length(x$eofs), layout$heights[2])) +
             patchwork::plot_annotation(
               title = glue::glue("EOF Analysis: {glue::glue_collapse(x$names, sep = ', ')}"),
               subtitle = glue::glue("k = {x$k} modes | {ifelse(x$scaled, 'scaled', 'unscaled')} | {ifelse(x$rotate, 'rotated', 'unrotated')}")
             ) &
             ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5)))
    }
  } else if(type == "eofs") {
    return(.plot_eofs_internal(x, scaled = scaled, rawdata = rawdata, overlay = overlay, overlay_color = overlay_color, overlay_fill = overlay_fill))
  } else if(type == "amplitudes") {
    return(.plot_amplitudes_internal(x, scale = scale, scale_y = scale_y, events = events))
  }
}

#' Simple smart layout function - KISS principle
#' @keywords internal
.simple_smart_layout <- function(k) {
  if(k <= 3) {
    # Single row: works for most cases, clean alignment
    list(eof = list(nrow = 1), pc = list(nrow = 1), heights = c(5, 1))
  } else {
    # k >= 4: Use square-ish grid
    ncol <- ceiling(sqrt(k))
    list(eof = list(ncol = ncol), pc = list(ncol = ncol), heights = c(5, 1))
  }
}

#' Internal function for EOF plotting
#' @keywords internal
.plot_eofs_internal <- function(x, scaled = FALSE, rawdata = NULL, layout = NULL,
                                 overlay = NULL, overlay_color = "grey30", overlay_fill = NA) {
  facet_args <- if(!is.null(layout)) layout else list()

  overlay_layer <- if(!is.null(overlay)) {
    ggplot2::geom_sf(data = overlay, fill = overlay_fill, color = overlay_color, inherit.aes = FALSE)
  } else {
    NULL
  }

  if(scaled) {
    if (length(x$eofs) > 1) {
      cli::cli_abort(
        c(
          "Correlation (scaled) EOF maps are only supported for single-variable patterns.",
          "i" = "Fit patterns on one variable, or plot with {.code scaled = FALSE}."
        ),
        class = "tidyeof_multivariate_unsupported"
      )
    }
    if(is.null(rawdata)) {
      rlang::abort("rawdata must be provided when scaled = TRUE for correlation calculation", class = "tidyeof_missing_rawdata")
    }
    return(
      ggplot2::ggplot() +
        stars::geom_stars(data = get_correlation(rawdata, x)) +
        overlay_layer +
        do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
        ggplot2::scale_fill_distiller(palette = 'RdBu', na.value = NA, limits = c(-1, 1)) +
        ggplot2::coord_sf() +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "right") +
        ggplot2::labs(fill = "Correlation")
    )
  }

  vars <- names(x$eofs)

  # One panel row per variable, each with its own fill scale: loadings of
  # different variables are not comparable on a shared color scale
  plot_one <- function(v) {
    ggplot2::ggplot() +
      stars::geom_stars(data = x$eofs[v]) +
      overlay_layer +
      do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
      scico::scale_fill_scico(palette = 'vik', midpoint = 0, na.value = NA) +
      ggplot2::coord_sf() +
      ggplot2::theme_void() +
      ggplot2::theme(legend.position = "right") +
      ggplot2::labs(fill = if (length(vars) == 1) "Loading" else v)
  }

  if (length(vars) == 1) {
    return(plot_one(vars))
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package needed for multivariate EOF plots. Showing first variable only.")
    return(plot_one(vars[1]))
  }

  patchwork::wrap_plots(purrr::map(vars, plot_one), ncol = 1)
}

#' Internal function for amplitude plotting
#' @keywords internal
.plot_amplitudes_internal <- function(x, scale = "standardized", scale_y = "fixed", events = NULL, layout = NULL) {

  # Use provided layout or let facet_wrap choose defaults
  facet_args <- if(!is.null(layout)) layout else list()


  # Get PC names directly from amplitudes to ensure matching

  pc_names <- names(x$amplitudes)[-1]  # exclude 'time' column

  # Calculate scaling factors based on method
  if(scale == "variance") {
    # Multiply by sqrt(eigenvalue) to show variance contribution
    eigs <- x$eigenvalues %>%
      dplyr::filter(PC <= x$k) %>%
      dplyr::select(PC, std.dev) %>%
      dplyr::mutate(PC = pc_names[PC])

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude * std.dev)
  } else if(scale == "raw") {
    if (length(x$eofs) > 1) {
      cli::cli_abort(
        "Raw amplitude scaling mixes units across variables and is not supported for multivariate patterns. Use scale = 'standardized' or 'variance'.",
        class = "tidyeof_multivariate_unsupported"
      )
    }
    # Use EOF loadings to convert back to original units
    eigs <- split(x$eofs) %>%
      as_tibble() %>%
      dplyr::summarise(dplyr::across(dplyr::starts_with('PC'), ~sqrt(sum(.x^2, na.rm = TRUE)))) %>%
      tidyr::pivot_longer(dplyr::everything(), names_to = 'PC', values_to = 'scale_factor')

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude * scale_factor)
  } else {
    # Default: standardized (divide by sdev to get unit variance)
    eigs <- x$eigenvalues %>%
      dplyr::filter(PC <= x$k) %>%
      dplyr::select(PC, std.dev) %>%
      dplyr::mutate(PC = pc_names[PC])

    amps <- x$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = 'PC', values_to = 'amplitude') %>%
      dplyr::left_join(eigs, by = 'PC') %>%
      dplyr::mutate(amplitude = amplitude / std.dev)
  }

  p <- ggplot2::ggplot(amps, ggplot2::aes(time, amplitude)) +
    ggplot2::geom_line() +
    do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
    ggplot2::labs(y = "Amplitude", x = NULL)

  if(scale_y == "fixed") {
    p <- p + ggplot2::coord_cartesian(ylim = range(amps$amplitude, na.rm = TRUE))
  }

  if(!is.null(events)) {
    p <- p + ggplot2::geom_vline(xintercept = events, linetype = 2, alpha = 0.5)
  }

  p + ggplot2::theme_bw()
}
