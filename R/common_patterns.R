#' Compute Common EOF Patterns from Multiple Sources
#'
#' Performs joint PCA on multiple datasets that share a spatial domain. Each
#' dataset is anomalized with its own climatology, then concatenated along time
#' for joint PCA. The resulting shared spatial patterns get source-specific
#' amplitudes and climatologies, enabling CCA coupling and cross-source prediction.
#'
#' @param datasets Named list of stars objects sharing the same spatial grid.
#'   Names become the source identifiers used for extraction.
#' @param k Number of EOF modes to retain
#' @param scale Logical, whether to scale anomalies by standard deviation (default TRUE)
#' @param rotate Logical, whether to apply varimax rotation (default FALSE)
#' @param monthly Logical, whether to use monthly climatology (default FALSE)
#' @param weight Logical, whether to apply area weighting (default TRUE)
#' @param irlba_threshold Minimum data elements to trigger IRLBA (default 500000)
#'
#' @return A `common_patterns` S3 object. Source-specific patterns are extracted
#'   with `$` or `[[` using source names (e.g., `cpat$era`). Each extracted
#'   element is a standard `patterns` object with shared EOFs but source-specific
#'   climatology and amplitudes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cpat <- common_patterns(
#'   list(era = era_coarse, phyda = phyda_coarse),
#'   k = 11, scale = TRUE
#' )
#'
#' # Extract source-specific patterns
#' cpat$era     # patterns object with ERA climatology + amplitudes
#' cpat$phyda   # patterns object with PHYDA climatology + amplitudes
#'
#' # Use with existing couple/predict workflow
#' fine_pat <- patterns(era_fine, k = 13)
#' coupled <- couple(cpat$era, fine_pat, k = 9)
#' predict(coupled, era_new)
#'
#' # Cross-source prediction
#' predict(coupled, phyda_new, predictor_patterns = cpat$phyda)
#' }
common_patterns <- function(datasets, k = 4, scale = TRUE, rotate = FALSE,
                            monthly = FALSE, weight = TRUE,
                            irlba_threshold = 500000) {

  # --- Input validation ---
  if (!is.list(datasets) || length(datasets) == 0) {
    cli::cli_abort(
      "{.arg datasets} must be a non-empty list of {.cls stars} objects.",
      class = "tidyeof_invalid_input"
    )
  }

  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    cli::cli_abort(
      "{.arg datasets} must be a named list. Each element needs a source name.",
      class = "tidyeof_invalid_input"
    )
  }

  for (nm in names(datasets)) {
    if (!inherits(datasets[[nm]], "stars")) {
      cli::cli_abort(
        "Element {.val {nm}} of {.arg datasets} must be a {.cls stars} object, not {.cls {class(datasets[[nm]])}}.",
        class = "tidyeof_invalid_input"
      )
    }
    check_single_attribute(datasets[[nm]], arg = nm)
  }

  if (isTRUE(rotate) && k <= 1) {
    cli::cli_abort(
      "Rotation requires k > 1.",
      class = "tidyeof_invalid_option"
    )
  }

  source_names <- names(datasets)
  first_dat <- datasets[[1]]

  # --- Validate spatial grids match ---
  if (length(datasets) > 1) {
    ref_dims <- stars::st_dimensions(first_dat)
    ref_spatial <- setdiff(names(ref_dims), "time")

    for (nm in source_names[-1]) {
      other_dims <- stars::st_dimensions(datasets[[nm]])
      other_spatial <- setdiff(names(other_dims), "time")

      # Check spatial dimension names
      if (!identical(ref_spatial, other_spatial)) {
        cli::cli_abort(
          "Spatial dimensions don't match: {.val {source_names[1]}} has {.val {ref_spatial}}, {.val {nm}} has {.val {other_spatial}}.",
          class = "tidyeof_grid_mismatch"
        )
      }

      # Check spatial dimension sizes and coordinates
      for (sdim in ref_spatial) {
        ref_vals <- stars::st_get_dimension_values(first_dat, sdim)
        other_vals <- stars::st_get_dimension_values(datasets[[nm]], sdim)

        if (length(ref_vals) != length(other_vals)) {
          cli::cli_abort(
            "Dimension {.field {sdim}} has {length(ref_vals)} cells in {.val {source_names[1]}} but {length(other_vals)} in {.val {nm}}.",
            class = "tidyeof_grid_mismatch"
          )
        }

        if (!isTRUE(all.equal(ref_vals, other_vals, tolerance = 1e-6))) {
          cli::cli_abort(
            "Dimension {.field {sdim}} coordinates differ between {.val {source_names[1]}} and {.val {nm}}. Datasets must share the same spatial grid.",
            class = "tidyeof_grid_mismatch"
          )
        }
      }
    }
  }

  # --- Per-source: compute climatology and anomalies ---
  source_climatologies <- purrr::map(datasets, get_climatology, monthly = monthly)
  source_anomalies <- purrr::imap(datasets, function(dat, nm) {
    get_anomalies(dat, clim = source_climatologies[[nm]], scale = scale, monthly = monthly)
  })

  # Track time steps per source for splitting amplitudes later
  source_times <- purrr::map(datasets, ~stars::st_get_dimension_values(.x, "time"))
  source_n_times <- purrr::map_int(source_times, length)

  if (length(unique(source_n_times)) > 1) {
    n_desc <- paste0(source_names, " (", source_n_times, ")", collapse = ", ")
    cli::cli_warn(
      "Datasets have different numbers of time steps: {n_desc}. Sources with more time steps will have greater influence on the shared EOFs.",
      class = "tidyeof_unequal_times"
    )
  }

  # --- Concatenate anomalies along time ---
  # Ensure consistent attribute names before concatenating
  ref_name <- names(first_dat)[[1]]
  source_anomalies <- purrr::map(source_anomalies, function(anom) {
    if (names(anom)[[1]] != ref_name) {
      setNames(anom, ref_name)
    } else {
      anom
    }
  })

  concat <- do.call(c, c(source_anomalies, along = "time"))

  # --- Compute spatial weights once from first dataset ---
  weights <- if (weight) area_weights(first_dat) else NULL

  # --- Run PCA on concatenated data (reuses existing get_eofs machinery) ---
  eofs <- get_eofs(concat, k = k, rotate = rotate,
                   irlba_threshold = irlba_threshold, weights = weights)

  # --- Compute shared signs from EOFs ---
  signs <- compute_eof_signs(eofs$spatial_patterns)

  # --- Split amplitudes by source and build per-source patterns ---
  source_patterns <- list()
  current_row <- 1L

  for (nm in source_names) {
    n <- source_n_times[[nm]]
    rows <- current_row:(current_row + n - 1L)
    current_row <- current_row + n

    # Split amplitudes and restore source-specific times
    source_amps <- eofs$amplitudes[rows, ]
    source_amps$time <- source_times[[nm]]

    # Capture source-specific units
    src_dat <- datasets[[nm]]
    src_units <- setNames(
      purrr::map(names(src_dat), ~tryCatch(units(src_dat[[.x]]), error = function(e) NULL)),
      names(src_dat)
    )

    pat <- new_patterns(
      eofs = eofs$spatial_patterns,
      amplitudes = source_amps,
      eigenvalues = eofs$eigenvalues,
      total_variance = eofs$total_variance,
      k = k,
      proj_matrix = eofs$proj_matrix,
      rotation = eofs$rotation_matrix,
      climatology = source_climatologies[[nm]],
      units = src_units,
      names = names(src_dat),
      scaled = scale,
      monthly = monthly,
      rotate = rotate,
      weight = weight,
      valid_pixels = eofs$valid_pixels,
      block_map = eofs$block_map
    )

    pat <- apply_sign_flips(pat, signs)
    source_patterns[[nm]] <- pat
  }

  # --- Build common_patterns S3 object ---
  structure(
    list(
      source_patterns = source_patterns,
      sources = source_names,
      k = k,
      shared_eigenvalues = eofs$eigenvalues,
      pattern_opts = list(
        scale = scale,
        rotate = rotate,
        monthly = monthly,
        weight = weight
      )
    ),
    class = "common_patterns"
  )
}

#' @export
`$.common_patterns` <- function(x, name) {
  sources <- .subset2(x, "sources")
  if (name %in% sources) {
    .subset2(.subset2(x, "source_patterns"), name)
  } else {
    .subset2(x, name)
  }
}

#' @export
`[[.common_patterns` <- function(x, i, ...) {
  sources <- .subset2(x, "sources")
  if (is.character(i) && length(i) == 1 && i %in% sources) {
    .subset2(.subset2(x, "source_patterns"), i)
  } else {
    .subset2(x, i)
  }
}

#' Plot method for common_patterns objects
#'
#' Shows shared EOFs on top and overlaid amplitude time series (colored by
#' source) on the bottom.
#'
#' @param x A common_patterns object
#' @param scale Amplitude scaling: "standardized" (default), "variance", or "raw"
#' @param scale_y Y-axis scaling: "fixed" (default) or "free"
#' @param overlay Optional sf object to overlay on EOF maps
#' @param overlay_color Color for overlay geometry (default "grey30")
#' @param overlay_fill Fill for overlay geometry (default NA)
#' @param ... Additional arguments (currently unused)
#' @return A patchwork object (EOFs + amplitudes)
#' @export
plot.common_patterns <- function(x,
                                 scale = c("standardized", "variance", "raw"),
                                 scale_y = c("fixed", "free"),
                                 overlay = NULL,
                                 overlay_color = "grey30",
                                 overlay_fill = NA,
                                 ...) {
  scale <- match.arg(scale)
  scale_y <- match.arg(scale_y)

  # Use the first source's patterns for shared EOFs and eigenvalues
  ref_pat <- .subset2(.subset2(x, "source_patterns"), .subset2(x, "sources")[1])
  k <- .subset2(x, "k")

  layout <- .simple_smart_layout(k)

  # --- Shared EOF plot (same as plot.patterns) ---
  p_eofs <- .plot_eofs_internal(ref_pat, layout = layout$eof,
                                overlay = overlay, overlay_color = overlay_color,
                                overlay_fill = overlay_fill)

  # --- Overlaid amplitude plot colored by source ---
  pc_names <- names(ref_pat$amplitudes)[-1]

  # Build combined amplitude data from all sources
  all_amps <- purrr::imap_dfr(.subset2(x, "source_patterns"), function(pat, nm) {
    eigs <- pat$eigenvalues %>%
      dplyr::filter(PC <= k) %>%
      dplyr::select(PC, std.dev) %>%
      dplyr::mutate(PC = pc_names[PC])

    amps <- pat$amplitudes %>%
      tidyr::pivot_longer(-time, names_to = "PC", values_to = "amplitude") %>%
      dplyr::left_join(eigs, by = "PC")

    if (scale == "variance") {
      amps <- dplyr::mutate(amps, amplitude = amplitude * std.dev)
    } else if (scale == "raw") {
      eof_scale <- split(pat$eofs) %>%
        as_tibble() %>%
        dplyr::summarise(dplyr::across(dplyr::starts_with("PC"), ~sqrt(sum(.x^2, na.rm = TRUE)))) %>%
        tidyr::pivot_longer(dplyr::everything(), names_to = "PC", values_to = "scale_factor")
      amps <- dplyr::left_join(amps, eof_scale, by = "PC") %>%
        dplyr::mutate(amplitude = amplitude * scale_factor)
    } else {
      amps <- dplyr::mutate(amps, amplitude = amplitude / std.dev)
    }

    dplyr::mutate(amps, source = nm)
  })

  facet_args <- if (!is.null(layout$pc)) layout$pc else list()

  p_amps <- ggplot2::ggplot(all_amps, ggplot2::aes(time, amplitude, color = source)) +
    ggplot2::geom_line(alpha = 0.8) +
    do.call(ggplot2::facet_wrap, c(list(~PC), facet_args)) +
    ggplot2::labs(y = "Amplitude", x = NULL, color = "Source") +
    ggplot2::theme_bw()

  if (scale_y == "fixed") {
    p_amps <- p_amps +
      ggplot2::coord_cartesian(ylim = range(all_amps$amplitude, na.rm = TRUE))
  }

  # --- Combine with patchwork ---
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork package needed for combined plots. Install with: install.packages('patchwork')\nShowing EOFs only.")
    return(p_eofs)
  }

  sources <- .subset2(x, "sources")
  patchwork::wrap_plots(p_eofs, p_amps, ncol = 1, heights = layout$heights) +
    patchwork::plot_annotation(
      title = glue::glue("Common EOF Analysis: {glue::glue_collapse(ref_pat$names, sep = ', ')}"),
      subtitle = glue::glue("k = {k} | Sources: {glue::glue_collapse(sources, sep = ', ')}")
    ) &
    ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))
}

#' @export
print.common_patterns <- function(x, ...) {
  cli::cli_h1("Common Patterns Object")
  cli::cli_text("Sources: {.val {x$sources}}")
  cli::cli_text("Shared EOF modes: {.field {x$k}}")

  cli::cli_h2("Time steps per source")
  for (nm in x$sources) {
    n_times <- nrow(x$source_patterns[[nm]]$amplitudes)
    cli::cli_text("{.field {nm}}: {n_times}")
  }

  # PC congruence: pairwise amplitude correlations on shared times
  if (length(x$sources) >= 2) {
    pairs <- utils::combn(x$sources, 2, simplify = FALSE)
    for (pair in pairs) {
      amps_a <- x$source_patterns[[pair[1]]]$amplitudes
      amps_b <- x$source_patterns[[pair[2]]]$amplitudes
      shared <- amps_a$time %in% amps_b$time
      if (sum(shared) < 3) next

      amps_a <- amps_a[shared, ]
      amps_b <- amps_b[amps_b$time %in% amps_a$time, ]
      # Align row order
      amps_a <- amps_a[order(amps_a$time), ]
      amps_b <- amps_b[order(amps_b$time), ]

      pc_cols <- setdiff(names(amps_a), "time")
      cors <- vapply(pc_cols, function(pc) {
        stats::cor(amps_a[[pc]], amps_b[[pc]], use = "pairwise.complete.obs")
      }, numeric(1))

      cli::cli_h2("PC congruence: {pair[1]} vs {pair[2]} (n = {nrow(amps_a)})")
      for (i in seq_along(pc_cols)) {
        cli::cli_text("{.field {pc_cols[i]}}: {formatC(cors[i], format = 'f', digits = 2)}")
      }
    }
  }

  cli::cli_h2("Processing Options")
  cli::cli_text("Scale: {.field {x$pattern_opts$scale}}")
  cli::cli_text("Rotate: {.field {x$pattern_opts$rotate}}")
  cli::cli_text("Monthly: {.field {x$pattern_opts$monthly}}")
  cli::cli_text("Weight: {.field {x$pattern_opts$weight}}")

  invisible(x)
}
