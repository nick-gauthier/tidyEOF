#' Create a patterns object with validation
#'
#' Low-level constructor for patterns objects. Validates that all required
#' components are present and have the correct types.
#'
#' @param eofs A stars object containing EOF spatial patterns
#' @param amplitudes A data.frame/tibble with time column and PC amplitudes
#' @param eigenvalues A data.frame/tibble with eigenvalue statistics
#' @param total_variance Total variance of the data (sum of all eigenvalues).
#'   Needed for variance summaries and Rule N when the stored eigenvalue
#'   table is truncated (e.g., IRLBA)
#' @param k Number of components retained
#' @param proj_matrix Projection matrix for new data
#' @param rotation Rotation matrix if varimax was applied (or NULL)
#' @param climatology List with `mean` and `sd` stars objects (from get_climatology)
#' @param units List of original data units
#' @param names Variable names
#' @param scaled Logical, whether data was scaled
#' @param monthly Logical, whether monthly climatology was used
#' @param rotate Logical, whether rotation was applied
#' @param weight Logical, whether area weighting was applied
#' @param valid_pixels Indices of valid (non-NA) pixels
#' @param block_map Named list mapping each variable to its column range in the
#'   concatenated space-time matrix
#'
#' @return A patterns object
#' @keywords internal
new_patterns <- function(eofs, amplitudes, eigenvalues, k, proj_matrix,
                         rotation = NULL, climatology = NULL, units = NULL,
                         names = NULL, scaled = FALSE, monthly = FALSE,
                         rotate = FALSE, weight = TRUE, valid_pixels = NULL,
                         total_variance = NULL, block_map = NULL) {
  stopifnot(inherits(eofs, "stars"))
  stopifnot(is.data.frame(amplitudes))
  stopifnot(is.data.frame(eigenvalues))
  stopifnot(is.numeric(k), length(k) == 1)
  stopifnot(is.matrix(proj_matrix))

  structure(
    list(
      eofs = eofs,
      amplitudes = amplitudes,
      climatology = climatology,
      eigenvalues = eigenvalues,
      total_variance = total_variance,
      rotation = rotation,
      units = units,
      names = names,
      scaled = scaled,
      monthly = monthly,
      rotate = rotate,
      k = k,
      weight = weight,
      valid_pixels = valid_pixels,
      block_map = block_map,
      proj_matrix = proj_matrix
    ),
    class = "patterns"
  )
}

#' Print method for patterns objects
#'
#' @param x A patterns object
#' @param ... Additional arguments passed to print
#' @export
print.patterns <- function(x, ...) {
  cli::cli_h1("Patterns Object")
  cli::cli_text("Modes: {.field {x$k}}")
  cli::cli_text("Variables: {.field {x$names}}")
  times <- x$amplitudes$time
  cli::cli_text("Time steps: {.field {length(times)}} ({times[1]} to {times[length(times)]})")

  cli::cli_h2("Processing Options")
  cli::cli_text("Scale: {.field {x$scaled}}")
  cli::cli_text("Monthly: {.field {x$monthly}}")
  cli::cli_text("Rotated: {.field {x$rotate}}")
  cli::cli_text("Area weighted: {.field {x$weight}}")

  cli::cli_h2("Eigenvalues (% variance)")
  top_eigs <- head(x$eigenvalues, x$k)
  variance_pct <- round(top_eigs$percent, 1)
  cli::cli_text("{.val {variance_pct}}")

  invisible(x)
}

#' Extract amplitude matrix from a patterns object or data frame
#'
#' Provides a unified way to get the numeric amplitude matrix (without the time
#' column) from either a patterns object or a bare amplitudes data frame. This
#' is the single dispatch point used throughout the package.
#'
#' @param x A patterns object or a data.frame/tibble with a time column
#' @param ... Additional arguments (currently unused)
#' @return A numeric matrix with rows = time steps, columns = PCs
#' @keywords internal
#' @export
as.matrix.patterns <- function(x, ...) {
  x$amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()
}

#' Extract time vector from a patterns object or data frame
#'
#' @param x A patterns object or a data.frame/tibble with a time column
#' @return A vector of time values
#' @keywords internal
get_times <- function(x) {
  if (inherits(x, "patterns")) {
    x$amplitudes$time
  } else {
    x$time
  }
}

#' Extract amplitude matrix from a patterns object or data frame
#'
#' Like [as.matrix.patterns()] but also works on bare data frames. Optionally
#' filters to specified times before extracting.
#'
#' @param x A patterns object or a data.frame/tibble with a time column
#' @param times Optional time vector to filter to before extraction
#' @return A numeric matrix
#' @keywords internal
extract_amplitudes_matrix <- function(x, times = NULL) {
  amps <- if (inherits(x, "patterns")) x$amplitudes else x

  if (!is.null(times)) {
    times_df <- tibble::tibble(time = times)
    # Sort so predictor and response matrices filtered to the same times are
    # row-aligned regardless of each object's storage order
    amps <- dplyr::semi_join(amps, times_df, by = "time") %>%
      dplyr::arrange(time)
  }

  amps %>%
    dplyr::select(-time) %>%
    as.matrix()
}

#' Extract subset of PCs from a patterns object
#'
#' @param x A patterns object
#' @param i Index or names of PCs to keep
#' @param ... Additional arguments (unused)
#' @return A patterns object with only selected PCs
#' @export
`[.patterns` <- function(x, i, ...) {
  # Handle different types of indexing
  if(is.character(i)) {
    pc_names <- setdiff(names(x$amplitudes), "time")
    i <- match(i, pc_names)
    if(any(is.na(i))) rlang::abort(glue("Invalid PC names. Available PCs: {glue_collapse(pc_names, sep = ', ')}"), class = "tidyeof_invalid_subset")
  }

  if(any(i > x$k)) rlang::abort(glue("Index out of bounds. Requested: {glue_collapse(i[i > x$k], sep = ', ')}, but only {x$k} components available."), class = "tidyeof_invalid_subset")

  # Only leading contiguous truncation is well-defined: the eigenvalue table
  # and the projection matrix have no meaning for non-contiguous subsets
  i <- as.integer(i)
  if(!identical(i, seq_len(length(i)))) {
    rlang::abort(
      "patterns objects can only be truncated to leading modes, i.e. x[1:k].",
      class = "tidyeof_invalid_subset"
    )
  }

  # Subset components while preserving scaling
  x$eofs <- dplyr::slice(x$eofs, i, along = "PC", drop = FALSE)
  x$amplitudes <- x$amplitudes[, c(1, i + 1)] # +1 because time is first column
  x$k <- length(i)

  if(!is.null(x$rotation)) {
    x$rotation <- x$rotation[, i, drop = FALSE]
  }

  if (!is.null(x$proj_matrix)) {
    x$proj_matrix <- x$proj_matrix[, i, drop = FALSE]
  }

  # Eigenvalues stay the same but we update k
  x
}

