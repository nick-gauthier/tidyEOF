#' Project New Data onto Existing Patterns
#'
#' This function projects new spatial-temporal data onto existing EOF patterns,
#' returning the corresponding principal component time series. This is a core
#' function used in pattern-based downscaling and reconstruction.
#'
#' @param patterns A patterns object containing EOFs, climatology, and other metadata
#' @param newdata A stars object with new spatial-temporal data to project
#'
#' @return A tibble with time column and PC amplitude columns
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get patterns from training data
#' pat <- patterns(training_data, k = 5)
#'
#' # Project new data onto these patterns
#' new_amplitudes <- project_patterns(pat, new_data)
#' }
project_patterns <- function(patterns, newdata) {
  # NOTE: No spatial compatibility validation currently performed.
  # Assumes newdata has same spatial structure as original data used to create patterns.
  # Mismatched data will cause errors downstream (matrix dimension mismatches).

  validate_patterns(patterns)

  # Univariate patterns accept any single-attribute newdata (name-agnostic,
  # e.g. cross-source prediction where the variable is named differently).
  # Multivariate projection needs the same variable set, reordered to match.
  if (length(patterns$names) > 1 || length(newdata) > 1) {
    if (!setequal(names(newdata), patterns$names)) {
      cli::cli_abort(
        c(
          "Attributes of {.arg newdata} ({.field {names(newdata)}}) must match the training variables ({.field {patterns$names}}).",
          "i" = "Multivariate patterns require the same set of variables."
        ),
        class = "tidyeof_attribute_mismatch"
      )
    }
    newdata <- newdata[patterns$names]
  }

  new_times <- st_get_dimension_values(newdata, 'time')

  anomalies <- newdata %>%
    get_anomalies(patterns$climatology,
                  scale = patterns$scaled,
                  monthly = patterns$monthly)

  # Apply area weighting if specified in patterns
  if (patterns$weight) {
    anomalies <- anomalies * area_weights(newdata)
  }

  anomalies <- units::drop_units(anomalies)

  flattened <- flatten_time_space(anomalies)
  data_matrix <- flattened$matrix

  valid_pixels <- if (!is.null(patterns$valid_pixels)) {
    patterns$valid_pixels
  } else {
    seq_len(ncol(data_matrix))
  }

  # Check that valid_pixels indices are within bounds
  if (max(valid_pixels) > ncol(data_matrix)) {
    cli::cli_abort(c(
      "valid_pixels indices exceed newdata dimensions.",
      "x" = "max(valid_pixels) = {max(valid_pixels)}, but newdata has only {ncol(data_matrix)} spatial columns.",
      "i" = "This suggests train and test data have different spatial extents."
    ))
  }

  data_matrix <- data_matrix[, valid_pixels, drop = FALSE]

  complete_rows <- stats::complete.cases(data_matrix)
  if (!all(complete_rows)) {
    data_matrix <- data_matrix[complete_rows, , drop = FALSE]
    times_complete <- new_times[complete_rows]
  } else {
    times_complete <- new_times
  }

  if (nrow(data_matrix) == 0) {
    cli::cli_abort("No complete time steps available after removing missing values in `newdata`.")
  }

  # Require proj_matrix (created by get_patterns)
 if (is.null(patterns$proj_matrix)) {
    cli::cli_abort(c(
      "Patterns object is missing {.field proj_matrix}.",
      "i" = "This patterns object may have been created with an older version of tidyeof.",
      "i" = "Please regenerate it using {.fn patterns}."
    ))
  }

  # Check dimension compatibility
  expected_cols <- nrow(patterns$proj_matrix)
  if (ncol(data_matrix) != expected_cols) {
    cli::cli_abort(c(
      "Spatial dimension mismatch between patterns and newdata.",
      "x" = "Patterns expect {expected_cols} spatial pixels, but newdata has {ncol(data_matrix)}.",
      "i" = "This can happen if newdata has a different spatial extent or NA pattern."
    ))
  }

  # Direct projection using pre-composed matrix
  # proj_matrix has shape (n_pixels x n_components), already includes rotation + flips
  # data_matrix has shape (n_times x n_pixels)
  # result has shape (n_times x n_components)
  projected <- data_matrix %*% patterns$proj_matrix

  # Convert to tibble with proper names
  result <- projected %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(names0(patterns$k, 'PC')) %>%
    mutate(time = times_complete, .before = 1)

  return(result)
}

#' Project Multiple Datasets onto Patterns
#'
#' Convenience function to project multiple datasets onto the same patterns
#'
#' @param patterns A patterns object
#' @param data_list List of stars objects to project
#' @param names Optional names for the datasets
#'
#' @return List of projected amplitude tibbles
#'
#' @export
project_patterns_multiple <- function(patterns, data_list, names = NULL) {

  if (!is.list(data_list)) {
    cli::cli_abort("`data_list` must be a list of {.cls stars} objects.")
  }

  # Apply projection to each dataset
  results <- purrr::map(data_list, ~ project_patterns(patterns, .x))

  # Add names if provided
  if (!is.null(names)) {
    if (length(names) != length(results)) {
      cli::cli_abort("Length of `names` must match length of `data_list`.")
    }
    names(results) <- names
  }

  return(results)
}
