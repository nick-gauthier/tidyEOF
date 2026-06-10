#' Compute spatial error metrics between predicted and observed fields
#'
#' @param predicted A stars object with predicted values
#' @param observed A stars object with observed values
#' @param metrics Character vector of metrics to compute. Options: "rmse", "cor_spatial", "cor_temporal"
#'
#' @return Named list of computed metrics
#' @keywords internal
compute_spatial_metrics <- function(predicted, observed, metrics = c("rmse", "cor_spatial", "cor_temporal")) {
  # Align times before comparison (handles dropped NA rows in predictions)
  pred_times <- stars::st_get_dimension_values(predicted, "time")
  obs_times <- stars::st_get_dimension_values(observed, "time")
  common_times <- pred_times[pred_times %in% obs_times]

  if (length(common_times) == 0) {
    cli::cli_abort("No common time steps between predicted and observed fields.")
  }

  if (length(common_times) < length(pred_times) || length(common_times) < length(obs_times)) {
    cli::cli_warn("Time mismatch: predicted has {length(pred_times)}, observed has {length(obs_times)}, using {length(common_times)} common times.")
  }


  # Filter to common times and sort both to ensure alignment
  common_times <- sort(common_times)
  predicted <- dplyr::filter(predicted, time %in% common_times)
  observed <- dplyr::filter(observed, time %in% common_times)
  predicted <- dplyr::slice(predicted, "time", order(stars::st_get_dimension_values(predicted, "time")))
  observed <- dplyr::slice(observed, "time", order(stars::st_get_dimension_values(observed, "time")))

  # Flatten to matrices for comparison
  pred_flat <- flatten_time_space(predicted)
  obs_flat <- flatten_time_space(observed)

  results <- list()

  if ("rmse" %in% metrics) {
    results$rmse <- calc_rmse(pred_flat$matrix, obs_flat$matrix)
  }

  if ("cor_spatial" %in% metrics) {
    results$cor_spatial <- calc_cor_spatial(pred_flat$matrix, obs_flat$matrix)
  }

  if ("cor_temporal" %in% metrics) {
    results$cor_temporal <- calc_cor_temporal(pred_flat$matrix, obs_flat$matrix)
  }

  results
}

#' Calculate root mean squared error
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#'
#' @return Scalar RMSE value
#' @keywords internal
calc_rmse <- function(pred_matrix, obs_matrix) {
  error <- pred_matrix - obs_matrix
  sqrt(mean(error^2, na.rm = TRUE))
}

#' Calculate spatial correlation (averaged over time)
#'
#' For each time step, compute correlation across spatial locations,
#' then average across all time steps.
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#'
#' @return Mean spatial correlation across time steps
#' @keywords internal
calc_cor_spatial <- function(pred_matrix, obs_matrix) {
  n_times <- nrow(pred_matrix)
  cors <- vapply(seq_len(n_times), function(t) {
    cor(pred_matrix[t, ], obs_matrix[t, ], use = "pairwise.complete.obs")
  }, numeric(1))
  mean(cors, na.rm = TRUE)
}

#' Calculate temporal correlation (averaged over space)
#'
#' For each spatial location, compute correlation across time,
#' then average across all locations.
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#'
#' @return Mean temporal correlation across spatial locations
#' @keywords internal
calc_cor_temporal <- function(pred_matrix, obs_matrix) {
  n_space <- ncol(pred_matrix)
  cors <- vapply(seq_len(n_space), function(s) {
    cor(pred_matrix[, s], obs_matrix[, s], use = "pairwise.complete.obs")
  }, numeric(1))
  mean(cors, na.rm = TRUE)
}
