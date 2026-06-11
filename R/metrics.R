#' Compute spatial error metrics between predicted and observed fields
#'
#' @param predicted A stars object with predicted values
#' @param observed A stars object with observed values
#' @param metrics Character vector of metrics to compute. Options: "rmse", "cor_spatial", "cor_temporal"
#'
#' @return Named list of computed metrics
#' @keywords internal
compute_spatial_metrics <- function(predicted, observed, metrics = c("rmse", "cor_spatial", "cor_temporal")) {
  if (length(predicted) > 1 || length(observed) > 1) {
    if (!setequal(names(predicted), names(observed))) {
      cli::cli_abort(
        "Predicted attributes ({.field {names(predicted)}}) must match observed attributes ({.field {names(observed)}}).",
        class = "tidyeof_attribute_mismatch"
      )
    }
    observed <- observed[names(predicted)]
  }

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

  compute_block_metrics(pred_flat$matrix, obs_flat$matrix,
                        pred_flat$block_map, metrics)
}

#' Compute per-variable and pooled metrics over block-structured matrices
#'
#' The plain metric name is always the pooled score. For a single block it is
#' the metric itself, preserving univariate behavior exactly. For multiple
#' blocks, RMSE is normalized by each block's observed standard deviation and
#' RMS-combined (raw pooling across different units is meaningless), and
#' correlations (already unitless) are averaged. Per-variable values get
#' suffixed names (e.g. rmse_tmean) only when there is more than one block.
#'
#' @details Pooled RMSE assumes each block's observed standard deviation is
#'   > 0; a constant (zero-variance) observed variable yields a non-finite
#'   pooled score.
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#' @param block_map Named list of column indices per variable
#' @param metrics Character vector of metric names
#' @return Named list of metric values
#' @keywords internal
compute_block_metrics <- function(pred_matrix, obs_matrix, block_map,
                                  metrics = c("rmse", "cor_spatial", "cor_temporal")) {
  calc <- list(rmse = calc_rmse, cor_spatial = calc_cor_spatial,
               cor_temporal = calc_cor_temporal)
  metrics <- intersect(metrics, names(calc))

  per_var <- purrr::map(block_map, function(cols) {
    p <- pred_matrix[, cols, drop = FALSE]
    o <- obs_matrix[, cols, drop = FALSE]
    vals <- purrr::map(calc[metrics], ~.x(p, o))
    vals$.obs_sd <- stats::sd(o, na.rm = TRUE)
    vals
  })

  results <- list()
  for (m in metrics) {
    vals <- purrr::map_dbl(per_var, m)
    if (length(per_var) == 1) {
      results[[m]] <- vals[[1]]
    } else {
      if (m == "rmse") {
        nrmse <- vals / purrr::map_dbl(per_var, ".obs_sd")
        results[[m]] <- sqrt(mean(nrmse^2))
      } else {
        results[[m]] <- mean(vals)
      }
      for (v in names(per_var)) {
        results[[paste0(m, "_", v)]] <- per_var[[v]][[m]]
      }
    }
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

#' Calculate spatial anomaly correlation (averaged over time)
#'
#' For each time step, compute correlation across spatial locations,
#' then average across all time steps. Each cell's temporal mean is removed
#' first so the correlation measures agreement of the anomaly patterns rather
#' than the shared climatology, which is constant in time and would otherwise
#' inflate the correlation toward one regardless of skill.
#'
#' @param pred_matrix Predicted values matrix (time x space)
#' @param obs_matrix Observed values matrix (time x space)
#'
#' @return Mean spatial anomaly correlation across time steps
#' @keywords internal
calc_cor_spatial <- function(pred_matrix, obs_matrix) {
  pred_anom <- sweep(pred_matrix, 2, colMeans(pred_matrix, na.rm = TRUE))
  obs_anom <- sweep(obs_matrix, 2, colMeans(obs_matrix, na.rm = TRUE))
  n_times <- nrow(pred_anom)
  cors <- vapply(seq_len(n_times), function(t) {
    cor(pred_anom[t, ], obs_anom[t, ], use = "pairwise.complete.obs")
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
