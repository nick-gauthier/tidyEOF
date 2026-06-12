#' Prepare cross-validation folds for CCA tuning
#'
#' Computes EOF patterns at maximum k once per fold, enabling cheap truncation
#' during grid search. This is the expensive step - call it once, then use
#' `tune_cca()` for fast hyperparameter exploration.
#'
#' @param predictor A stars object with predictor data
#' @param response A stars object with response data
#' @param kfolds Number of cross-validation folds (default 5)
#' @param max_k_pred Maximum number of predictor EOFs to compute (default 10)
#' @param max_k_resp Maximum number of response EOFs to compute (default 10)
#' @param scale Logical, whether to scale data before EOF extraction (default FALSE).
#'   Sets the default for both predictor and response; use `scale_pred` and/or
#'   `scale_resp` to override individually.
#' @param scale_pred Logical, whether to scale predictor data before EOF extraction.
#'   Overrides `scale` for the predictor when not NULL (default NULL).
#' @param scale_resp Logical, whether to scale response data before EOF extraction.
#'   Overrides `scale` for the response when not NULL (default NULL).
#' @param rotate Logical, whether to apply varimax rotation (default FALSE)
#' @param monthly Logical, whether to compute monthly climatology (default FALSE)
#' @param weight Logical, whether to apply area weighting (default TRUE)
#' @param common_with Optional named list of additional stars objects to include in
#'   common EOF computation via [common_patterns()]. When provided, predictor patterns
#'   are computed jointly with these datasets. The primary predictor is included under
#'   the name `.primary`. Test times are excluded from all datasets to prevent leakage.
#'
#' @return A cv_folds S3 object containing:
#'   \item{folds}{List of fold data, each containing train patterns and test data}
#'   \item{max_k_pred}{Maximum predictor k used}
#'   \item{max_k_resp}{Maximum response k used}
#'   \item{kfolds}{Number of folds}
#'   \item{common_times}{Vector of overlapping time steps}
#'   \item{pattern_opts}{List of pattern extraction options}
#'   \item{common_with_sources}{Names of common EOF sources (if used)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Standard folds
#' cv <- prep_cv_folds(coarse_data, fine_data,
#'                     kfolds = 5, max_k_pred = 10, max_k_resp = 10)
#'
#' # With common EOFs from additional sources
#' cv <- prep_cv_folds(coarse_data, fine_data,
#'                     common_with = list(phyda = phyda_coarse),
#'                     kfolds = 5, max_k_pred = 10, max_k_resp = 10)
#'
#' # Then tune over k_pred and k_resp (unchanged)
#' results <- tune_cca(cv, k_pred = 1:10, k_resp = 1:10)
#' }
prep_cv_folds <- function(predictor, response,
                          kfolds = 5,
                          max_k_pred = 10,
                          max_k_resp = 10,
                          scale = FALSE,
                          scale_pred = NULL,
                          scale_resp = NULL,
                          rotate = FALSE,
                          monthly = FALSE,
                          weight = TRUE,
                          common_with = NULL) {
  if (isTRUE(rotate)) {
    cli::cli_abort(
      "rotate = TRUE is not supported in cross-validation. Rotation depends on k, so truncating max-k patterns invalidates the rotated basis.",
      class = "tidyeof_invalid_option"
    )
  }

  # Resolve per-field scale settings
  scale_pred <- scale_pred %||% scale
  scale_resp <- scale_resp %||% scale

  # Validate common_with if provided
  if (!is.null(common_with)) {
    if (!is.list(common_with) || is.null(names(common_with)) || any(names(common_with) == "")) {
      cli::cli_abort(
        "{.arg common_with} must be a named list of {.cls stars} objects.",
        class = "tidyeof_invalid_input"
      )
    }
    if (".primary" %in% names(common_with)) {
      cli::cli_abort(
        "The name {.val .primary} is reserved. Please use a different name for your {.arg common_with} sources.",
        class = "tidyeof_invalid_input"
      )
    }
  }

  # Find common times between predictor and response
  pred_times <- stars::st_get_dimension_values(predictor, "time")
  resp_times <- stars::st_get_dimension_values(response, "time")
  common_times <- pred_times[pred_times %in% resp_times]

  if (length(common_times) == 0) {
    cli::cli_abort("No common time steps found between predictor and response.",
                   class = "tidyeof_no_common_times")
  }

  if (length(common_times) < kfolds * 2) {
    cli::cli_abort("Not enough common time steps ({length(common_times)}) for {kfolds}-fold CV.",
                   class = "tidyeof_insufficient_data")
  }

  # Filter to common times
  predictor <- dplyr::filter(predictor, time %in% common_times)
  response <- dplyr::filter(response, time %in% common_times)

  # Create fold assignments
  fold_times <- prep_folds(common_times, kfolds = kfolds)

  cli::cli_progress_step("Computing patterns for {kfolds} folds{if (!is.null(common_with)) ' (with common EOFs)' else ''}")

  # Build each fold
  folds <- purrr::imap(fold_times, function(test_times, fold_id) {
    # Training data: all times NOT in this fold
    train_pred <- dplyr::filter(predictor, !(time %in% test_times))
    train_resp <- dplyr::filter(response, !(time %in% test_times))

    # Test data: times in this fold
    test_pred <- dplyr::filter(predictor, time %in% test_times)
    test_resp <- dplyr::filter(response, time %in% test_times)

    # Compute predictor patterns - with or without common EOFs
    if (!is.null(common_with)) {
      # Filter common_with datasets to exclude test times (prevent leakage)
      common_train <- purrr::map(common_with, function(cw) {
        cw_times <- stars::st_get_dimension_values(cw, "time")
        if (any(cw_times %in% test_times)) {
          dplyr::filter(cw, !(time %in% test_times))
        } else {
          cw
        }
      })

      cpat <- common_patterns(
        c(list(.primary = train_pred), common_train),
        k = max_k_pred, scale = scale_pred, rotate = rotate,
        monthly = monthly, weight = weight
      )
      train_pred_patterns <- cpat$.primary
    } else {
      train_pred_patterns <- patterns(train_pred, k = max_k_pred,
                                          scale = scale_pred, rotate = rotate,
                                          monthly = monthly, weight = weight)
    }

    train_resp_patterns <- patterns(train_resp, k = max_k_resp,
                                        scale = scale_resp, rotate = rotate,
                                        monthly = monthly, weight = weight)

    list(
      train_pred_patterns = train_pred_patterns,
      train_resp_patterns = train_resp_patterns,
      test_pred_data = test_pred,
      test_resp_data = test_resp,
      fold_id = as.integer(fold_id)
    )
  })

  cli::cli_progress_done()

  result <- list(
    folds = folds,
    max_k_pred = max_k_pred,
    max_k_resp = max_k_resp,
    kfolds = kfolds,
    common_times = common_times,
    pattern_opts = list(
      scale_pred = scale_pred,
      scale_resp = scale_resp,
      rotate = rotate,
      monthly = monthly,
      weight = weight
    ),
    common_with_sources = if (!is.null(common_with)) names(common_with) else NULL
  )

  class(result) <- "cv_folds"
  result
}

#' Tune CCA hyperparameters via cross-validation
#'
#' Performs grid search over k_pred, k_resp, and optionally k_cca using
#' precomputed patterns from `prep_cv_folds()`. Pattern truncation is cheap,
#' so this runs quickly even with large grids.
#'
#' @param cv_folds A cv_folds object from `prep_cv_folds()`
#' @param k_pred Vector of predictor EOF counts to try (default 1:10)
#' @param k_resp Vector of response EOF counts to try (default 1:10)
#' @param k_cca Vector of CCA mode counts to try, or NULL (default) to use
#'   `min(k_pred, k_resp)` for each combination. Using fewer CCA modes than
#'   the maximum can act as regularization.
#' @param metrics Character vector of metrics to compute. Options:
#'   "rmse", "cor_spatial", "cor_temporal" (default: all three). For
#'   multivariate fields each metric also gets per-variable columns (e.g.
#'   `rmse_tmean`); the plain name is the pooled score (sd-normalized RMS for
#'   `rmse`, mean for correlations).
#' @param parallel Logical, whether to use furrr for parallel execution (default FALSE)
#'
#' @return A tibble with columns: k_pred, k_resp, k_cca, fold, and one column
#'   per metric requested.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cv <- prep_cv_folds(coarse_data, fine_data, kfolds = 5,
#'                     max_k_pred = 10, max_k_resp = 10)
#'
#' # Grid search over k_pred and k_resp (k_cca = min automatically)
#' results <- tune_cca(cv, k_pred = 1:10, k_resp = 1:10)
#'
#' # Also tune k_cca for regularization
#' results <- tune_cca(cv, k_pred = 1:10, k_resp = 1:10, k_cca = 1:5)
#'
#' # Summarize and find best params
#' summary <- summarize_cv(results, metric = "rmse")
#' }
tune_cca <- function(cv_folds,
                     k_pred = 1:10,
                     k_resp = 1:10,
                     k_cca = NULL,
                     method = "cca",
                     metrics = c("rmse", "cor_spatial", "cor_temporal"),
                     parallel = FALSE) {

  if (!inherits(cv_folds, "cv_folds")) {
    cli::cli_abort("cv_folds must be a cv_folds object from prep_cv_folds()",
                   class = "tidyeof_invalid_input")
  }
  if (isTRUE(cv_folds$pattern_opts$rotate)) {
    cli::cli_abort(
      "cv_folds were created with rotate = TRUE, which is not supported for cross-validation truncation.",
      class = "tidyeof_invalid_option"
    )
  }

  # Validate k values against max_k
  if (max(k_pred) > cv_folds$max_k_pred) {
    cli::cli_abort("max(k_pred) = {max(k_pred)} exceeds max_k_pred = {cv_folds$max_k_pred} from prep_cv_folds()",
                   class = "tidyeof_invalid_k")
  }
  if (max(k_resp) > cv_folds$max_k_resp) {
    cli::cli_abort("max(k_resp) = {max(k_resp)} exceeds max_k_resp = {cv_folds$max_k_resp} from prep_cv_folds()",
                   class = "tidyeof_invalid_k")
  }

  # Build parameter grid

  if (is.null(k_cca)) {
    # Default: k_cca = min(k_pred, k_resp) for each combination
    param_grid <- tidyr::expand_grid(
      k_pred = k_pred,
      k_resp = k_resp
    ) %>%
      dplyr::mutate(k_cca = pmin(k_pred, k_resp))
  } else {
    # Tune k_cca independently (useful for regularization)
    param_grid <- tidyr::expand_grid(
      k_pred = k_pred,
      k_resp = k_resp,
      k_cca = k_cca
    ) %>%
      dplyr::filter(k_cca <= pmin(k_pred, k_resp))
  }

  cli::cli_inform("Evaluating {nrow(param_grid)} parameter combinations across {cv_folds$kfolds} folds")

  # Choose mapper function
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE)) {
      cli::cli_warn("furrr not available, falling back to sequential execution")
      map_fn <- purrr::map
    } else {
      map_fn <- furrr::future_map
    }
  } else {
    map_fn <- purrr::map
  }

  # Evaluate each parameter combination across all folds
  results <- map_fn(seq_len(nrow(param_grid)), function(i) {
    params <- param_grid[i, ]

    # Evaluate across folds
    fold_results <- purrr::map(cv_folds$folds, function(fold) {
      evaluate_fold(
        fold = fold,
        k_pred = params$k_pred,
        k_resp = params$k_resp,
        k_cca = params$k_cca,
        method = method,
        metrics = metrics
      )
    })

    # Combine fold results
    dplyr::bind_rows(fold_results) %>%
      dplyr::mutate(
        k_pred = params$k_pred,
        k_resp = params$k_resp,
        k_cca = params$k_cca,
        .before = 1
      )
  }, .progress = TRUE)

  dplyr::bind_rows(results)
}

#' Evaluate a single fold with given parameters
#'
#' @param fold A single fold list from cv_folds$folds
#' @param k_pred Number of predictor EOFs
#' @param k_resp Number of response EOFs
#' @param k_cca Number of CCA modes
#' @param metrics Metrics to compute
#'
#' @return Tibble with fold_id and metric values
#' @keywords internal
evaluate_fold <- function(fold, k_pred, k_resp, k_cca, method = "cca", metrics) {
  # Truncate patterns to requested k (cheap operation using [.patterns)
  pred_patterns <- fold$train_pred_patterns[1:k_pred]
  resp_patterns <- fold$train_resp_patterns[1:k_resp]

  # Couple patterns (k_cca is inert for method = "pcr")
  coupled <- couple(pred_patterns, resp_patterns, k = k_cca,
                    method = method, validate = FALSE)

  # Predict on test data
  predicted <- predict(coupled, fold$test_pred_data)

  # Compute metrics
  metric_values <- compute_spatial_metrics(predicted, fold$test_resp_data, metrics)

  # Return as single-row tibble
  result <- tibble::tibble(fold = fold$fold_id)
  for (m in names(metric_values)) {
    result[[m]] <- metric_values[[m]]
  }
  result
}

#' Summarize cross-validation results
#'
#' Aggregates results across folds and identifies best hyperparameters.
#'
#' @param cv_results A tibble from `tune_cca()`
#' @param metric Which metric to optimize (default "rmse")
#' @param minimize Logical, whether to minimize (TRUE for RMSE) or maximize
#'   (FALSE for correlations). Default TRUE.
#'
#' @return A tibble with mean and sd of each metric per parameter combination,
#'   sorted by the target metric. Best parameters are attached as attribute "best_params".
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- tune_cca(cv, k_pred = 1:5, k_resp = 1:5)
#' summary <- summarize_cv(results, metric = "rmse", minimize = TRUE)
#'
#' # Get best parameters (k_pred, k_resp, k_cca)
#' attr(summary, "best_params")
#' }
summarize_cv <- function(cv_results, metric = "rmse", minimize = TRUE) {
  if (!metric %in% names(cv_results)) {
    available <- setdiff(names(cv_results), c("k_pred", "k_resp", "k_cca", "fold"))
    cli::cli_abort("Metric '{metric}' not found. Available: {available}",
                   class = "tidyeof_invalid_metric")
  }

  # Find all metric columns
  metric_cols <- setdiff(names(cv_results), c("k_pred", "k_resp", "k_cca", "fold"))

  # Summarize across folds
  summary <- cv_results %>%
    dplyr::group_by(k_pred, k_resp, k_cca) %>%
    dplyr::summarize(
      dplyr::across(dplyr::all_of(metric_cols),
                    list(mean = ~mean(., na.rm = TRUE),
                         sd = ~sd(., na.rm = TRUE))),
      n_folds = dplyr::n(),
      .groups = "drop"
    )

  # Sort by target metric
  mean_col <- paste0(metric, "_mean")
  if (minimize) {
    summary <- dplyr::arrange(summary, .data[[mean_col]])
  } else {
    summary <- dplyr::arrange(summary, dplyr::desc(.data[[mean_col]]))
  }

  # Extract best params
  best <- summary[1, c("k_pred", "k_resp", "k_cca")]
  attr(summary, "best_params") <- as.list(best)

  summary
}

#' Print method for cv_folds objects
#' @param x A cv_folds object
#' @param ... Additional arguments (ignored)
#' @export
print.cv_folds <- function(x, ...) {
  cli::cli_h1("Cross-Validation Folds")
  cli::cli_text("Folds: {.field {x$kfolds}}")
  cli::cli_text("Common time steps: {.field {length(x$common_times)}}")
  cli::cli_text("Max predictor EOFs: {.field {x$max_k_pred}}")
  cli::cli_text("Max response EOFs: {.field {x$max_k_resp}}")

  if (!is.null(x$common_with_sources)) {
    cli::cli_text("Common EOF sources: {.val {x$common_with_sources}}")
  }

  cli::cli_h2("Pattern Options")
  cli::cli_text("Scale (predictor): {.field {x$pattern_opts$scale_pred}}")
  cli::cli_text("Scale (response): {.field {x$pattern_opts$scale_resp}}")
  cli::cli_text("Rotate: {.field {x$pattern_opts$rotate}}")
  cli::cli_text("Monthly: {.field {x$pattern_opts$monthly}}")
  cli::cli_text("Weight: {.field {x$pattern_opts$weight}}")

  invisible(x)
}

#' Cross-validate EOF truncation for a single field
#'
#' Evaluates reconstruction skill for different numbers of EOFs using k-fold
#' cross-validation with a speckled holdout. For each held-out fold a random
#' scatter of grid cells is hidden, mode amplitudes are estimated from the
#' visible cells, and the hidden cells are predicted. Because the hidden cells
#' are not used to estimate the amplitudes, prediction error is genuinely
#' out-of-sample and stops improving once `k` exceeds the field's effective
#' rank -- so the RMSE-minimising `k` is a meaningful estimate of how many
#' modes the data support (Bro et al. 2008).
#'
#' Naively projecting the full held-out field onto the EOFs and scoring the
#' reconstruction (the approach used before tidyeof 0.1.0) does NOT work: the
#' projection is least-squares optimal for the very data being scored, so error
#' decreases monotonically with `k` and the "best" `k` is always the largest.
#'
#' @param data A stars object with spatial-temporal data
#' @param k Vector of EOF counts to evaluate (default 1:10)
#' @param kfolds Number of cross-validation folds (default 5)
#' @param max_k Maximum EOFs to compute per fold (default max(k))
#' @param metrics Character vector of metrics to compute. Options:
#'   "rmse", "cor_spatial", "cor_temporal" (default: all three). Metrics are
#'   computed on the hidden cells only. For multivariate fields each metric
#'   also gets per-variable columns (e.g. `rmse_tmean`); the plain name is the
#'   pooled score (sd-normalized RMS for `rmse`, mean for correlations).
#'   Per-variable cross-validation scores for very small variable blocks are
#'   noisier and may be NA in some replicates, because the hidden cells are
#'   sampled across the concatenated space, proportional to block size.
#' @param scale Logical, whether to scale data before EOF extraction (default FALSE)
#' @param monthly Logical, whether to compute monthly climatology (default FALSE)
#' @param weight Logical, whether to apply area weighting (default TRUE)
#' @param hidden_fraction Fraction of grid cells to hide in each held-out fold
#'   (default 0.2). Hidden cells are predicted from the visible ones.
#' @param n_reps Number of random hidden-cell masks to average over per fold
#'   (default 5). More replicates give smoother, more stable estimates.
#' @param seed Base random seed for hidden-cell masks (default 1). Masks depend
#'   only on the fold and replicate, not on `k`, so all `k` are compared on the
#'   same hidden cells. The global RNG is left undisturbed.
#'
#' @return A tibble with columns: k, fold, and one column per metric.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Find optimal k for precipitation field
#' results <- tune_eof(precip_data, k = 1:15, kfolds = 5)
#' summary <- summarize_eof_cv(results, metric = "rmse")
#'
#' # Plot reconstruction skill vs k
#' library(ggplot2)
#' results %>%
#'   group_by(k) %>%
#'   summarize(rmse = mean(rmse)) %>%
#'   ggplot(aes(k, rmse)) + geom_line() + geom_point()
#' }
tune_eof <- function(data,
                     k = 1:10,
                     kfolds = 5,
                     max_k = max(k),
                     metrics = c("rmse", "cor_spatial", "cor_temporal"),
                     scale = FALSE,
                     monthly = FALSE,
                     weight = TRUE,
                     hidden_fraction = 0.2,
                     n_reps = 5,
                     seed = 1L) {

  times <- stars::st_get_dimension_values(data, "time")

  if (length(times) < kfolds * 2) {
    cli::cli_abort("Not enough time steps ({length(times)}) for {kfolds}-fold CV.",
                   class = "tidyeof_insufficient_data")
  }

  # Create fold assignments

  fold_times <- prep_folds(times, kfolds = kfolds)

  cli::cli_progress_step("Computing patterns for {kfolds} folds")

  # Build folds with pre-computed patterns at max_k
  folds <- purrr::imap(fold_times, function(test_times, fold_id) {
    train_data <- dplyr::filter(data, !(time %in% test_times))
    test_data <- dplyr::filter(data, time %in% test_times)

    train_patterns <- patterns(train_data, k = max_k,
                                   scale = scale, rotate = FALSE,
                                   monthly = monthly, weight = weight)

    list(
      train_patterns = train_patterns,
      test_data = test_data,
      fold_id = as.integer(fold_id)
    )
  })

  cli::cli_progress_done()

  cli::cli_inform("Evaluating {length(k)} k values across {kfolds} folds")

  # Evaluate each k value across all folds
  results <- purrr::map(k, function(k_val) {
    fold_results <- purrr::map(folds, function(fold) {
      evaluate_eof_fold(
        fold = fold,
        k = k_val,
        metrics = metrics,
        hidden_fraction = hidden_fraction,
        n_reps = n_reps,
        seed = seed
      )
    })

    dplyr::bind_rows(fold_results) %>%
      dplyr::mutate(k = k_val, .before = 1)
  }, .progress = TRUE)

  dplyr::bind_rows(results)
}

#' Evaluate one tuning argument with the global RNG temporarily seeded
#'
#' Runs `code` after `set.seed(seed)` and restores the previous RNG state on
#' exit, so cross-validation masks are reproducible without disturbing the
#' caller's random stream.
#' @keywords internal
with_seed <- function(seed, code) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }
  set.seed(seed)
  force(code)
}

#' Evaluate EOF reconstruction for a single fold via speckled holdout
#'
#' Hides a random scatter of grid cells in the held-out data, estimates mode
#' amplitudes from the visible cells by least squares, and scores the
#' prediction on the hidden cells. This makes reconstruction skill a genuine
#' out-of-sample quantity, so over-fitting with too many EOFs is penalised
#' (Bro et al. 2008). Masks depend only on the fold and replicate, not on `k`,
#' so every `k` is scored on the same hidden cells.
#'
#' @param fold A fold list containing train_patterns and test_data
#' @param k Number of EOFs to use
#' @param metrics Metrics to compute
#' @param hidden_fraction Fraction of valid grid cells to hide per replicate
#' @param n_reps Number of random hidden-cell masks to average over
#' @param seed Base seed; combined with the fold id so masks are reproducible
#'
#' @return Tibble with fold_id and metric values. For multivariate fields the
#'   metrics include pooled scores plus per-variable scores (e.g. `rmse_tmean`).
#' @keywords internal
evaluate_eof_fold <- function(fold, k, metrics, hidden_fraction = 0.2,
                              n_reps = 5, seed = 1L) {
  patterns_k <- fold$train_patterns[1:k]

  # Held-out anomalies in the space the patterns were fit in (own climatology)
  anomalies <- get_anomalies(fold$test_data, patterns_k$climatology,
                             scale = patterns_k$scaled, monthly = patterns_k$monthly)
  anom_mat <- flatten_time_space(units::drop_units(anomalies))$matrix

  valid <- patterns_k$valid_pixels
  n_valid <- length(valid)
  obs <- anom_mat[, valid, drop = FALSE]               # time x valid cells

  # Area weights applied during fitting must also weight the LS estimate;
  # one weight per grid cell, replicated across variable blocks
  n_vars <- length(patterns_k$names)
  w <- if (isTRUE(patterns_k$weight)) {
    rep(area_weights(fold$test_data), n_vars)[valid]
  } else {
    rep(1, n_valid)
  }

  block_map <- patterns_k$block_map
  if (is.null(block_map)) {
    block_map <- setNames(list(seq_len(ncol(anom_mat))), patterns_k$names[[1]])
  }

  weighted_loadings <- eof_loading_matrix(patterns_k)[valid, , drop = FALSE] * w
  weighted_obs <- sweep(obs, 2, w, `*`)                # time x valid cells

  n_hidden <- max(2L, round(hidden_fraction * n_valid))
  if (n_valid - n_hidden < k) {
    cli::cli_abort(
      "hidden_fraction = {hidden_fraction} leaves fewer than k = {k} visible cells.",
      class = "tidyeof_insufficient_cells"
    )
  }

  rep_metrics <- purrr::map(seq_len(n_reps), function(rep) {
    hidden <- with_seed(seed + fold$fold_id * 1000L + rep,
                        sample.int(n_valid, n_hidden))
    visible <- setdiff(seq_len(n_valid), hidden)

    # Estimate amplitudes from visible cells, predict the hidden ones
    amps <- t(qr.solve(weighted_loadings[visible, , drop = FALSE],
                       t(weighted_obs[, visible, drop = FALSE])))
    pred_hidden <- sweep(amps %*% t(weighted_loadings[hidden, , drop = FALSE]),
                         2, w[hidden], `/`)
    obs_hidden <- obs[, hidden, drop = FALSE]

    # Map hidden columns back to variable blocks for per-variable metrics
    hidden_cols <- valid[hidden]
    hidden_blocks <- purrr::map(block_map, ~which(hidden_cols %in% .x))
    hidden_blocks <- hidden_blocks[lengths(hidden_blocks) > 0]

    compute_block_metrics(pred_hidden, obs_hidden, hidden_blocks, metrics)
  })

  result <- tibble::tibble(fold = fold$fold_id)
  metric_names <- unique(unlist(purrr::map(rep_metrics, names)))
  for (m in metric_names) {
    vals <- vapply(rep_metrics,
                   function(v) if (is.null(v[[m]])) NA_real_ else v[[m]],
                   numeric(1))
    result[[m]] <- mean(vals, na.rm = TRUE)
  }
  result
}

#' Summarize EOF cross-validation results
#'
#' Aggregates results across folds and identifies optimal k.
#'
#' @param cv_results A tibble from `tune_eof()`
#' @param metric Which metric to optimize (default "rmse")
#' @param minimize Logical, whether to minimize (TRUE for RMSE) or maximize
#'   (FALSE for correlations). Default TRUE.
#'
#' @return A tibble with mean and sd of each metric per k value,
#'   sorted by the target metric. Best k is attached as attribute "best_k".
#'
#' @export
summarize_eof_cv <- function(cv_results, metric = "rmse", minimize = TRUE) {
  if (!metric %in% names(cv_results)) {
    available <- setdiff(names(cv_results), c("k", "fold"))
    cli::cli_abort("Metric '{metric}' not found. Available: {available}",
                   class = "tidyeof_invalid_metric")
  }

  # Find all metric columns
  metric_cols <- setdiff(names(cv_results), c("k", "fold"))

  # Summarize across folds
  summary <- cv_results %>%
    dplyr::group_by(k) %>%
    dplyr::summarize(
      dplyr::across(dplyr::all_of(metric_cols),
                    list(mean = ~mean(., na.rm = TRUE),
                         sd = ~sd(., na.rm = TRUE))),
      n_folds = dplyr::n(),
      .groups = "drop"
    )

  # Sort by target metric
  mean_col <- paste0(metric, "_mean")
  if (minimize) {
    summary <- dplyr::arrange(summary, .data[[mean_col]])
  } else {
    summary <- dplyr::arrange(summary, dplyr::desc(.data[[mean_col]]))
  }

  # Extract best k
  best_k <- summary$k[1]
  attr(summary, "best_k") <- best_k

  summary
}
