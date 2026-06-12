# Coupling functions for EOF patterns via CCA
# Consolidated from couple_patterns.R and predict_coupled_patterns.R

#' Fit the OLS coefficient matrix for PCR coupling
#'
#' Ordinary least squares of the predictand PC amplitudes on the predictor PC
#' amplitudes. The EOF truncation already performed the dimension reduction, so
#' this is plain (optionally centered) multivariate regression. The fit is
#' rank-safe: a rank-deficient predictor matrix (e.g. `k_pred >= n_times` or
#' collinear amplitudes) aborts rather than returning `NA` coefficients.
#'
#' @param pred_amps Predictor amplitude matrix (time x k_pred)
#' @param resp_amps Predictand amplitude matrix (time x k_resp)
#' @param center Logical; center both sides before fitting (default TRUE)
#' @return List with `coefficients` (k_pred x k_resp), `xcenter`, `ycenter`
#'   (each a named numeric vector when centered, or `FALSE`)
#' @keywords internal
fit_pcr <- function(pred_amps, resp_amps, center = TRUE) {
  if (isTRUE(center)) {
    xcenter <- colMeans(pred_amps)
    ycenter <- colMeans(resp_amps)
    Xc <- sweep(pred_amps, 2, xcenter, "-")
    Yc <- sweep(resp_amps, 2, ycenter, "-")
  } else {
    xcenter <- FALSE
    ycenter <- FALSE
    Xc <- pred_amps
    Yc <- resp_amps
  }

  qrX <- qr(Xc)
  if (qrX$rank < ncol(Xc)) {
    cli::cli_abort(
      c(
        "Predictor amplitudes are rank-deficient ({qrX$rank} < {ncol(Xc)} columns); OLS is not identifiable.",
        "i" = "Reduce the number of predictor EOFs so it is below the number of time steps and free of collinearity."
      ),
      class = "tidyeof_rank_deficient"
    )
  }

  coefficients <- qr.coef(qrX, Yc)
  list(coefficients = coefficients, xcenter = xcenter, ycenter = ycenter)
}

#' Couple Pattern Relationships Using CCA
#'
#' This function couples predictor and response patterns using Canonical Correlation Analysis (CCA)
#' as the primary method. CCA finds linear combinations of predictor and response patterns that
#' maximize correlation between them.
#'
#' @param predictor_patterns A patterns object containing predictor patterns (e.g., from patterns())
#' @param response_patterns A patterns object containing response patterns (e.g., from patterns())
#' @param k Number of CCA modes to retain. If NULL, uses min(ncol(predictor), ncol(response))
#' @param method Coupling method. Currently only "cca" is supported
#' @param center Logical, whether to center the amplitudes before CCA
#'   (default: TRUE). Centering is the statistically standard choice and makes
#'   retaining all modes equivalent to multivariate regression with an
#'   intercept. It is a no-op when the amplitudes are already zero-mean over
#'   the coupled period (the usual case), but is essential when the predictor
#'   and response patterns were fit on different periods and then filtered to a
#'   common one, which leaves the common-period amplitudes with a nonzero mean.
#' @param validate Logical, whether to validate input patterns compatibility
#'
#' @return A coupled_patterns object containing:
#'   \item{cca}{The CCA results from cancor()}
#'   \item{predictor_patterns}{The original predictor patterns}
#'   \item{response_patterns}{The original response patterns}
#'   \item{k}{Number of CCA modes retained}
#'   \item{method}{Coupling method used}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Get patterns from your data
#' pred_patterns <- patterns(predictor_data, k = 5)
#' resp_patterns <- patterns(response_data, k = 5)
#'
#' # Couple the patterns
#' coupled <- couple(pred_patterns, resp_patterns, k = 3)
#'
#' # Make predictions
#' predictions <- predict(coupled, new_predictor_data)
#' }
couple <- function(predictor_patterns, response_patterns, k = NULL,
                  method = "cca", center = TRUE, validate = TRUE) {

  if (!method %in% c("cca", "pcr")) {
    cli::cli_abort(
      "Unsupported coupling method {.val {method}}. Use {.val cca} or {.val pcr}.",
      class = "tidyeof_unsupported_method"
    )
  }

  # Validate inputs and get common times
  common_times <- if (validate) {
    validate_patterns_compatibility(predictor_patterns, response_patterns)
  } else {
    pred_times <- get_times(predictor_patterns)
    resp_times <- get_times(response_patterns)
    pred_times[pred_times %in% resp_times]
  }

  # Extract amplitude matrices filtered to common times
  pred_amps <- extract_amplitudes_matrix(predictor_patterns, common_times)
  resp_amps <- extract_amplitudes_matrix(response_patterns, common_times)

  # Notify user if filtering occurred
  pred_n <- length(get_times(predictor_patterns))
  resp_n <- length(get_times(response_patterns))
  if (length(common_times) < pred_n || length(common_times) < resp_n) {
    cli::cli_inform("Filtered to {length(common_times)} common time steps (predictor had {pred_n}, response had {resp_n}).")
  }

  coupled <- list(
    predictor_patterns = predictor_patterns,
    response_patterns = response_patterns,
    method = method,
    center = center
  )

  if (method == "cca") {
    # Determine and validate k (number of canonical modes)
    if (is.null(k)) {
      k <- min(ncol(pred_amps), ncol(resp_amps))
    }
    max_k <- min(ncol(pred_amps), ncol(resp_amps))
    if (k > max_k) {
      warning("k = ", k, " exceeds maximum possible (", max_k, "). Setting k = ", max_k)
      k <- max_k
    }
    coupled$cca <- cancor(pred_amps, resp_amps, xcenter = center, ycenter = center)
    coupled$k <- k
  } else {  # method == "pcr"
    # The coupling-k is inert for PCR; regularization is the predictor
    # truncation (k_pred). k is stored only for reporting.
    coupled$pcr <- fit_pcr(pred_amps, resp_amps, center = center)
    coupled$k <- ncol(pred_amps)
  }

  class(coupled) <- "coupled_patterns"
  return(coupled)
}

#' Print method for coupled_patterns
#' @param x A coupled_patterns object
#' @param ... Additional arguments (ignored)
#' @export
print.coupled_patterns <- function(x, ...) {
  cli::cli_h1("Coupled Patterns Object")
  cli::cli_text("Method: {.field {x$method}}")
  if (x$method == "cca") {
    cli::cli_text("CCA modes retained: {.field {x$k}}")
    cli::cli_text("Canonical correlations: {.val {round(x$cca$cor[1:x$k], 3)}}")
  } else {
    cli::cli_text("Predictor PCs used: {.field {x$k}}")
  }
  cli::cli_text("Predictor patterns: {.field {ncol(extract_amplitudes_matrix(x$predictor_patterns))}} PCs")
  cli::cli_text("Response patterns: {.field {ncol(extract_amplitudes_matrix(x$response_patterns))}} PCs")
  invisible(x)
}

#' Summary method for coupled_patterns
#' @param object A coupled_patterns object
#' @param ... Additional arguments (ignored)
#' @export
summary.coupled_patterns <- function(object, ...) {
  cli::cli_h1("Coupled Patterns Summary")
  cli::cli_text("Method: {.field {object$method}}")
  cli::cli_text("Centered: {.field {object$center}}")
  if (object$method == "cca") {
    cli::cli_text("CCA modes retained: {.field {object$k}}")
    cli::cli_h2("Canonical Correlations")
    print(get_canonical_correlations(object))
  } else {
    cli::cli_text("Predictor PCs used: {.field {object$k}}")
    cli::cli_text("Response PCs predicted: {.field {ncol(object$pcr$coefficients)}}")
  }
  invisible(object)
}

# Helper function to validate pattern compatibility and return common times
validate_patterns_compatibility <- function(predictor_patterns, response_patterns) {
  if (!inherits(predictor_patterns, "patterns") && !is.data.frame(predictor_patterns)) {
    cli::cli_abort("predictor_patterns must be a patterns object or data frame with amplitudes",
                   class = "tidyeof_invalid_input")
  }

  if (!inherits(response_patterns, "patterns") && !is.data.frame(response_patterns)) {
    cli::cli_abort("response_patterns must be a patterns object or data frame with amplitudes",
                   class = "tidyeof_invalid_input")
  }

  pred_times <- get_times(predictor_patterns)
  resp_times <- get_times(response_patterns)

  # Use match-based intersection to preserve Date/POSIXct class
  common_times <- pred_times[pred_times %in% resp_times]

  if (length(common_times) == 0) {
    cli::cli_abort("No common time steps found between predictor and response patterns",
                   class = "tidyeof_no_common_times")
  }

  if (length(common_times) < 10) {
    cli::cli_warn("Only {length(common_times)} common time steps found. Consider using more data.")
  }

  common_times
}

# Prediction methods ----

#' Predict Method for Coupled Patterns
#'
#' This function makes predictions using a coupled_patterns object created by couple().
#' It applies the learned CCA relationship to new predictor data to predict response patterns.
#'
#' @param object A coupled_patterns object from couple()
#' @param newdata New predictor data (stars object) for making predictions
#' @param k Number of CCA modes to use for prediction. If NULL, uses all available modes
#' @param reconstruct Logical, whether to reconstruct the full spatial field (default: TRUE)
#' @param predictor_patterns Optional patterns object to use instead of the one stored
#'   in the coupled object. Useful for cross-source prediction with common EOFs: the
#'   override patterns share the same EOF space but carry a different climatology.
#' @param ... Additional arguments (currently unused)
#'
#' @return If reconstruct=TRUE, returns a stars object with reconstructed spatial fields.
#'         If reconstruct=FALSE, returns a tibble with predicted amplitudes.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create coupled patterns
#' coupled <- couple(pred_patterns, resp_patterns, k = 3)
#'
#' # Make predictions on new data
#' predictions <- predict(coupled, new_predictor_data)
#'
#' # Just get predicted amplitudes without spatial reconstruction
#' amplitudes <- predict(coupled, new_predictor_data, reconstruct = FALSE)
#'
#' # Cross-source prediction with common EOFs
#' cpat <- common_patterns(list(era = era, phyda = phyda), k = 5)
#' coupled <- couple(cpat$era, fine_patterns, k = 3)
#' predict(coupled, phyda_new, predictor_patterns = cpat$phyda)
#' }
predict.coupled_patterns <- function(object, newdata, k = NULL, reconstruct = TRUE,
                                   predictor_patterns = NULL, ...) {

  if (!inherits(object, "coupled_patterns")) {
    cli::cli_abort("object must be a coupled_patterns object from couple()",
                   class = "tidyeof_invalid_input")
  }

  if (!object$method %in% c("cca", "pcr")) {
    cli::cli_abort("Unsupported coupling method {.val {object$method}}.",
                   class = "tidyeof_unsupported_method")
  }

  # Use override patterns if provided, otherwise use stored patterns
  proj_patterns <- predictor_patterns %||% object$predictor_patterns

  # Validate compatibility if overriding
  if (!is.null(predictor_patterns)) {
    if (!identical(dim(predictor_patterns$proj_matrix),
                   dim(object$predictor_patterns$proj_matrix))) {
      cli::cli_abort(
        "{.arg predictor_patterns} must share the same EOF space as the training predictor patterns.",
        class = "tidyeof_incompatible_patterns"
      )
    }
  }

  # Project new data onto predictor patterns to get PC amplitudes
  new_amplitudes <- project_patterns(proj_patterns, newdata)

  predicted_amplitudes <- if (object$method == "cca") {
    if (is.null(k)) {
      k <- object$k
    }
    if (k > object$k) {
      warning("Requested k (", k, ") exceeds available modes (", object$k, "). Using k = ", object$k)
      k <- object$k
    }
    apply_cca_prediction(new_amplitudes = new_amplitudes, cca_result = object$cca, k = k)
  } else {
    # k is inert for PCR (regularization is the predictor truncation at couple() time)
    if (!is.null(k)) {
      cli::cli_warn(
        "{.arg k} is ignored for a PCR coupling; predictor truncation is fixed at {.fn couple} time.",
        class = "tidyeof_k_ignored"
      )
    }
    apply_pcr_prediction(new_amplitudes = new_amplitudes, pcr = object$pcr)
  }

  if (!reconstruct) {
    return(predicted_amplitudes)
  }

  reconstruct(target_patterns = object$response_patterns,
              amplitudes = predicted_amplitudes)
}

#' Apply CCA Prediction Transform
#'
#' Internal function that applies the CCA transformation to predict response amplitudes
#' from predictor amplitudes.
#'
#' @param new_amplitudes Tibble with time and predictor PC amplitudes
#' @param cca_result CCA result object from cancor()
#' @param k Number of CCA modes to use
#'
#' @return Tibble with time and predicted response PC amplitudes
#'
#' @keywords internal
apply_cca_prediction <- function(new_amplitudes, cca_result, k) {

  # Extract times
  new_times <- new_amplitudes$time

  # Convert to matrix for CCA transformation
  pred_matrix <- new_amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  # Apply training centering if CCA was fit with centering
  if (!identical(cca_result$xcenter, FALSE)) {
    xcenter <- cca_result$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(pred_matrix))) {
      xcenter <- xcenter[colnames(pred_matrix)]
    }
    pred_matrix <- sweep(pred_matrix, 2, xcenter, "-")
  }

  # Apply CCA transformation:
  # 1. Transform predictors to canonical variables
  # 2. Apply canonical correlations
  # 3. Transform back to response space
  canonical_predictors <- pred_matrix %*% cca_result$xcoef[, 1:k, drop = FALSE]
  canonical_responses <- canonical_predictors %*% diag(cca_result$cor[1:k], nrow = k)

  # Transform canonical responses back to PC space: regression of response
  # amplitudes on the retained canonical variates (Glahn 1968). The variates
  # are orthonormal, so the regression coefficients are the leading k rows of
  # the pseudo-inverse of the FULL ycoef. (A pseudo-inverse of the truncated
  # ycoef would give a minimum-norm preimage instead of the regression,
  # degrading predictions whenever k < ncol(ycoef).)
  response_amplitudes <- canonical_responses %*%
    MASS::ginv(cca_result$ycoef)[seq_len(k), , drop = FALSE]

  # Add back response centering if used during training
  if (!identical(cca_result$ycenter, FALSE)) {
    ycenter <- cca_result$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(response_amplitudes))) {
      ycenter <- ycenter[colnames(response_amplitudes)]
    }
    response_amplitudes <- sweep(response_amplitudes, 2, ycenter, "+")
  }

  # Convert back to tibble with proper column names
  n_response_pcs <- ncol(response_amplitudes)
  pc_names <- paste0("PC", 1:n_response_pcs)

  result <- response_amplitudes %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(pc_names) %>%
    mutate(time = new_times, .before = 1)

  return(result)
}

#' Apply PCR Prediction Transform
#'
#' Internal function that maps predictor PC amplitudes to predicted predictand
#' PC amplitudes via the fitted OLS coefficient matrix.
#'
#' @param new_amplitudes Tibble with `time` and predictor PC amplitudes
#' @param pcr The `pcr` slot of a coupled object: `coefficients`, `xcenter`, `ycenter`
#' @return Tibble with `time` and predicted response PC amplitudes
#' @keywords internal
apply_pcr_prediction <- function(new_amplitudes, pcr) {
  new_times <- new_amplitudes$time

  pred_matrix <- new_amplitudes %>%
    dplyr::select(-time) %>%
    as.matrix()

  # Apply training centering if the fit was centered (mirrors apply_cca_prediction)
  if (!identical(pcr$xcenter, FALSE)) {
    xcenter <- pcr$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(pred_matrix))) {
      xcenter <- xcenter[colnames(pred_matrix)]
    }
    pred_matrix <- sweep(pred_matrix, 2, xcenter, "-")
  }

  response_amplitudes <- pred_matrix %*% pcr$coefficients

  if (!identical(pcr$ycenter, FALSE)) {
    ycenter <- pcr$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(response_amplitudes))) {
      ycenter <- ycenter[colnames(response_amplitudes)]
    }
    response_amplitudes <- sweep(response_amplitudes, 2, ycenter, "+")
  }

  n_response_pcs <- ncol(response_amplitudes)
  pc_names <- paste0("PC", 1:n_response_pcs)

  response_amplitudes %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(pc_names) %>%
    mutate(time = new_times, .before = 1)
}

# CCA accessors ----

#' Get Canonical Variables from Coupled Patterns
#'
#' Extract canonical variables from either predictor or response patterns
#'
#' @param object A coupled_patterns object
#' @param data Original data (patterns object or amplitudes tibble)
#' @param type Either "predictor" or "response"
#' @param k Number of canonical modes to extract
#'
#' @return Tibble with canonical variables
#'
#' @export
get_canonical_variables <- function(object, data, type = c("predictor", "response"), k = NULL) {

  type <- match.arg(type)

  if (is.null(k)) {
    k <- object$k
  }

  # Get transformation coefficients
  if (type == "predictor") {
    coef_matrix <- object$cca$xcoef[, 1:k, drop = FALSE]
  } else {
    coef_matrix <- object$cca$ycoef[, 1:k, drop = FALSE]
  }

  # Apply transformation
  times <- get_times(data)
  amp_matrix <- extract_amplitudes_matrix(data)

  # Apply training centering so canonical variables match cancor() inputs
  if (type == "predictor" && !identical(object$cca$xcenter, FALSE)) {
    xcenter <- object$cca$xcenter
    if (!is.null(names(xcenter)) && !is.null(colnames(amp_matrix))) {
      xcenter <- xcenter[colnames(amp_matrix)]
    }
    amp_matrix <- sweep(amp_matrix, 2, xcenter, "-")
  }
  if (type == "response" && !identical(object$cca$ycenter, FALSE)) {
    ycenter <- object$cca$ycenter
    if (!is.null(names(ycenter)) && !is.null(colnames(amp_matrix))) {
      ycenter <- ycenter[colnames(amp_matrix)]
    }
    amp_matrix <- sweep(amp_matrix, 2, ycenter, "-")
  }

  canonical_vars <- amp_matrix %*% coef_matrix

  # Return as tibble
  canonical_names <- paste0("CV", 1:k)
  result <- canonical_vars %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(canonical_names) %>%
    mutate(time = times, .before = 1)

  return(result)
}

#' Get Canonical Spatial Patterns from Coupled Patterns
#'
#' Computes the spatial patterns corresponding to each canonical mode by
#' taking linear combinations of the original EOFs weighted by CCA coefficients.
#' These are the spatial patterns that, when projected onto the data, yield
#' the canonical variates.
#'
#' @param object A coupled_patterns object
#' @param type Either "predictor" or "response"
#' @param k Number of canonical modes to extract (default: all available)
#'
#' @return A stars object with canonical spatial patterns (dimension "CV" instead of "PC")
#'
#' @export
#'
#' @examples
#' \dontrun{
#' coupled <- couple(pred_patterns, resp_patterns, k = 3)
#'
#' # Get canonical patterns for response side
#' resp_canonical <- get_canonical_patterns(coupled, type = "response")
#' plot(resp_canonical)
#'
#' # Compare to original EOFs
#' plot(coupled$response_patterns$eofs)
#' }
get_canonical_patterns <- function(object, type = c("predictor", "response"), k = NULL) {

  type <- match.arg(type)

  if (is.null(k)) {
    k <- object$k
  }

  # Get the patterns object and CCA coefficients
  if (type == "predictor") {
    patterns <- object$predictor_patterns
    coef_matrix <- object$cca$xcoef[, 1:k, drop = FALSE]
  } else {
    patterns <- object$response_patterns
    coef_matrix <- object$cca$ycoef[, 1:k, drop = FALSE]
  }

  # EOF attribute arrays are constructed spatial-first with PC last
  # (see get_eofs), one attribute per variable
  eof_stars <- patterns$eofs
  eof_dims <- stars::st_dimensions(eof_stars)
  spatial_dim_names <- setdiff(names(eof_dims), "PC")
  n_pcs <- length(stars::st_get_dimension_values(eof_stars, "PC"))

  new_dims <- eof_dims[spatial_dim_names]
  cv_dim <- list(
    from = 1L,
    to = k,
    offset = NA_real_,
    delta = NA_real_,
    refsys = NA_character_,
    point = FALSE,
    values = paste0("CV", 1:k)
  )
  class(cv_dim) <- "dimension"
  new_dims$CV <- cv_dim
  class(new_dims) <- "dimensions"

  result_list <- purrr::map(names(eof_stars), function(v) {
    eof_array <- eof_stars[[v]]
    spatial_shape <- dim(eof_array)[-length(dim(eof_array))]
    eof_matrix <- matrix(eof_array, nrow = prod(spatial_shape), ncol = n_pcs)

    # canonical_pattern[i] = sum_j EOF[j] * coef[j, i]
    canonical_array <- array(eof_matrix %*% coef_matrix,
                             dim = c(spatial_shape, k))
    setNames(stars::st_as_stars(canonical_array, dimensions = new_dims), v)
  })

  do.call(c, result_list)
}

#' Get Canonical Correlations from Coupled Patterns
#'
#' Extract the canonical correlations and related statistics
#'
#' @param object A coupled_patterns object
#' @param k Number of modes to return (default: all available)
#'
#' @return Data frame with canonical correlation statistics
#'
#' @export
get_canonical_correlations <- function(object, k = NULL) {

  if (is.null(k)) {
    k <- object$k
  }

  correlations <- object$cca$cor[1:k]

  result <- data.frame(
    mode = 1:k,
    correlation = correlations,
    correlation_squared = correlations^2
  )

  return(result)
}
