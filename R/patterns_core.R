# Input validation helpers ----

#' Check if input is a stars object
#' @param x Object to check
#' @param arg Argument name for error messages
#' @param call Calling environment for error messages
#' @keywords internal
check_stars_object <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (!inherits(x, "stars")) {
    cli::cli_abort(
      "Argument {.arg {arg}} must be a {.cls stars} object, not {.cls {class(x)}}.",
      class = "tidyeof_invalid_input",
      call = call
    )
  }
}

#' Check that a stars object has exactly one attribute
#'
#' EOF analysis operates on a single variable; silently using the first
#' attribute of a multi-attribute object would hide the others.
#' @param x A stars object
#' @param arg Argument name for error messages
#' @param call Calling environment for error messages
#' @keywords internal
check_single_attribute <- function(x, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
  if (length(x) > 1) {
    cli::cli_abort(
      c(
        "Argument {.arg {arg}} has {length(x)} attributes ({.field {names(x)}}), but only single-variable analysis is supported.",
        "i" = "Subset to one variable first, e.g. {.code {arg}[\"{names(x)[1]}\"]}."
      ),
      class = "tidyeof_multiple_attributes",
      call = call
    )
  }
}

#' Check that multivariate input is standardized
#'
#' PCA is variance-driven, so joint EOFs of variables with different units
#' require per-pixel standardization to contribute comparably.
#' @param dat A stars object
#' @param scale The scale argument passed to the caller
#' @param call Calling environment for error messages
#' @keywords internal
check_multivariate_scale <- function(dat, scale, call = rlang::caller_env()) {
  if (length(dat) > 1 && !isTRUE(scale)) {
    cli::cli_abort(
      c(
        "Multivariate input ({length(dat)} attributes: {.field {names(dat)}}) requires {.code scale = TRUE}.",
        "i" = "PCA is variance-driven; variables with different units must be standardized to contribute comparably."
      ),
      class = "tidyeof_multivariate_scale",
      call = call
    )
  }
}

#' Check if k is valid
#' @param k Number of components
#' @param max_k Maximum allowed components
#' @param arg Argument name for error messages
#' @param call Calling environment for error messages
#' @keywords internal
check_k_valid <- function(k, max_k, arg = rlang::caller_arg(k), call = rlang::caller_env()) {
  if (!is.numeric(k) || length(k) != 1 || k < 1 || k > max_k || k != floor(k)) {
    cli::cli_abort(
      "Argument {.arg {arg}} must be a single integer between 1 and {max_k}.",
      class = "tidyeof_invalid_k",
      call = call
    )
  }
}

#' Get EOFs and PCs from spatiotemporal data
#'
#' This function performs Empirical Orthogonal Function (EOF) analysis on
#' spatial-temporal data. For large datasets, it automatically uses IRLBA
#' (Implicitly Restarted Lanczos Bidiagonalization Algorithm) for efficient
#' computation when the irlba package is available.
#'
#' @param dat A `stars` object containing spatial and temporal dimensions.
#'   Multiple attributes (e.g. temperature and precipitation on the same grid)
#'   are analyzed jointly as combined EOFs: each variable contributes a block
#'   of the space dimension, modes share one amplitude time series, and
#'   `scale = TRUE` is required so variables with different units contribute
#'   comparably.
#' @param k The number of PC/EOF modes to retain
#' @param scale Logical, whether to scale before PCA
#' @param rotate Logical, whether to apply Varimax rotation. Rotation follows
#'   the standard REOF convention (Hannachi et al. 2007): varimax operates on
#'   sqrt(eigenvalue)-scaled EOFs without Kaiser row-normalization
#'   (`normalize = FALSE`), the stored patterns are the rotated loadings
#'   (unit norm, not mutually orthogonal), and amplitudes remain
#'   uncorrelated with sd = sqrt(rotated eigenvalue). Requires k > 1.
#' @param monthly Logical, whether to use monthly climatology
#' @param weight Logical, whether to apply area weighting
#' @param irlba_threshold Minimum number of data elements to trigger IRLBA usage
#'   (default: 500000). Set to Inf to always use base prcomp().
#'
#' @details
#' The eigenvalue table includes North et al. (1982) sampling-error bars
#' (`low`/`hi`), computed with a relative error of sqrt(2/n) where n is the
#' number of time steps. This assumes temporally independent samples: for
#' autocorrelated data (e.g., monthly anomalies) the effective sample size is
#' smaller and the bars are too narrow, so modes that appear well-separated
#' may not be.
#'
#' For multivariate input, per-pixel standardization is undefined where the
#' climatological standard deviation is ~0 (e.g. precipitation in arid cells);
#' such cells are dropped automatically like NA cells. Strongly skewed
#' variables such as precipitation often benefit from a sqrt or log transform
#' before analysis.
#'
#' @return A `patterns` object containing EOFs, amplitudes, and metadata
#' @export
patterns <- function(dat, k = 4, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE, irlba_threshold = 500000){

  # Input validation
  check_stars_object(dat)
  check_multivariate_scale(dat, scale)

  if (isTRUE(rotate) && k <= 1) {
    cli::cli_abort(
      "Rotation requires k > 1.",
      class = "tidyeof_invalid_option"
    )
  }

  # Capture units from original data before any modifications
  original_units <- setNames(purrr::map(names(dat), ~tryCatch(units(dat[[.x]]), error = function(e) NULL)), names(dat))

  # Get climatology and anomalies
  climatology <- get_climatology(dat, monthly = monthly)
  anomalies <- get_anomalies(dat, clim = climatology, scale = scale, monthly = monthly)

  # Compute spatial weights once so we can apply them just for the PCA step
  weights <- if (weight) area_weights(dat) else NULL

  # Get EOFs, which includes PCA/rotation handling internally
  eofs <- get_eofs(
    anomalies,
    k = k,
    rotate = rotate,
    irlba_threshold = irlba_threshold,
    weights = weights
  )

  # Create output pattern object using formal constructor
  patterns <- new_patterns(
    eofs = eofs$spatial_patterns,
    amplitudes = eofs$amplitudes,
    eigenvalues = eofs$eigenvalues,
    total_variance = eofs$total_variance,
    k = k,
    proj_matrix = eofs$proj_matrix,
    rotation = eofs$rotation_matrix,
    climatology = climatology,
    units = original_units,
    names = names(dat),
    scaled = scale,
    monthly = monthly,
    rotate = rotate,
    weight = weight,
    valid_pixels = eofs$valid_pixels,
    block_map = eofs$block_map
  )

  # Align patterns so EOFs have roughly similar dominant signs
  patterns <- flip_patterns(patterns)

  return(patterns)
}

#' Matrix-centric rotation helper that keeps all components synchronized
#'
#' Standard REOF convention (Hannachi et al. 2007): varimax operates on the
#' sqrt(eigenvalue)-scaled loadings, and the rotated scaled loadings ARE the
#' rotated patterns. They are stored unit-norm with the variance carried by
#' the scores, mirroring the unrotated convention (score sd = sqrt(eigenvalue));
#' rotated patterns are not mutually orthogonal, but scores stay uncorrelated.
#' @keywords internal
rotate_pca_components <- function(loadings_matrix, scores_matrix, sdev_vector) {
  # Scale eigenvectors by sqrt(eigenvalue) before rotation (Hannachi et al.
  # 2007). normalize = FALSE: the varimax criterion is applied to the scaled
  # loadings directly, without Kaiser row-normalization, per that convention.
  scaled_loadings <- loadings_matrix %*% diag(sdev_vector)
  rot <- varimax(scaled_loadings, normalize = FALSE)

  rotation_matrix <- rot$rotmat
  rotated_scaled_loadings <- unclass(rot$loadings)

  # Compute explained variance in rotated space for ordering
  rotated_eigenvals <- colSums(rotated_scaled_loadings^2)
  rotated_sdev <- sqrt(rotated_eigenvals)
  ev_order <- order(rotated_eigenvals, decreasing = TRUE)

  rotated_loadings <- sweep(rotated_scaled_loadings, 2, rotated_sdev, `/`)

  # Scores transform as Z R diag(rotated_sdev), with Z the standardized
  # scores: reconstruction scores %*% t(loadings) is preserved, and the same
  # matrix maps (weighted) anomalies to amplitudes — it is the least-squares
  # dual basis of the non-orthogonal rotated patterns.
  amplitude_transform <- sweep(sweep(rotation_matrix, 1, sdev_vector, `/`),
                               2, rotated_sdev, `*`)
  rotated_scores <- scores_matrix %*% amplitude_transform

  list(
    loadings = rotated_loadings[, ev_order, drop = FALSE],
    scores = rotated_scores[, ev_order, drop = FALSE],
    sdev = rotated_sdev[ev_order],
    eigenvalues = rotated_eigenvals[ev_order],
    rotation_matrix = rotation_matrix[, ev_order, drop = FALSE],
    amplitude_transform = amplitude_transform[, ev_order, drop = FALSE]
  )
}

#' Internal function to calculate EOFs and related components
#' @param weights Optional numeric vector of spatial weights (one per spatial
#'   location) to be applied to the anomalies before PCA
#' @keywords internal
get_eofs <- function(dat, k, rotate = FALSE, irlba_threshold, weights = NULL) {
  times <- stars::st_get_dimension_values(dat, "time")

  pc_names <- names0(k, 'PC')

  var_names <- names(dat)
  n_vars <- length(var_names)

  flattened <- flatten_time_space(units::drop_units(dat))
  anomaly_matrix_full <- flattened$matrix
  n_pixels <- ncol(anomaly_matrix_full)   # n_vars * n_space
  n_space <- flattened$n_space
  spatial_shape <- flattened$spatial_shape
  block_map <- flattened$block_map

  # Valid pixels must be finite everywhere: this drops NA-masked cells and
  # also Inf/NaN cells produced by sd ~ 0 standardization (e.g. arid cells
  # for precipitation)
  valid_pixels <- which(apply(anomaly_matrix_full, 2,
                              function(col) all(is.finite(col))))

  # Validate k value
  max_k <- min(length(times) - 1, length(valid_pixels))
  check_k_valid(k, max_k)

  # Extract matrix for valid pixels (time x V*space)
  anomaly_matrix <- anomaly_matrix_full[, valid_pixels, drop = FALSE]

  # Apply spatial weights column-wise if provided; one weight per grid cell,
  # replicated across variable blocks
  if (!is.null(weights)) {
    if (length(weights) != n_space) {
      cli::cli_abort(
        "Length of {.arg weights} ({length(weights)}) must match number of spatial points per variable ({n_space}).",
        class = "tidyeof_weight_mismatch"
      )
    }
    weights_valid <- rep(weights, n_vars)[valid_pixels]
    anomaly_matrix <- sweep(anomaly_matrix, 2, weights_valid, `*`)
  } else {
    weights_valid <- rep(1, length(valid_pixels))
  }

  # Perform PCA on (optionally) weighted anomalies without centering (already anomalies)
  pca_result <- perform_pca_smart(
    anomaly_matrix,
    k = k,
    center = FALSE,
    size_threshold = irlba_threshold
  )

  # Keep everything as matrices initially for synchronized operations
  loadings_matrix <- pca_result$rotation[, 1:k, drop = FALSE]
  scores_matrix <- pca_result$x[, 1:k, drop = FALSE]
  sdev_vector <- pca_result$sdev[1:k]

  rotation_matrix <- NULL

  # Handle rotation with synchronized matrix operations
  if (rotate && k > 1) {
    rotation_result <- rotate_pca_components(loadings_matrix, scores_matrix, sdev_vector)
    loadings_weighted <- rotation_result$loadings
    amplitudes <- rotation_result$scores
    rotation_matrix <- rotation_result$rotation_matrix
    component_sdev <- rotation_result$sdev
    component_variance <- rotation_result$eigenvalues
    # Rotated patterns are not orthogonal, so projection onto them uses the
    # least-squares dual basis rather than the patterns themselves
    proj_weighted <- loadings_matrix %*% rotation_result$amplitude_transform
  } else {
    loadings_weighted <- loadings_matrix
    amplitudes <- scores_matrix
    component_sdev <- sdev_vector
    component_variance <- sdev_vector^2
    proj_weighted <- loadings_matrix
  }

  # Convert loadings back to physical space so EOFs carry interpretable units
  loadings <- sweep(loadings_weighted, 1, weights_valid, `/`)

  # Create EOF spatial patterns
  full_patterns <- array(NA, dim = c(k, n_pixels))
  full_patterns[, valid_pixels] <- t(loadings)

  # Build a multi-attribute template (one attribute per variable) with the
  # time dimension relabeled as PC, then fill each attribute with its block.
  # Attribute names come from the data, so univariate EOFs are named after
  # their variable (previously "weight").
  if (has_geometry_dimension(dat)) {
    template <- dat[, , 1:k, drop = FALSE] %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC')
    for (i in seq_along(var_names)) {
      template[[i]] <- t(full_patterns[, block_map[[i]], drop = FALSE])  # geometry x PC
    }
  } else {
    template <- dat[, , , 1:k, drop = FALSE] %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC')
    for (i in seq_along(var_names)) {
      pattern_array <- array(full_patterns[, block_map[[i]], drop = FALSE],
                             dim = c(k, spatial_shape))
      template[[i]] <- aperm(pattern_array, c(2, 3, 1))  # x, y, PC
    }
  }
  spatial_patterns <- template

  # Format amplitudes (already rotated and reordered if needed)
  colnames(amplitudes) <- pc_names
  amplitudes <- amplitudes %>%
    as_tibble() %>%
    mutate(time = times, .before = 1)

  # Calculate eigenvalues - always use original unrotated values for scree plot
  # prcomp_irlba returns only the leading k singular values, so percent/low/hi
  # must use the true total variance ($totalvar), and the North et al. (1982)
  # sampling error needs the number of temporal samples, not length(sdev)
  total_var <- pca_result$totalvar
  if (is.null(total_var)) total_var <- sum(pca_result$sdev^2)
  n_times <- length(times)
  eigenvalues <- tidy_pca_sdev(pca_result, total_var) |>
    mutate(eigenvalues = std.dev ^ 2,
           percent = percent * 100,
           cumulative = cumulative * 100,
           error = sqrt(2 / n_times),
           low =  eigenvalues * (1 - error) * 100 / total_var,
           hi = eigenvalues * (1 + error) * 100 / total_var)

  if (rotate && k > 1) {
    # Replace variance stats for the retained modes with the rotated values so
    # downstream scaling/plots stay in sync with the reordered amplitudes.
    # (See Hannachi et al. 2007 re. Kaiser-normalised rotation preserving total variance.)
    rotated_percent <- component_variance / total_var * 100
    rotated_cumulative <- cumsum(rotated_percent)

    eigenvalues <- eigenvalues %>%
      mutate(
        std.dev = if_else(PC <= k, component_sdev[PC], std.dev),
        eigenvalues = if_else(PC <= k, component_variance[PC], eigenvalues),
        percent = if_else(PC <= k, rotated_percent[PC], percent),
        cumulative = if_else(PC <= k, rotated_cumulative[PC], cumulative)
      )
  }

  list(
    spatial_patterns = spatial_patterns,
    amplitudes = amplitudes,
    eigenvalues = eigenvalues,
    total_variance = total_var,
    rotation_matrix = rotation_matrix,
    valid_pixels = valid_pixels,
    block_map = block_map,
    spatial_dims = flattened$spatial_dims,
    spatial_shape = flattened$spatial_shape,
    # proj_matrix maps weighted anomalies to amplitudes (the dual basis of the
    # patterns; equal to the patterns themselves only when unrotated)
    proj_matrix = proj_weighted
  )
}

#' Check if stars object has sf geometry
#' @param dat A stars object
#' @return Logical indicating if object has geometry dimension
#' @keywords internal
has_geometry_dimension <- function(dat) {
  "geometry" %in% names(stars::st_dimensions(dat))
}

#' Get spatial dimension names for stars object
#' @param dat A stars object
#' @return Character vector of spatial dimension names
#' @keywords internal
get_spatial_dimensions <- function(dat) {
  if (has_geometry_dimension(dat)) {
    return("geometry")
  }

  dim_names <- names(stars::st_dimensions(dat))
  x_patterns <- c("x", "lon", "longitude", "easting")
  y_patterns <- c("y", "lat", "latitude", "northing")

  x_dim <- dim_names[tolower(dim_names) %in% x_patterns][1]
  y_dim <- dim_names[tolower(dim_names) %in% y_patterns][1]

  if (is.na(x_dim) || is.na(y_dim)) {
    cli::cli_abort(
      "Could not find spatial dimensions. Expected 'x'/'lon'/'longitude'/'easting' and 'y'/'lat'/'latitude'/'northing', or 'geometry'.",
      class = "tidyeof_invalid_dimensions"
    )
  }

  c(x_dim, y_dim)
}

#' Compute area-based weights for spatial data
#'
#' Calculates area weights for spatial data using st_area(). Works uniformly
#' for raster grids, irregular geometries, and different coordinate systems.
#'
#' @param dat A stars object with spatial dimensions
#' @return Numeric weights (sqrt of normalized areas)
#' @export
area_weights <- function(dat) {
  if (has_geometry_dimension(dat)) {
    # For sf geometry, extract the geometry and compute areas directly
    geom <- sf::st_geometry(dat)
    areas <- sf::st_area(geom)
    area_values <- as.numeric(areas)
  } else {
    # For raster grids, use st_area on the stars object
    areas <- sf::st_area(dat)
    area_values <- as.numeric(areas[[1]])
  }

  # Normalize by mean area and take sqrt
  # sqrt so the covariance matrix is weighted by area
  mean_area <- mean(area_values, na.rm = TRUE)
  if (!is.finite(mean_area) || mean_area <= 0) {
    cli::cli_abort(
      c(
        "Cannot compute area weights: geometries have zero or undefined area.",
        "i" = "This is expected for POINT or LINESTRING geometries (e.g., station data).",
        "i" = "Use {.code weight = FALSE} for data without areal grid cells."
      ),
      class = "tidyeof_zero_area"
    )
  }
  sqrt(area_values / mean_area)
}

#' Tidy PCA standard deviations into a tibble
#'
#' Replaces `broom::tidy(pca, matrix = "pcs")` with a dependency-free version.
#' Returns a tibble with columns PC, std.dev, percent, cumulative.
#'
#' @param pca_result A prcomp (or prcomp_irlba) result
#' @param total_var Total variance of the data (sum of all eigenvalues).
#'   Defaults to `sum(pca_result$sdev^2)`, which is only correct when the
#'   full spectrum is present; pass `pca_result$totalvar` for truncated
#'   (IRLBA) results.
#' @return Tibble with PC, std.dev, percent, cumulative
#' @keywords internal
tidy_pca_sdev <- function(pca_result, total_var = NULL) {
  sdev <- pca_result$sdev
  variance <- sdev^2
  if (is.null(total_var)) total_var <- sum(variance)
  pct <- variance / total_var
  tibble::tibble(
    PC = seq_along(sdev),
    std.dev = sdev,
    percent = pct,
    cumulative = cumsum(pct)
  )
}

# from tidymodels/recipes
names0 <- function(num, prefix = "PC") {
  if (num < 1) {
    rlang::abort("`k` should be > 0.")
  }
  ind <- format(seq_len(num))
  ind <- gsub(" ", "0", ind)
  paste0(prefix, ind)
}

#' Smart PCA Selection with IRLBA Support
#'
#' Automatically selects between base `prcomp()` and `prcomp_irlba()` based on
#' data size. For large datasets, IRLBA provides significant computational
#' savings when only the first few components are needed.
#'
#' @param x A numeric matrix for PCA computation
#' @param k Number of components to compute
#' @param center Logical, whether to center the data
#' @param scale. Logical, whether to scale the data
#' @param size_threshold Minimum number of elements to trigger IRLBA
#' @param ... Additional arguments passed to the PCA function
#'
#' @return A PCA result object compatible with `prcomp()` output
#' @keywords internal
perform_pca_smart <- function(x, k = NULL, center = TRUE, scale. = FALSE,
                              size_threshold, ...) {
  data_size <- nrow(x) * ncol(x)

  if (data_size >= size_threshold && !is.null(k)) {
    cli::cli_inform(
      "Using IRLBA for efficient PCA computation on large dataset ({format(data_size, big.mark = ',')} elements)."
    )
    irlba::prcomp_irlba(x, n = k, center = center, scale. = scale., ...)
  } else {
    prcomp(x, center = center, scale. = scale., ...)
  }
}
