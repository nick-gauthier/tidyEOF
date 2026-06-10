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
#' @param dat A `stars` object containing spatial and temporal dimensions
#' @param k The number of PC/EOF modes to retain
#' @param scale Logical, whether to scale before PCA
#' @param rotate Logical, whether to apply Varimax rotation. Rotation follows
#'   the standard REOF convention (Hannachi et al. 2007): varimax operates on
#'   sqrt(eigenvalue)-scaled EOFs, the stored patterns are the rotated
#'   loadings (unit norm, not mutually orthogonal), and amplitudes remain
#'   uncorrelated with sd = sqrt(rotated eigenvalue)
#' @param monthly Logical, whether to use monthly climatology
#' @param weight Logical, whether to apply area weighting
#' @param irlba_threshold Minimum number of data elements to trigger IRLBA usage
#'   (default: 50000). Set to Inf to always use base prcomp().
#'
#' @return A `patterns` object containing EOFs, amplitudes, and metadata
#' @export
patterns <- function(dat, k = 4, scale = FALSE, rotate = FALSE, monthly = FALSE, weight = TRUE, irlba_threshold = 500000){

  # Input validation
  check_stars_object(dat)

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
    valid_pixels = eofs$valid_pixels
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
  # Kaiser-normalise eigenvectors before rotation (Hannachi et al. 2007)
  scaled_loadings <- loadings_matrix %*% diag(sdev_vector)
  rot <- varimax(scaled_loadings)

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
  dims <- dim(dat)

  pc_names <- names0(k, 'PC')

  # Handle spatial dimensions differently for geometry vs raster
  flattened <- flatten_time_space(units::drop_units(dat)[1])
  anomaly_matrix_full <- flattened$matrix
  n_pixels <- ncol(anomaly_matrix_full)

  # Get valid pixel indices (non-NA)
  valid_pixels <- which(!apply(anomaly_matrix_full, 2, anyNA))

  # Validate k value
  max_k <- min(length(times) - 1, length(valid_pixels))
  check_k_valid(k, max_k)

  # Extract matrix for valid pixels (time x space)
  anomaly_matrix <- anomaly_matrix_full[, valid_pixels, drop = FALSE]

  # Apply spatial weights column-wise if provided
  if (!is.null(weights)) {
    if (length(weights) != n_pixels) {
      cli::cli_abort(
        "Length of {.arg weights} ({length(weights)}) must match number of spatial points ({n_pixels}).",
        class = "tidyeof_weight_mismatch"
      )
    }
    weights_valid <- weights[valid_pixels]
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

  if (has_geometry_dimension(dat)) {
    # For geometry: (geometry, time) → (geometry, PC)
    template <- dat[,,1:k, drop = FALSE]  # geometry, first k time slices
    spatial_patterns <- template %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")

    # Fill with EOF pattern data (geometry x PC)
    spatial_patterns[[1]] <- t(full_patterns)  # transpose to get geometry x PC

  } else {
    # For raster: (x, y, time) → (x, y, PC)
    template <- dat[,,,1:k, drop = FALSE]  # x, y, first k time slices
    spatial_patterns <- template %>%
      stars::st_set_dimensions('time', values = pc_names, names = 'PC') %>%
      setNames("weight")

    # For raster, reshape back to x,y,PC structure
    pattern_array <- array(full_patterns, dim = c(k, dims[[1]], dims[[2]]))
    reordered_patterns <- aperm(pattern_array, c(2, 3, 1))  # x, y, PC

    # Fill with EOF pattern data (x, y, PC)
    spatial_patterns[[1]] <- reordered_patterns
  }

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
           hi = eigenvalues * (1 + error) * 100 / total_var,
           cumvar_line = hi + 0.02 * max(hi))

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
  sqrt(area_values / mean(area_values, na.rm = TRUE))
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

#' Smart PCA Selection with Optional IRLBA Support
#'
#' Automatically selects between base `prcomp()` and `prcomp_irlba()` based on
#' data size and package availability. For large datasets, IRLBA provides
#' significant computational savings when only the first few components are needed.
#'
#' @param x A numeric matrix for PCA computation
#' @param k Number of components to compute
#' @param center Logical, whether to center the data
#' @param scale. Logical, whether to scale the data
#' @param size_threshold Minimum number of elements to trigger IRLBA (default: 50000)
#' @param ... Additional arguments passed to the PCA function
#'
#' @return A PCA result object compatible with `prcomp()` output
#' @keywords internal
perform_pca_smart <- function(x, k = NULL, center = TRUE, scale. = FALSE,
                              size_threshold, ...) {

  # Calculate data size
  data_size <- nrow(x) * ncol(x)
  use_irlba <- FALSE

  # Check if we should use IRLBA
  if (data_size >= size_threshold && !is.null(k)) {
    # Check if irlba package is available
    if (requireNamespace("irlba", quietly = TRUE)) {
      use_irlba <- TRUE
      cli::cli_inform(
        "Using IRLBA for efficient PCA computation on large dataset ({format(data_size, big.mark = ',')} elements)."
      )
    } else {
      cli::cli_inform(
        "Large dataset detected ({format(data_size, big.mark = ',')} elements) but package {.pkg irlba} is unavailable. Consider installing it for faster computation."
      )
    }
  }

  # Perform PCA using selected method
  if (use_irlba) {
    # Use IRLBA for efficient computation of first k components
      irlba::prcomp_irlba(x, n = k, center = center, scale. = scale., ...)
  } else {
    # Use base prcomp
    prcomp(x, center = center, scale. = scale., ...)
  }
}
