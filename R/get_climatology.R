#' Extract calendar month numbers from a time vector
#'
#' Locale-independent month extraction: `"%m"` is numeric, unlike `"%B"`, and
#' `format()` dispatches on the time class, so `Date`, `POSIXct` (respecting
#' its tzone attribute), and CF-calendar classes with format methods all work.
#'
#' @param times A vector of time values
#' @return Integer vector of calendar months (1-12)
#' @keywords internal
month_index <- function(times) {
  as.integer(format(times, "%m"))
}

#' Flatten a stars object to a (dim x space) numeric matrix
#'
#' Moves `dim_name` to the first dimension and flattens the remaining
#' (spatial) dimensions. For multi-attribute objects, attributes are
#' concatenated variable-major — all spatial columns for the first attribute,
#' then all columns for the second, etc. — matching the column ordering of
#' [flatten_time_space()].
#'
#' Units are dropped: the result is always a plain `double` matrix with no
#' `"units"` class. No block map is returned; callers that need per-variable
#' column ranges should derive them from the spatial shape or use
#' [flatten_time_space()].
#'
#' @param x A stars object (one or more attributes)
#' @param dim_name Name of the dimension to keep as rows
#' @return A plain numeric matrix with rows = `dim_name` and columns =
#'   flattened space (concatenated variable-major for multi-attribute objects).
#'   Units are stripped; no block map is attached.
#' @keywords internal
flatten_dim_space <- function(x, dim_name) {
  spatial <- setdiff(names(stars::st_dimensions(x)), dim_name)
  permuted <- aperm(units::drop_units(x), c(dim_name, spatial))
  blocks <- purrr::map(seq_along(names(x)), function(i) {
    arr <- permuted[[i]]
    matrix(arr, nrow = dim(arr)[1])
  })
  do.call(cbind, blocks)
}

#' Build a calendar month dimension (values 1:12)
#' @keywords internal
month_dimension <- function() {
  month_dim <- list(
    from = 1L,
    to = 12L,
    offset = NA_real_,
    delta = NA_real_,
    refsys = NA_character_,
    point = FALSE,
    values = 1:12
  )
  class(month_dim) <- "dimension"
  month_dim
}

#' Calculate climatological mean and standard deviation for spatial data
#'
#' Computes climatological statistics (mean and standard deviation) for a spatial field,
#' either annually or monthly. Preserves spatial dimensions and units from the input data.
#'
#' @param dat A stars object containing a spatial field with dimensions (x, y, time)
#' @param monthly Logical. If TRUE, computes monthly climatology. If FALSE (default),
#'   computes statistics over the entire period.
#'
#' @return A list with two stars objects:
#'   \item{mean}{Climatological mean with original spatial dimensions and units}
#'   \item{sd}{Climatological standard deviation with same structure}
#'
#'   For monthly climatologies, the `month` dimension always has values 1:12 in
#'   calendar order, regardless of which month the data starts in. Months not
#'   present in the data are NA. Complete years are not required, but a message
#'   is emitted when months have unequal sample sizes.
#'
#' @examples
#' # Create sample data
#' library(stars)
#' times <- seq(as.Date("2000-01-01"), as.Date("2002-12-31"), by = "month")
#' x <- seq(0, 1, length.out = 10)
#' y <- seq(0, 1, length.out = 10)
#' dat <- stars::st_as_stars(array(rnorm(10*10*36), c(10, 10, 36))) %>%
#'   st_set_dimensions(1, values = x, name = "x") %>%
#'   st_set_dimensions(2, values = y, name = "y") %>%
#'   st_set_dimensions(3, values = times, name = "time")
#'
#' # Calculate annual climatology
#' clim <- get_climatology(dat)
#'
#' # Calculate monthly climatology
#' monthly_clim <- get_climatology(dat, monthly = TRUE)
#'
#' @export
get_climatology <- function(dat, monthly = FALSE) {
  # Basic validation
  if (!inherits(dat, "stars")) {
    rlang::abort("Input must be a stars object", class = "tidyeof_invalid_input")
  }

  if (monthly) {
    times <- st_get_dimension_values(dat, 'time')
    m <- month_index(times)

    counts <- tabulate(m, nbins = 12L)
    if (length(unique(counts)) > 1) {
      cli::cli_inform(
        "Months have unequal sample sizes ({min(counts)}-{max(counts)} time steps per month); the climatology will be noisier for some months.",
        class = "tidyeof_unbalanced_months"
      )
    }

    # One path for raster and geometry: flatten to (time x space), compute
    # per-calendar-month statistics, and rebuild with a month dimension that
    # is always 1:12 so month number doubles as the array index.
    spatial <- setdiff(names(st_dimensions(dat)), "time")
    mat <- flatten_dim_space(dat[1], "time")
    n_space <- ncol(mat)

    month_stat <- function(stat_fn) {
      vapply(1:12, function(mm) {
        idx <- which(m == mm)
        if (length(idx) == 0) {
          rep(NA_real_, n_space)
        } else {
          stat_fn(mat[idx, , drop = FALSE])
        }
      }, numeric(n_space))
    }

    mean_mat <- month_stat(function(x) colMeans(x, na.rm = TRUE))
    sd_mat <- month_stat(function(x) apply(x, 2, sd, na.rm = TRUE))

    new_dims <- st_dimensions(dat)[spatial]
    new_dims$month <- month_dimension()
    class(new_dims) <- "dimensions"
    spatial_shape <- dim(dat)[spatial]

    to_stars <- function(values) {
      stars::st_as_stars(
        array(values, dim = c(spatial_shape, month = 12L)),
        dimensions = new_dims
      ) %>%
        setNames(names(dat)[1])
    }

    mean_result <- to_stars(mean_mat)
    sd_result <- to_stars(sd_mat)
  } else {
    # Annual climatology calculation
    spatial_dims <- get_spatial_dimensions(dat)
    mean_result <- st_apply(dat, spatial_dims, mean, na.rm = TRUE, rename = FALSE)
    sd_result <- st_apply(dat, spatial_dims, sd, na.rm = TRUE, rename = FALSE)
  }

  list(
    mean = restore_units(mean_result, dat),
    sd = restore_units(sd_result, dat)
  )
}

#' Apply or remove a monthly climatology by calendar month index
#'
#' Shared engine for [get_anomalies()] and [restore_climatology()] with
#' `monthly = TRUE`. Because the climatology's month dimension is always 1:12
#' in calendar order, each time step's climatology is looked up by direct
#' indexing with its month number -- no ordering convention to maintain, no
#' complete-years requirement, and identical code for raster and geometry data.
#'
#' @param dat A stars object with a time dimension (data or anomalies)
#' @param clim Climatology list from [get_climatology()] with `monthly = TRUE`
#' @param scale Logical, whether to divide/multiply by the climatological sd
#' @param direction "anomalize" (subtract climatology) or "restore" (add it back)
#' @return A stars object with the same structure as `dat`, units dropped
#' @keywords internal
apply_monthly_climatology <- function(dat, clim, scale,
                                      direction = c("anomalize", "restore")) {
  direction <- match.arg(direction)

  clim_dims <- names(st_dimensions(clim$mean))
  if (!"month" %in% clim_dims ||
      !identical(as.integer(st_get_dimension_values(clim$mean, "month")), 1:12)) {
    cli::cli_abort(
      c(
        "Monthly climatology must have a {.field month} dimension with values 1:12.",
        "i" = "Regenerate the climatology with {.fn get_climatology}."
      ),
      class = "tidyeof_invalid_climatology"
    )
  }

  times <- st_get_dimension_values(dat, 'time')
  m <- month_index(times)

  mn_mat <- flatten_dim_space(clim$mean, "month")
  sd_mat <- if (scale) flatten_dim_space(clim$sd, "month") else NULL

  spatial <- setdiff(names(st_dimensions(dat)), "time")
  permuted <- aperm(units::drop_units(dat), c("time", spatial))
  arr <- permuted[[1]]
  mat <- matrix(arr, nrow = length(times))

  if (ncol(mat) != ncol(mn_mat)) {
    cli::cli_abort(
      "Spatial size mismatch: data has {ncol(mat)} cells but the climatology has {ncol(mn_mat)}.",
      class = "tidyeof_grid_mismatch"
    )
  }

  available <- which(rowSums(!is.na(mn_mat)) > 0)
  missing_months <- setdiff(unique(m), available)
  if (length(missing_months) > 0) {
    cli::cli_abort(
      "Data contains month{?s} {month.name[sort(missing_months)]} not present in the climatology.",
      class = "tidyeof_missing_month"
    )
  }

  if (direction == "anomalize") {
    mat <- mat - mn_mat[m, , drop = FALSE]
    if (scale) mat <- mat / sd_mat[m, , drop = FALSE]
  } else {
    if (scale) mat <- mat * sd_mat[m, , drop = FALSE]
    mat <- mat + mn_mat[m, , drop = FALSE]
  }

  permuted[[1]] <- array(mat, dim = dim(arr))
  aperm(permuted, names(st_dimensions(dat)))
}

#' Calculate anomalies from a climatological mean
#'
#' @param dat A stars object with dimensions (x, y, time)
#' @param clim Optional climatology from get_climatology(). If NULL, computed internally
#' @param scale Logical. If TRUE, divide by standard deviation
#' @param monthly Logical. If TRUE, compute monthly anomalies. Each time step is
#'   matched to its calendar month, so the data may start in any month, span
#'   partial years, or cover a single year. Aborts if the data contains a month
#'   absent from the climatology.
#' @return A stars object with anomalies
#' @export
get_anomalies <- function(dat, clim = NULL, scale = FALSE, monthly = FALSE) {
  # Basic validation
  if (!inherits(dat, "stars")) {
    rlang::abort("Input must be a stars object", class = "tidyeof_invalid_input")
  }

  # Get or validate climatology
  if (is.null(clim)) {
    clim <- get_climatology(dat, monthly = monthly)
  } else if (!is.list(clim) || !all(c("mean", "sd") %in% names(clim))) {
    rlang::abort("climatology must be a list with 'mean' and 'sd' stars objects (from get_climatology())",
                 class = "tidyeof_invalid_input")
  }

  if (monthly) {
    out <- apply_monthly_climatology(dat, clim, scale = scale, direction = "anomalize")
    # Scaled anomalies are dimensionless; unscaled keep the data's units
    if (!scale) out <- restore_units(out, dat)
    return(out)
  }

  out <- dat - clim$mean
  if (scale) {
    out <- out / clim$sd
  }
  out
}


#' Restore original field from anomalies and climatology
#'
#' Reverses the operation of \code{get_anomalies()}, adding the climatological
#' mean (and optionally multiplying by standard deviation) back to anomaly fields.
#'
#' @param anomalies A stars object containing anomalies (from \code{get_anomalies()})
#' @param clim A climatology list with \code{mean} and \code{sd} stars objects
#'   (from \code{get_climatology()})
#' @param scale Logical. If TRUE, multiply by standard deviation before adding mean
#'   (use when anomalies were standardized)
#' @param monthly Logical. If TRUE, restore using monthly climatology. Time
#'   steps are matched to the climatology by calendar month, so the anomalies
#'   may start in any month or span partial years.
#'
#' @return A stars object with the original field restored
#'
#' @examples
#' \dontrun{
#' clim <- get_climatology(dat)
#' anom <- get_anomalies(dat, clim)
#' restored <- restore_climatology(anom, clim)
#' }
#'
#' @export
restore_climatology <- function(anomalies, clim, scale = FALSE, monthly = FALSE) {
  # Basic validation
  if (!inherits(anomalies, "stars")) {
    rlang::abort("Anomalies must be a stars object", class = "tidyeof_invalid_input")
  }
  if (!is.list(clim) || !all(c("mean", "sd") %in% names(clim))) {
    rlang::abort("Climatology must be a list with 'mean' and 'sd' stars objects (from get_climatology())",
                 class = "tidyeof_invalid_input")
  }

  if (monthly) {
    out <- apply_monthly_climatology(anomalies, clim, scale = scale, direction = "restore")
    return(restore_units(out, clim$mean))
  }

  target_mean <- clim$mean
  target_sd <- clim$sd

  # Check if anomalies and climatology have matching units
  anomalies_has_units <- any(sapply(names(anomalies), function(var) !is.null(tryCatch(units(anomalies[[var]]), error = function(e) NULL))))
  clim_has_units <- any(sapply(names(target_mean), function(var) !is.null(tryCatch(units(target_mean[[var]]), error = function(e) NULL))))

  # If units mismatch, drop units from climatology to match anomalies
  if(!anomalies_has_units && clim_has_units) {
    target_mean <- units::drop_units(target_mean)
    target_sd <- units::drop_units(target_sd)
  }

  # Restore climatology
  if (scale) {
    anomalies <- anomalies * target_sd
  }
  out <- anomalies + target_mean

  # Restore units
  restore_units(out, clim$mean)
}


# Restore units from reference stars object to new stars object
restore_units <- function(new, ref) {
  for (var in names(new)) {
    ref_units <- tryCatch(units(ref[[var]]), error = function(e) NULL)
    if (!is.null(ref_units)) {
      new[[var]] <- units::set_units(new[[var]], ref_units, mode = 'standard')
    }
  }
  new
}
