#' Flatten a time-indexed stars object into a matrix
#'
#' Internal helper to turn a stars object (one or more attributes) with a
#' `time` dimension into a matrix of shape time x (V * space) along with the
#' metadata needed to reconstruct the original spatial layout.  When the object
#' has multiple attributes they are concatenated variable-major: all columns for
#' the first attribute, then all columns for the second, etc.  The returned
#' `block_map` is a named list mapping each attribute name to its column
#' indices, and `n_space` is the number of spatial cells per attribute.
#'
#' @param dat A stars object with a `time` dimension
#' @return A list containing the flattened matrix, block_map (named list of
#'   column index ranges per attribute), n_space (spatial cells per attribute),
#'   spatial dimension names, their shape, and coordinate values
#' @keywords internal
flatten_time_space <- function(dat) {
  check_stars_object(dat)

  dims <- stars::st_dimensions(dat)
  if (!"time" %in% names(dims)) {
    cli::cli_abort("Object must contain a {.field time} dimension to be flattened.")
  }

  spatial_dims <- setdiff(names(dims), "time")
  permuted <- aperm(dat, c("time", spatial_dims))

  var_names <- names(dat)
  blocks <- purrr::map(seq_along(var_names), function(i) {
    arr <- permuted[[i]]
    matrix(arr, nrow = dim(arr)[1], ncol = prod(dim(arr)[-1]))
  })
  mat <- do.call(cbind, blocks)

  n_space <- ncol(blocks[[1]])
  block_map <- setNames(
    purrr::map(seq_along(var_names), ~((.x - 1L) * n_space + 1L):(.x * n_space)),
    var_names
  )

  spatial_values <- purrr::map(spatial_dims, ~stars::st_get_dimension_values(permuted, .x))
  names(spatial_values) <- spatial_dims

  list(
    matrix = mat,
    block_map = block_map,
    n_space = n_space,
    spatial_dims = spatial_dims,
    spatial_shape = dim(permuted[[1]])[-1],
    spatial_values = spatial_values
  )
}

#' Convert a time-by-space matrix back to a stars object
#'
#' @param mat Matrix with rows = time, columns = flattened spatial cells
#' @param template_eofs EOF stars object providing spatial metadata
#' @param spatial_template Stars object supplying spatial dimension metadata
#'   (e.g., the climatology)
#' @param valid_pixels Integer indices of spatial cells with valid data
#' @param times Vector of time values
#' @param var_names Character vector of attribute names for the result
#' @keywords internal
matrix_to_spacetime <- function(mat,
                                template_eofs,
                                spatial_template,
                                valid_pixels,
                                times,
                                var_names) {
  dims <- stars::st_dimensions(template_eofs)
  if (!"PC" %in% names(dims)) {
    cli::cli_abort("Template EOFs must include a {.field PC} dimension.")
  }

  spatial_dims <- setdiff(names(dims), "PC")
  spatial_sizes <- purrr::map_int(spatial_dims, ~dimension_size(dims[[.x]]))
  total_space <- prod(spatial_sizes)

  full_mat <- matrix(NA_real_, nrow = nrow(mat), ncol = total_space)
  full_mat[, valid_pixels] <- mat

  out_array <- array(t(full_mat), dim = c(spatial_sizes, nrow(mat)))

  spatial_template_dims <- stars::st_dimensions(spatial_template)
  spatial_dim_names <- intersect(names(spatial_template_dims), spatial_dims)

  new_dims <- spatial_template_dims[spatial_dim_names]

  time_dim <- list(
    from = 1,
    to = nrow(mat),
    offset = NA_real_,
    delta = NA_real_,
    refsys = if (inherits(times, "POSIXct")) attr(times, "tzone") else if (inherits(times, "Date")) "Date" else NA_character_,
    point = FALSE,
    values = times
  )
  class(time_dim) <- "dimension"

  new_dims$time <- time_dim
  class(new_dims) <- "dimensions"

  out <- stars::st_as_stars(out_array, dimensions = new_dims)

  setNames(out, var_names)
}

#' Extract the EOF loading matrix (concatenated variable-space x PC) from a patterns object
#'
#' Stacks every variable's loadings into the concatenated variable-major
#' layout used by [flatten_time_space()], so rows align with
#' `patterns$valid_pixels` and `patterns$block_map`.
#' @param patterns A patterns object
#' @return A numeric matrix with prod(spatial) * n_vars rows and k columns
#' @keywords internal
eof_loading_matrix <- function(patterns) {
  do.call(rbind, purrr::map(names(patterns$eofs), function(v) {
    arr <- patterns$eofs[[v]]
    matrix(arr, nrow = prod(dim(arr)[-length(dim(arr))]), ncol = patterns$k)
  }))
}

#' Obtain the size of a stars dimension definition
#' @keywords internal
dimension_size <- function(dimension) {
  if (!is.null(dimension$values)) {
    length(dimension$values)
  } else if (!is.null(dimension$from) && !is.null(dimension$to)) {
    dimension$to - dimension$from + 1
  } else {
    cli::cli_abort("Unable to determine dimension length from metadata.")
  }
}
