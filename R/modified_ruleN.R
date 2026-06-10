#' Test EOF significance using modified Rule N
#'
#' Tests whether the k-th eigenvalue is significantly different from noise
#' using a modified Rule N approach based on the Tracy-Widom distribution.
#'
#' Follows the Rule N significance framework for EOFs of Overland &
#' Preisendorfer (1982), comparing each eigenvalue against the noise null given
#' by the Tracy-Widom Type-1 law for the largest eigenvalue of a white-noise
#' covariance matrix (Johnstone 2001). The Tracy-Widom CDF is evaluated through
#' the gamma approximation of Chiani (2014); its constants (shape = 46.4,
#' scale factor = 0.186, location = 9.85) come from fitting a gamma CDF to the
#' Tracy-Widom Type-1 distribution. Because the null assumes independent noise,
#' the test tends to be liberal for spatially correlated geophysical fields.
#'
#' @references
#' Overland, J.E. & Preisendorfer, R.W. (1982). A significance test for
#' principal components applied to a cyclone climatology. \emph{Monthly Weather
#' Review}, 110(1), 1-4.
#'
#' Johnstone, I.M. (2001). On the distribution of the largest eigenvalue in
#' principal components analysis. \emph{Annals of Statistics}, 29(2), 295-327.
#'
#' Chiani, M. (2014). Distribution of the largest eigenvalue for real Wishart
#' and Gaussian random matrices and a simple approximation for the Tracy-Widom
#' distribution. \emph{Journal of Multivariate Analysis}, 129, 69-81.
#'
#' @param lambdas Vector of eigenvalues from PCA
#' @param k Index of eigenvalue to test
#' @param M Number of spatial points (grid cells)
#' @param n Number of time steps
#' @param p Significance level (default 0.05)
#' @param total_var Total variance of the data (sum of all eigenvalues).
#'   Required when `lambdas` contains only the leading modes (e.g., from
#'   IRLBA); the remaining noise variance is then `total_var` minus the
#'   eigenvalues above `k`. If NULL (default), all eigenvalues must be present.
#'
#' @return Logical, TRUE if eigenvalue is significant at level p
#' @export
eigen_test <- function(lambdas, k, M, n, p = 0.05, total_var = NULL){
  nrank <- min(n - 1, M)

  if (is.null(total_var) && length(lambdas) < nrank) {
    cli::cli_abort(
      "{.arg lambdas} has only {length(lambdas)} of {nrank} eigenvalues. Supply {.arg total_var} when the spectrum is truncated (e.g., IRLBA).",
      class = "tidyeof_truncated_spectrum"
    )
  }

  kstar <- M - k + 1
  nk <- n - k + 1
  mu <- (sqrt(nk - 0.5) + sqrt(kstar - 0.5))^2
  sigma <- sqrt(mu) * (1 / sqrt(nk - 0.5) + 1 / sqrt(kstar - 0.5)) ^ (1/3)
  shape <- 46.4
  beta <- (0.186 * sigma) / max(nk, kstar)
  zeta <- (mu - 9.85 * sigma) / max(nk, kstar)
  noise_sum <- if (is.null(total_var)) {
    sum(lambdas[k:nrank])
  } else {
    total_var - sum(lambdas[seq_len(k - 1)])
  }
  lambda_star <- lambdas[k] / (noise_sum / (nrank - k + 1))

  (1 - pgamma(((lambda_star - zeta) / beta), shape)) < p
}

#' Find Rule N significance cutoff for a patterns object
#'
#' Tests each eigenvalue in sequence and returns the index of the last
#' significant mode. This is the maximum k supported by the data according
#' to the modified Rule N test.
#'
#' @param x A patterns object
#' @param p Significance level (default 0.05)
#'
#' @return Integer: index of last significant eigenvalue, or 0 if none are significant
#' @keywords internal
rule_n_cutoff <- function(x, p = 0.05) {
  lambdas <- x$eigenvalues$eigenvalues
  n_times <- nrow(x$amplitudes)
  n_valid <- length(x$valid_pixels)
  max_test <- min(length(lambdas), n_times - 1)

  significant <- vapply(seq_len(max_test), function(k) {
    eigen_test(lambdas, k = k, M = n_valid, n = n_times, p = p,
               total_var = x$total_variance)
  }, logical(1))

  # Last TRUE in contiguous sequence from the start
  if (!significant[1]) return(0L)
  # Find where significance first drops off
  first_nonsig <- which(!significant)[1]
  if (is.na(first_nonsig)) max_test else first_nonsig - 1L
}
