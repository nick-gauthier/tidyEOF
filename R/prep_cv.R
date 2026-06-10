#' Prepare contiguous cross-validation folds
#'
#' Divides time steps into k contiguous folds for temporal cross-validation.
#' Folds are kept contiguous to preserve temporal structure.
#'
#' @param times Vector of time values to split
#' @param kfolds Number of folds (default 5)
#'
#' @return List of k vectors, each containing the time values for that fold
#'
#' @seealso [prep_cv_folds()] for preparing complete CV folds with pre-computed patterns
#' @seealso [tune_cca()] for hyperparameter grid search
#' @export
prep_folds <- function(times, kfolds = 5){
  if (kfolds < 2) cli::cli_abort("{.arg kfolds} must be at least 2, not {kfolds}.")
  # divide years into kfolds contiguous folds
  n <- length(times)
  r <- n %% kfolds
  fold_times <- rep(n %/% kfolds, kfolds)
  if(r > 0) fold_times[1:r] <- fold_times[1:r] + 1
  tibble(time = times, fold = rep(1:kfolds, times = fold_times)) %>%
    group_nest(fold) %>%
    pull(data) %>%
    purrr::map(pull)
}
