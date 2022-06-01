#' Title
#'
#' @param times
#' @param kfolds
#'
#' @return
#' @export
#'
#' @examples
prep_folds <- function(times, kfolds = 5){
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

#' @export
prep_cca <- function(preds, obs, k_preds, k_obs, kfolds = 5, scale = FALSE, rotate = FALSE) {
  # find the years of overlap between the response and predictor fields
  time_steps <- intersect(st_get_dimension_values(preds, 'time'),
                          st_get_dimension_values(obs, 'time'))
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps, kfolds = kfolds)

  # preprocess the training data for each fold
  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .))  %>%
                            get_patterns(k = k_obs, scale = scale, rotate = rotate))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                              get_patterns(k = k_preds, scale = scale, rotate = rotate))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_obs, train_preds, test)
}


#' @export
prep_eot <- function(preds, obs, k_preds, k_obs, k_cca){
  # find the years of overlap between the response and predictor fields
  time_steps <- intersect(st_get_dimension_values(preds, 'time'),
                          st_get_dimension_values(obs, 'time'))
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps) # this does 5 fold by default, but could change

  # preprocess the training data for each fold
  obs_train_rast <- purrr::map(folds, ~ filter(obs, !(time %in% .)) %>%
                                 as('Raster'))

  pred_train_rast <-  purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                                   as('Raster'))

  obs_train <- purrr::map(obs_train_rast, ~ anomalize(.) %>%
                            denoise(k = k_obs, weighted = TRUE, verbose = FALSE))

  pred_train <-  purrr::map(pred_train_rast, ~ anomalize(.) %>%
                              denoise(k = k_preds, weighted = TRUE, verbose = FALSE))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .) %>%
                       as('Raster')) %>% # test years
    map2(pred_train_rast, ~.x - mean(.y)) # subtract the training predictor mean from the test predictors

  means <- purrr::map(obs_train_rast, raster::mean)

  # fit model to training data
  eots <- map2(pred_train, obs_train, ~eot(.x, .y, n = k_cca, standardised = FALSE, verbose = FALSE))

  list(eot = eots, # this hard codes 10 patterns, but could be changed to min(k_preds, k_obs)
       test = test,
       mean = means)
}

#' @export
prep_delta <- function(preds, obs) {
  # find the years of overlap between the response and predictor fields
  time_steps <- intersect(st_get_dimension_values(preds, 'time'),
                          st_get_dimension_values(obs, 'time'))
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps) # this does 5 fold by default, but could change

  # preprocess the training data for each fold
  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .)))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_preds, train_obs, test)
}