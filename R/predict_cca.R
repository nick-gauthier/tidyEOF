#' Title
#'
#' @param preds
#' @param obs
#' @param newdata
#' @param k
#'
#' @return
#' @export
#'
#' @examples
predict_cca <- function(preds, obs, newdata, k) { # is weight used anywhere?
  #should check if both preds and obs are both scaled/rotated and fail if not
  obs_amps <- obs$amplitudes %>%
    select(-time) %>%
    as.matrix()
  preds_amps <- preds$amplitudes %>%
    select(-time) %>%
    as.matrix()

  new_times <- st_get_dimension_values(newdata, 'time')

  # train CCA on training data
  cca <- cancor(preds_amps, obs_amps, xcenter = FALSE, ycenter = FALSE)

  # add a check here if k <= min(k_preds, k_obs)

    # get PCs from newdata
  project_patterns(preds, newdata) %>%
    tibble::column_to_rownames('time') %>%
    # predict CCA on newdata
    as.matrix() %*%
    cca$xcoef[,1:k, drop = FALSE] %*% # drop = FALSE prevents a vector returning when k = 1
    diag(cca$cor, nrow = k) %*%  # nrow likewise prevents weird behavior when k = 1
    solve(cca$ycoef)[1:k,,drop = FALSE] %>%
    as_tibble() %>%
    mutate(time = new_times, .before = 1) %>%
    reconstruct_field(obs, amplitudes = ., nonneg = FALSE) # better nonneg handling needed!
}

#' @export
project_patterns <- function(patterns, newdata) {
  new_times <- st_get_dimension_values(newdata, 'time')

  if(patterns$weight) newdata <- newdata * lat_weights(newdata) # weight by sqrt cosine latitude, in radians

  newdata %>%
    get_anomalies(patterns$climatology, scale = patterns$scaled, monthly = patterns$monthly) %>%
    units::drop_units() %>%
    split('time') %>% # split along the time dimension # this fails when folds = 1
    setNames(new_times) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>%
    predict(patterns$pca, .) %>%
    scale(center = FALSE, scale = patterns$pca$sdev) %>% # this should be from the training pcs
    .[,1:patterns$k, drop = FALSE] %>%
    {if(is.matrix(patterns$rotation)) . %*% patterns$rotation else .} %>%
    as_tibble() %>%
    setNames(names0(patterns$k, 'PC')) %>%
    mutate(time = new_times, .before = 1) # add times this way rather than with rownames above to preserve POSIXct class
}

predict_eot <- function(cv, k) {
  pmap(list(cv$eot, cv$test, cv$mean), ~predict(..1, ..2, n = k) + ..3) %>%
    brick() %>%
    stars::st_as_stars() %>%
    stars::st_set_dimensions('band', names = 'time') %>%
    transmute(value = units::set_units(layer.1.1, mm))
}


# special predict method for cv that scales predictors (should fold into original predict method in future)
predict_gam <- function(preds, obs, newdata, k) {

  k_preds <- preds$amplitudes %>%
    select(-time) %>%
    ncol()

  new_times <-  st_get_dimension_values(newdata, 'time')

  mod <- prep_data(preds, obs) %>% # in the cca code this is done internally, but here it calls another function . . . choose one?
    fit_model()

  new_pcs <- newdata %>%
    get_anomalies(preds$climatology, scale = FALSE, monthly = preds$monthly) %>% # should this scale too?
    units::drop_units() %>%
    area_weight() %>% # weight by sqrt cosine latitude, in radians
    split('time') %>% # split along the time dimension
    setNames(new_times) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # this should be from the training pcs
    .[,1:k_preds, drop = FALSE] %>%
    as_tibble(rownames = 'time')


  purrr::map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(time, PC, amplitude) %>%
    mutate(time = as.numeric(time)) %>%
    tidyr::pivot_wider(names_from = PC, names_prefix = 'PC', values_from = amplitude) %>% # names_prefix here because the gam predictor doesn't do pcs with pc attached, only the number. could fix upstream
    reconstruct_field(obs, .) # use train_obs to get training pcs
}



# same as above but for pcr
predict_pcr <- function(preds, obs, newdata, k) {
  # k is just here as a placeholder so the pmap command in fit_cv will run, but it doesn't have to be and should be removed
  k_preds <- preds$amplitudes %>%
    select(-time) %>%
    ncol()

  new_times <-  st_get_dimension_values(newdata, 'time')

  mod <- prep_data(preds, obs) %>%
    fit_pcr()

  new_pcs <- newdata %>%
    units::drop_units() %>%
    get_anomalies(preds$clmatology, monthly = preds$monthly) %>%
    area_weight() %>% # weight by sqrt cosine latitude, in radians
    split('time') %>% # split along the time dimension
    setNames(new_times) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # this should be from the training pcs
    .[,1:k_preds, drop = FALSE] %>%
    as_tibble(rownames = 'time')

  purrr::map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(time, PC, amplitude) %>%
    mutate(time = as.numeric(time)) %>%
    tidyr::pivot_wider(names_from = PC, names_prefix = 'PC', values_from = amplitude) %>% # names_prefix here because the gam predictor doesn't do pcs with pc attached, only the number. could fix upstream
    reconstruct_field(obs, .) # use train_obs to get training pcs
}
