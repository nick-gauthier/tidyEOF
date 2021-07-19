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
predict_cca <- function(preds, obs, newdata, k) {
  #c should heck if both preds and obs are both scaled/rotated and fail if not
  obs_amps <- obs$amplitudes %>%
    select(-time) %>%
    as.matrix()
  preds_amps <- preds$amplitudes %>%
    select(-time) %>%
    as.matrix()

  new_times <-  st_get_dimension_values(newdata, 'time')

  k_preds <- ncol(preds_amps)

  # train CCA on training data
  cca <- cancor(preds_amps, obs_amps, xcenter = FALSE, ycenter = FALSE)

  # get PCs from newdata
  new_pcs <- newdata %>%
    units::drop_units() %>%
    `-`(preds$climatology['mean']) %>%
    {if(preds$scaled) . / preds$climatology['sd'] else .} %>%
    area_weight() %>% # weight by sqrt cosine latitude, in radians
    split('time') %>% # split along the time dimension
    setNames(new_times) %>%
    as_tibble() %>%
    dplyr::select(-c(x,y)) %>%
    na.omit() %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # this should be from the training pcs
    {if(!is.na(preds$rotation)) . %*% preds$rotation else .} %>%
    .[,1:k_preds, drop = FALSE]

  # add a check here if k <= min(k_preds, k_obs)

  # predict CCA on newdata
  new_pcs %*%
    cca$xcoef[,1:k, drop = FALSE] %*% # drop = FALSE prevents a vector returning when k = 1
    diag(cca$cor, nrow = k) %*%  # nrow likewise prevents weird behavior when k = 1
    solve(cca$ycoef)[1:k,,drop = FALSE] %>%
    `rownames<-`(new_times) %>%
    as_tibble(rownames = 'time') %>%
    mutate(time = as.numeric(time)) %>%
    reconstruct_field(obs, amplitudes = .)
}

predict_eot <- function(cv, k) {
  pmap(list(cv$eot, cv$test, cv$mean), ~predict(..1, ..2, n = k) + ..3) %>%
    brick() %>%
    st_as_stars() %>%
    st_set_dimensions('band', names = 'time') %>%
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
    units::drop_units() %>%
    `-`(preds$climatology['mean']) %>%
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


  map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
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
    `-`(preds$climatology['mean']) %>%
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

  map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(time, PC, amplitude) %>%
    mutate(time = as.numeric(time)) %>%
    tidyr::pivot_wider(names_from = PC, names_prefix = 'PC', values_from = amplitude) %>% # names_prefix here because the gam predictor doesn't do pcs with pc attached, only the number. could fix upstream
    reconstruct_field(obs, .) # use train_obs to get training pcs
}
