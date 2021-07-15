prep_folds <- function(times, kfolds = 5){
  # divide years into kfolds contiguous folds
  n <- length(times)
  r <- n %% kfolds
  fold_times <- rep(n %/% kfolds, kfolds)
  fold_times[1:r] <- fold_times[1:r] + 1
  tibble(time = times, fold = rep(1:kfolds, times = fold_times)) %>%
    group_nest(fold) %>%
    pull(data) %>%
    purrr::map(pull)
}

prep_cca <- function(preds, obs, k_preds, k_obs) {
  # find the years of overlap between the response and predictor fields
  time_steps <- intersect(st_get_dimension_values(preds, 'time'),
                          st_get_dimension_values(obs, 'time'))
  preds <- filter(preds, time %in% time_steps)
  obs <- filter(obs, time %in% time_steps)
  folds <- prep_folds(time_steps) # this does 5 fold by default, but could change

  # preprocess the training data for each fold
  train_obs <- purrr::map(folds, ~ filter(obs, !(time %in% .))  %>%
                     get_patterns(k = k_obs))

  train_preds <- purrr::map(folds, ~ filter(preds, !(time %in% .)) %>%
                      get_patterns(k = k_preds))

  # preprocess test data for each fold
  test <- purrr::map(folds, ~ filter(preds, time %in% .))

  tibble(train_obs, train_preds, test)
}


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

  means <- purrr::map(obs_train_rast, mean)

  # fit model to training data
  eots <- map2(pred_train, obs_train, ~eot(.x, .y, n = k_cca, standardised = FALSE, verbose = FALSE))

  list(eot = eots, # this hard codes 10 patterns, but could be changed to min(k_preds, k_obs)
       test = test,
       mean = means)
}

predict_eot <- function(cv, k) {
  pmap(list(cv$eot, cv$test, cv$mean), ~predict(..1, ..2, n = k) + ..3) %>%
    brick() %>%
    st_as_stars() %>%
    st_set_dimensions('band', names = 'time') %>%
    transmute(value = units::set_units(layer.1.1, mm))
}

cv_eot_error <- function(recon, obs) {
  # the dates shouldn't be hardcoded!
  error <- recon - filter(obs, between(time, 1982, 2010))

  rmse <- sqrt(mean(pull(error, 1) ^ 2, na.rm = TRUE)) %>% as.numeric()

  corr <- obs %>%
    filter(between(time, 1982, 2010)) %>%
    get_total_num() %>%
    cor(get_total_num(recon))

  c(rmse = rmse,
    corr = corr)
}

get_total_num <- function(dat) {
  (dat * st_area(dat)) %>%
    st_apply(3, function(x) sum(x, na.rm = TRUE)) %>%
    pull() %>%
    as.numeric()
}

get_total <- function(dat) {
(dat * st_area(dat)) %>%
  st_apply(3, function(x) sum(x, na.rm = TRUE)) %>%
    as_tibble() %>%
    mutate(SWE = units::set_units(SWE, m^2*mm) %>%
             units::set_units(Tl))
}

fit_cv <- function(dat, fun, k, obs) {
  recon <- pmap(list(dat$train_preds, dat$train_obs, dat$test), fun, k = k) %>%
    do.call('c', .)

  error <- dplyr::filter(recon, between(time, 1982, 2010)) - dplyr::filter(obs, between(time, 1982, 2010))

  rmse <- sqrt(mean(pull(error, 1) ^ 2, na.rm = TRUE)) %>% as.numeric()

  corr <- obs %>%
    filter(between(time, 1982, 2010)) %>%
    get_total_num() %>%
    cor(get_total_num(filter(recon, between(time, 1982, 2010))))

  c(rmse = rmse,
    corr = corr)
}

total_swe_corr <- function(errors) {
  errors %>%
    left_join(areas, by = c("x", "y")) %>%
    group_by(year) %>%
    summarise(SWE_obs = sum(SWE_obs * area, na.rm = TRUE),
              SWE_recon = sum(SWE_recon * area, na.rm = TRUE)) %>%
    summarise(correlation = cor(SWE_recon, SWE_obs)) %>%
    pull(correlation)
}

#' @export
predict_cca <- function(preds, obs, newdata, k) {
  #check if both preds and obs or scaled or not?
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

# check that the newdata argument is handled the same as not
# identical(cca_fit(prism_dat, cera_dat, k_preds = 4, k_obs = 4) ,
#          cca_fit(prism_dat, cera_dat, cera_dat, k_preds = 4, k_obs = 4) %>% filter(year >= 1982) )

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

#####

get_errors <- function(x) {
  x %>%
    inner_join(prism_dat, by = c('x', 'y', 'year'), suffix = c('_recon', '_obs')) %>%
    mutate(error = SWE_recon - SWE_obs)
        #   relative_error = error / SWE_obs,
         #  accuracy_ratio = if_else(SWE_recon < 0.03, 0.03, SWE_recon) /  if_else(SWE_obs < 0.03, 0.03, SWE_obs),
        #   log_q = log(accuracy_ratio))
}

get_scores <- function(x) {
  x %>%
    #filter(sd > 0) %>%
    summarise(xbar = mean(SWE_obs),
              sd = sd(SWE_obs),
              me = mean(error),
              mse = mean(error ^ 2),
              rmse = sqrt(mse),
              srmse = rmse/sd,
              mae = mean(abs(error)),
              mdae = median(abs(error)),
              mse_clim = sum((SWE_obs - xbar)^2) * (1 / (n() - 1)) * ((n() - 1) / n()),
              msss = 1 - mse / mse_clim,
              mape = mean(abs(relative_error)) * 100,
              mpe = mean(relative_error) * 100,
              mdsa = 100 * (exp(median(abs(log_q))) - 1), # median symmetric accuracy
              sspb = 100 * sign(median(log_q)) * (exp(abs(median(log_q))) - 1),
              rmsle = sqrt(mean(log(SWE_obs / SWE_recon)^2))) %>%
    ungroup() %>%
    dplyr::select(-xbar, -sd, -mse_clim)
}

# for calculating domain wide totals
get_areas <- function(dat, areas){
  dat %>%
    left_join(areas) %>%
    group_by(year) %>%
    summarise(SWE = sum(SWE * area))
}


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

delta <- function(pred, obs, newdata = NULL, k = NULL) { # k is just aplaceholder
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred)
  obs_clim <- get_climatology(obs)

  (newdata / pred_clim['mean']) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    `*`(obs_clim['mean'])
}

delta_add <- function(pred, obs, newdata = NULL, k = NULL) { # k is just aplaceholder
  if(is.null(newdata)) newdata <- pred
  pred_clim <- get_climatology(pred)
  obs_clim <- get_climatology(obs)

  # simplify units here!
  (units::drop_units(newdata) - pred_clim['mean']) %>%
    st_warp(slice(obs, 'time', 1), use_gdal = TRUE, method = 'bilinear') %>%
    setNames(names(newdata)) %>%
    mutate(across(everything(), ~units::set_units(.x, units(newdata[[1]]), mode = 'standard'))) %>%
    st_set_dimensions('band', values = st_get_dimension_values(newdata, 'time'), names = 'time') %>%
    units::drop_units() %>%
    `+`(obs_clim['mean']) %>%
    # make varname generic!
    mutate(SWE = units::set_units(if_else(SWE < 0, 0, SWE), mm))
}

