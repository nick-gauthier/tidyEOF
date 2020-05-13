fit_kfold <- function(k_preds, k_obs, preds, obs){
  # this should be refactored ...

  # find the years of overlap between the response and predictor fields
  years <- intersect(preds$year, obs$year)
  preds <- filter(preds, year %in% years)
  obs <- filter(obs, year %in% years)

  # divide years into 5 contiguous folds
  n <- length(years)
  r <- n %% 5
  fold_times <- rep(n %/% 5, 5)
  fold_times[1:r] <- fold_times[1:r] + 1
  folds <- tibble(year = years, fold = rep(1:5, times = fold_times))

  # preprocess the training data for each fold
  obs_train <- map(1:5, ~ filter(folds, fold == .) %>%
                     anti_join(obs, ., by = 'year') %>% # training years
                     filter(year %in% years) %>% # remove years outside overlapping interval
                     get_patterns(k = k_obs))

  pred_train <- map(1:5, ~ filter(folds, fold == .) %>%
                      anti_join(preds, ., by = 'year') %>% # training years
                      get_patterns(k = k_preds))

  # preprocess test data for each fold
  test <- map(1:5, ~ filter(folds, fold == .) %>%
                semi_join(preds, ., by = 'year')) # test years

  # fit model to training data and use to predict new fields
  pmap_dfr(list(pred_train, obs_train, test), fit_cv)
}

# special predict method for cv that scales predictors (should fold into original predict method in future)
fit_cv <- function(pred_train, obs_train, test) {

  mod <- prep_data(pred_train, obs_train) %>%
    fit_model()

  preds <- left_join(test, pred_train$climatology, by = c("x", "y")) %>%
    mutate(SWE = SWE - swe_mean) %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>%
    dplyr::select(-swe_mean, -swe_sd) %>%
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    predict(pred_train$pca, .) %>%
    scale() %>%
    as_tibble(rownames = 'year')

  map(mod$mod, ~add_predictions(preds, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(year, PC, amplitude) %>%
    mutate(year = as.numeric(year)) %>%
    reconstruct_field(obs_train, .)  # use obs_train to get training pcs
}

get_errors <- function(x) {
  x %>%
    inner_join(prism_dat, by = c('x', 'y', 'year'), suffix = c('_recon', '_obs')) %>%
    mutate(SWE_recon = if_else(SWE_recon < 0.03, 0.03, SWE_recon), #for the metrics that use logq
           SWE_obs = if_else(SWE_obs < 0.03, 0.03, SWE_obs),
           error = SWE_recon - SWE_obs,
           relative_error = error / SWE_obs,
           accuracy_ratio = SWE_recon / SWE_obs,
           log_q = log(accuracy_ratio))
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
