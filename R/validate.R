fit_kfold <- function(k_preds, k_obs){
  preds <- exp_test[[k_preds, "cera_patterns"]][[1]]
  obs <- exp_test[[k_obs, "prism_patterns"]][[1]]

  # prepare data frame for cross validation by dividing into k folds
  folds <- prep_data(preds, obs) %>%
    group_by(PC) %>%
    mutate(fold = rep(1:5, each = 6, length.out = n())) %>%
    ungroup()

  train <- map(1:5, ~filter(folds, fold != .) %>% dplyr::select(-fold))
  test <- map(1:5, ~filter(folds, fold == .) %>% dplyr::select(-fold, -PC, -amplitude) %>% unique)

  map2_dfr(train, test, ~ fit_model(.x) %>%
             predict_patterns(.y) %>%
             reconstruct_field(obs, .))
}

get_errors <- function(x) {
  x %>%
    left_join(prism_dat, by = c('x', 'y', 'year'), suffix = c('_recon', '_obs')) %>%
    mutate(SWE_recon = if_else(SWE_recon < 0.01, 0.01, SWE_recon),
           SWE_obs = if_else(SWE_obs < 0.01, 0.01, SWE_obs),
           error = SWE_recon - SWE_obs,
           relative_error = error / SWE_obs,
           accuracy_ratio = SWE_recon / SWE_obs,
           log_q = log(accuracy_ratio))
}

get_scores <- function(x) {
  x%>%
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
    ungroup()
}