prep_folds <- function(years, kfolds = 5){
  # divide years into 5 contiguous folds
  n <- length(years)
  r <- n %% kfolds
  fold_times <- rep(n %/% kfolds, kfolds)
  fold_times[1:r] <- fold_times[1:r] + 1
  tibble(year = years, fold = rep(1:kfolds, times = fold_times))
}

prep_cca <- function(preds, obs, k_preds, k_obs) {
  # find the years of overlap between the response and predictor fields
  years <- intersect(preds$year, obs$year)
  preds <- filter(preds, year %in% years)
  obs <- filter(obs, year %in% years)
  folds <- prep_folds(years)

  # preprocess the training data for each fold
  train_obs <- map(1:5, ~ filter(folds, fold == .) %>%
                     anti_join(obs, ., by = 'year') %>% # training years
                     get_patterns(k = k_obs))

  train_preds <- map(1:5, ~ filter(folds, fold == .) %>%
                      anti_join(preds, ., by = 'year') %>% # training years
                      get_patterns(k = k_preds))

  # preprocess test data for each fold
  test <- map(1:5, ~ filter(folds, fold == .) %>%
                semi_join(preds, ., by = 'year')) # test years

  tibble(train_obs, train_preds, test)
}


prep_eot <- function(preds, obs, k_preds, k_obs, k_cca){
  # find the years of overlap between the response and predictor fields
  years <- intersect(preds$year, obs$year)
  preds <- filter(preds, year %in% years)
  obs <- filter(obs, year %in% years)
  folds <- prep_folds(years)

  # preprocess the training data for each fold
  obs_train_rast <- purrr::map(1:5, ~ filter(folds, fold == .) %>%
                                 anti_join(obs, ., by = 'year') %>% # training years
                                 spread(year, SWE) %>%
                                 rasterFromXYZ())

  pred_train_rast <- purrr::map(1:5, ~ filter(folds, fold == .) %>%
                                  anti_join(preds, ., by = 'year') %>% # training years
                                  spread(year, SWE) %>%
                                  rasterFromXYZ())

  obs_train <- purrr::map(obs_train_rast, ~ anomalize(.) %>%
                            denoise(k = k_obs, weighted = TRUE, verbose = FALSE))

  pred_train <-  purrr::map(pred_train_rast, ~ anomalize(.) %>%
                             denoise(k = k_preds, weighted = TRUE, verbose = FALSE))

  # preprocess test data for each fold
  test <- purrr::map(1:5, ~ filter(folds, fold == .) %>%
                       semi_join(preds, ., by = 'year') %>%
                       spread(year, SWE) %>%
                       rasterFromXYZ())%>% # test years
    map2(pred_train_rast, ~.x - mean(.y)) # subtract the training predictor mean from the test predictors

  means <- purrr::map(obs_train_rast, mean)

  eots <- map2(pred_train, obs_train, ~eot(.x, .y, n = k_cca, standardised = FALSE, verbose = FALSE))

  # fit model to training data and use to predict new fields
  list(eot = eots, # this hard codes 10 patterns, but could be changed to min(k_preds, k_obs)
       test = test,
       mean = means)
}

predict_eot <- function(cv, k) {
  pmap(list(cv$eot, cv$test, cv$mean), ~predict(..1, ..2, n = k) + ..3) %>%
    brick
}


fit_cv <- function(dat, fun, k) {
  errors <- pmap_dfr(list(dat$train_preds, dat$train_obs, dat$test), fun, k = k) %>%
    arrange(x, y, year) %>%
    get_errors()

  c(rmse = summarise(errors, rmse = sqrt(mean(error^2))) %>% pull(rmse), # a na.rm =TRUE could go with mean, but best not to now
    corr = total_swe_corr(errors))
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


cv_eot_error <- function(rast) {
  dat <- rast %>%
    setNames(1982:2010) %>%
    as.data.frame(xy = TRUE, na.rm = TRUE, long = TRUE) %>%
    mutate(year = parse_number(layer)) %>%
    mutate_at(c('x', 'y'), round, digits = 5) %>%
    rename(SWE = value) %>%
    dplyr::select(x, y, year, SWE) %>%
    inner_join(  mutate_at(prism_dat, c('x', 'y'), round, digits = 5) , by = c('x', 'y', 'year'), suffix = c('_recon', '_obs'))

  rmse_dat<- dat %>%
    mutate(error = SWE_recon - SWE_obs) %>%
    summarise(rmse = sqrt(mean(error^2)))%>%
    pull(rmse)

  corr_dat<- dat %>%
    left_join( mutate_at(areas, c('x', 'y'), round, digits = 5), by = c("x", "y")) %>%
    group_by(year) %>%
    summarise(SWE_obs = sum(SWE_obs * area),
              SWE_recon = sum(SWE_recon * area)) %>%
    summarise(correlation = cor(SWE_recon, SWE_obs)) %>%
    pull(correlation)

  c(rmse = rmse_dat,
    corr = corr_dat)
}



# special predict method for cv that scales predictors (should fold into original predict method in future)
predict_gam <- function(preds, obs, newdata, k) {

  k_preds <- unique(preds$amplitudes$PC) %>% length()

  mod <- prep_data(preds, obs) %>%
    fit_model()

  new_pcs <- left_join(newdata, preds$climatology, by = c("x", "y")) %>%
    mutate(SWE = SWE - swe_mean) %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>%
    dplyr::select(-swe_mean, -swe_sd) %>%
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # scale according to training pcs
    .[,1:k_preds, drop = FALSE] %>%
    as_tibble(rownames = 'year')

  map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(year, PC, amplitude) %>%
    mutate(year = as.numeric(year)) %>%
    reconstruct_field(obs, .)  # use train_obs to get training pcs
}

# same as above but for pcr
predict_pcr <- function(preds, obs, newdata, k) {
# k is just here as a placeholder so the pmap command in fit_cv will run, but it doesn't have to be and should be removed
  k_preds <- unique(preds$amplitudes$PC) %>% length()

  mod <- prep_data(preds, obs) %>%
    fit_pcr()

  new_pcs <- left_join(newdata, preds$climatology, by = c("x", "y")) %>%
    mutate(SWE = SWE - swe_mean) %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>%
    dplyr::select(-swe_mean, -swe_sd) %>%
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # scale according to training pcs
    .[,1:k_preds, drop = FALSE] %>%
    as_tibble(rownames = 'year')

  map(mod$mod, ~add_predictions(new_pcs, ., var = 'amplitude', type = 'response')) %>%
    bind_rows(.id = 'PC') %>%
    dplyr::select(year, PC, amplitude) %>%
    mutate(year = as.numeric(year)) %>%
    reconstruct_field(obs, .)  # use train_obs to get training pcs
}



predict_cca <- function(preds, obs, newdata, k) {
  obs_amps <- obs$amplitudes %>%
                          mutate(PC = as.numeric(PC)) %>% # so spread doesn't lose order of columns
                          spread(PC, amplitude) %>%
                          select(-year) %>%
                          as.matrix()
  preds_amps <- preds$amplitudes %>%
                            mutate(PC = as.numeric(PC)) %>%
                            spread(PC, amplitude) %>%
                            select(-year) %>%
                            as.matrix()

  k_preds <- ncol(preds_amps)

  # train CCA on training data
  cca <- cancor(preds_amps, obs_amps, xcenter = FALSE, ycenter = FALSE)

  # predict CCA on new data
  new_amps <- left_join(newdata, preds$climatology, by = c("x", "y")) %>%
    mutate(SWE = SWE - swe_mean) %>%
    dplyr::select(-swe_mean, -swe_sd) %>%
    dplyr::mutate(SWE = SWE * sqrt(cos(y * pi / 180))) %>%
    tidyr::spread(year, SWE) %>%
    dplyr::select(-x, -y) %>%
    t() %>%
    predict(preds$pca, .) %>%
    scale(center = FALSE, scale = preds$pca$sdev) %>% # this should be from the training pcs
    .[,1:k_preds, drop = FALSE]

  # add a check here if k <= min(k_preds, k_obs)

  new_years <- unique(newdata$year)

  new_amps %*%
    cca$xcoef[,1:k, drop = FALSE] %*% # drop = FALSE prevents a vector returning when k = 1
    diag(cca$cor, nrow = k) %*%  # nrow likewise prevents weird behavior when k = 1
    solve(cca$ycoef)[1:k,,drop = FALSE] %>%
    `rownames<-`(new_years) %>%
    as_tibble(rownames = 'year') %>%
    gather(PC, amplitude, -year) %>%
    mutate(year = as.numeric(year)) %>%
    reconstruct_field(obs, .)
}

# check that the newdata argument is handled the same as not
# identical(cca_fit(prism_dat, cera_dat, k_preds = 4, k_obs = 4) ,
#          cca_fit(prism_dat, cera_dat, cera_dat, k_preds = 4, k_obs = 4) %>% filter(year >= 1982) )

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
