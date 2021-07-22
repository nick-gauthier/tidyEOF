



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

#make generic!
#' @export
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


# check that the newdata argument is handled the same as not
# identical(cca_fit(prism_dat, cera_dat, k_preds = 4, k_obs = 4) ,
#          cca_fit(prism_dat, cera_dat, cera_dat, k_preds = 4, k_obs = 4) %>% filter(year >= 1982) )



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






