library(survival)
## simulate survival data that follows a lognormal distribution...
## ...with a very specific structure:
##
## log T ~ alpha_0 + alpha_1 + alpha_1 * alpha_2 * x + epsilon
## 
## where x ~ N(0,1) and epsilon ~ N(0, sigma^2)

rlnorm_trial <- function(n,
                         alpha_0 = 2,
                         alpha_1 = 1,
                         alpha_2 = 0.5,
                         sigma = 1,
                         rec_period, 
                         max_t){
  
  x <- rnorm(n)
  
  time_star <- rnorm(n, alpha_0 + alpha_1 + alpha_1 * alpha_2 * x, sigma) |> exp()
  
  rec <- rec_period * runif(n)
  
  event <- time_star + rec < max_t
  
  time <- pmin(time_star, max_t - rec)
  
  data.frame(time = time,
             event = event,
             x = x,
             rec = rec,
             max_t = max_t)
}


## sim one trial 
sim_one_trial <- function(n){ 
  
  dat_0 <- rlnorm_trial(n = n, alpha_1 = 0, rec_period = 6, max_t = 26) |> mutate(treat = 0)
  dat_1 <- rlnorm_trial(n = n, alpha_1 = 1, rec_period = 6, max_t = 26) |> mutate(treat = 1)
  
  dat <- rbind(dat_0, dat_1)
  
  return(dat)
  
}


## sim one trial and estimate using 3 different methods
sim_study <- function(n = 50, tau = 24){ 
  
  
  dat <- sim_one_trial(n)
  
  ## apply RMST estimates
  res_1 <- RMST_stratCox(time=dat$time, 
                         status=dat$event, 
                         arm=dat$treat, 
                         covariates=dat$x, 
                         tau=tau) |> mutate(method = "stratCox")

  res_2 <- RMST_unadj(time=dat$time, 
                      status=dat$event,
                      arm=dat$treat, 
                      tau=tau) |> mutate(method = "unadj")

  res_3 <- RMST_Cox(time=dat$time,
                    status=dat$event, 
                    arm=dat$treat,
                    covariates=dat$x, 
                    tau=tau) |> mutate(method = "Cox")
  # 
  #   res_3 <- ARMST_iptw(time=dat$time, 
  #                       status=dat$event,
  #                       arm=dat$treat,
  #                       covariates=dat$strata,
  #                       tau=tau) |> mutate(method = "iptw")
  
  
  return(rbind(res_1, res_2, res_3))
  
}

##-----------------------------------------------------------


sim_one_exp_trial <- function(n_per_arm,
                              lambda_c,
                              hr,
                              rec_period,
                              max_t,
                              k = 1){
  
  time_star_c <- rexp(n_per_arm, rate = lambda_c)
  time_star_e <- rexp(n_per_arm, rate = lambda_c * hr)
  
  rec_c <- rec_period * runif(n_per_arm) ^ (1/k)
  rec_e <- rec_period * runif(n_per_arm) ^ (1/k)
  
  event_c <- as.numeric(time_star_c + rec_c < max_t)
  event_e <- as.numeric(time_star_e + rec_e < max_t)
  
  time_c <- pmin(time_star_c, max_t - rec_c)
  time_e <- pmin(time_star_e, max_t - rec_e)
  
  data.frame(time = c(time_c, time_e),
             event = c(event_c, event_e),
             treat = rep(c("c", "e"), each = n_per_arm))
}


##-----------------------------------------------------------


sim_one_weibull_trial <- function(n_per_arm,
                                  s_0,
                                  time_0,
                                  hr,
                                  shape,
                                  rec_period,
                                  max_t,
                                  k = 1){
  
  scale_0 <- time_0 / ((-log(s_0)) ^ (1 / shape))
  scale_1 <- scale_0 / (hr ^ (1 / shape))
  
  time_star_c <- rweibull(n_per_arm, shape = shape, scale = scale_0)
  time_star_e <- rweibull(n_per_arm, shape = shape, scale = scale_1)
  
  rec_c <- rec_period * runif(n_per_arm) ^ (1/k)
  rec_e <- rec_period * runif(n_per_arm) ^ (1/k)
  
  event_c <- as.numeric(time_star_c + rec_c < max_t)
  event_e <- as.numeric(time_star_e + rec_e < max_t)
  
  time_c <- pmin(time_star_c, max_t - rec_c)
  time_e <- pmin(time_star_e, max_t - rec_e)
  
  cen_c <- rexp(n_per_arm, rate = 0.0006)
  cen_e <- rexp(n_per_arm, rate = 0.0006)
  
  event_c <- as.numeric((event_c == 1) & (time_c < cen_c))
  event_e <- as.numeric((event_e == 1) & (time_e < cen_e))
  
  time_c <- pmin(time_c, cen_c)
  time_e <- pmin(time_e, cen_e)
  
  
  
  data.frame(time = c(time_c, time_e),
             event = c(event_c, event_e),
             treat = rep(c("c", "e"), each = n_per_arm))
}

##------------------------------------------------------------------------------

one_sim_zs <- function(n_per_arm,
                       s_0,
                       hr = 1 / 1.5,
                       rec_period = 0.0001,
                       max_t = 3.001,
                       tau_target = 3){
  
  dat <- sim_one_exp_trial(n_per_arm = n_per_arm,
                           lambda_c = -log(s_0) / 3,
                           hr = hr,
                           rec_period = rec_period,
                           max_t = max_t,
                           k = 1)
  
  dat$treat <- ifelse(dat$treat == "c", 0, 1)
  max_t_0 <- max(dat[dat$treat == 0, "time"])
  max_t_1 <- max(dat[dat$treat == 1, "time"])
  
  res_rmst <- rmst2(time = dat$time,
                    status = dat$event,
                    arm = dat$treat,
                    tau = min(tau_target, min(max_t_0, max_t_1)))
  
  tau = min(tau_target, min(max_t_0, max_t_1))
  n_events_tau = sum((dat$event == 1) & (dat$time > tau))
  n_events = sum(dat$event == 1)
  
  rmst_est <- res_rmst$unadjusted.result[1,1]
  rmst_se <- (res_rmst$unadjusted.result[1,2:3] |> diff()) / (2 * qnorm(0.975))
  z_rmst <- rmst_est / rmst_se
  
  res_cox <- coxph(Surv(time, event) ~ treat, data = dat)
  sum_cox <- summary(res_cox)
  wald_cox_z <- sum_cox$coef[,"z"]
  z_cox <- sqrt(sum_cox$sctest["test"]) * ifelse(wald_cox_z < 0, -1, 1)
  
  return(data.frame(z_rmst = z_rmst,
                    z_cox = z_cox,
                    tau = tau,
                    n_events = n_events,
                    n_events_tau = n_events_tau))
}















