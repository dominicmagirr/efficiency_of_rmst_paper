library(survival)
##---------------------------------------------------------
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















