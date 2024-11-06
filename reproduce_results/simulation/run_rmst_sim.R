library(survival)
source("src/sim_helper_fns.R")
source("src/survRM2_fns.R")


##---------------------------------------------------------------------------------------
post_sim <- function(res){
  
  pow_rmst <- mean(res$z_rmst > qnorm(0.975))
  pow_ph <- mean(res$z_cox < qnorm(0.025))
  
  eff_rmst_cox <- ((qnorm(0.975) + qnorm(mean(res$z_rmst > qnorm(0.975)))) / (qnorm(0.975) + qnorm(mean(res$z_cox < qnorm(0.025))))) ^ 2
  
  frac_post_tau <- mean(res$n_events_tau / res$n_events)
  
  return(data.frame(pow_rmst = pow_rmst |> round(2),
                    pow_ph = pow_ph|> round(2),
                    eff_rmst_cox = eff_rmst_cox|> round(2),
                    frac_post_tau = frac_post_tau |> round(2)))
}



## Simulate 10,000 runs
library(clustermq)
run_sim_wrapper <- function(n_per_arm, s_0, rec_period, max_t){
  
  if(!exists('worker_setup_complete', mode='logical')) {
    cat("Calling setup function:\n")
    
    source("src/sim_helper_fns.R")
    source("src/survRM2_fns.R")
    
    
    worker_setup_complete <<- TRUE
  } else {
    cat("Calling GC:\n")
    print(gc())
  }
  
  one_sim_zs(n_per_arm = n_per_arm, s_0 = s_0, rec_period = rec_period, max_t = max_t)
  
}




##------------------------- Reproduce Freidlin et al. ------------------------------------------------
N_runs <- 1e4

set.seed(32624)

## option to use clustermq or not...

#res_low_instant <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 0.0001, max_t = 3.001)
res_low_instant <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 0.0001, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_mod_instant <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 0.0001, max_t = 3.001)
res_mod_instant <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 0.0001, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_high_instant <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 0.0001, max_t = 3.001)
res_high_instant <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 0.0001, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_instant)
post_sim(res_mod_instant)
post_sim(res_high_instant)


##-----------------------------   Fast recruitment --------------------------------------------

#res_low_fast <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 0.5, max_t = 3.001)
res_low_fast <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 0.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_mod_fast <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 0.5, max_t = 3.001)
res_mod_fast <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 0.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_high_fast <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 0.5, max_t = 3.001)
res_high_fast <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 0.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_fast)
post_sim(res_mod_fast)
post_sim(res_high_fast)
##------------------------------------ Moderate recruitment ---------------------------------------------------

#res_low_mod <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 1.5, max_t = 3.001)
res_low_mod <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 1.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_mod_mod <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 1.5, max_t = 3.001)
res_mod_mod <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 1.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_high_mod <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 1.5, max_t = 3.001)
res_high_mod <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 1.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_mod)
post_sim(res_mod_mod)
post_sim(res_high_mod)


##------------------------------------ Slow recruitment ---------------------------------------------------

#res_low_slow <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 2.5, max_t = 3.001)
res_low_slow <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 2.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_mod_slow <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 2.5, max_t = 3.001)
res_mod_slow <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 2.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)
#res_high_slow <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 2.5, max_t = 3.001)
res_high_slow <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 2.5, max_t = 3.001), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_slow)
post_sim(res_mod_slow)
post_sim(res_high_slow)


##------------------------------------- Extra Data ----------------------------------------------
set.seed(32624)

#res_low_instant_2 <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 0.0001, max_t = 3.5)
res_low_instant_2 <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 0.0001, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_mod_instant_2 <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 0.0001, max_t = 3.5)
res_mod_instant_2 <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 0.0001, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_high_instant_2 <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 0.0001, max_t = 3.5)
res_high_instant_2 <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 0.0001, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------


post_sim(res_low_instant_2)
post_sim(res_mod_instant_2)
post_sim(res_high_instant_2)

##---------------------------------------------------------------------------------------

#res_low_fast_2 <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 0.5, max_t = 3.5)
res_low_fast_2 <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 0.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_mod_fast_2 <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 0.5, max_t = 3.5)
res_mod_fast_2 <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 0.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_high_fast_2 <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 0.5, max_t = 3.5)
res_high_fast_2 <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 0.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_fast_2)
post_sim(res_mod_fast_2)
post_sim(res_high_fast_2)

##----------------------------------------------------------------------------------------

#res_low_mod_2 <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 1.5, max_t = 3.5)
res_low_mod_2 <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 1.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_mod_mod_2 <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 1.5, max_t = 3.5)
res_mod_mod_2 <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 1.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_high_mod_2 <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 1.5, max_t = 3.5)
res_high_mod_2 <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 1.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_mod_2)
post_sim(res_mod_mod_2)
post_sim(res_high_mod_2)

##---------------------------------------------------------------------------------------

#res_low_slow_2 <- purrr::map_df(rep(1000, 1e4), one_sim_zs, s_0 = exp(log(0.9) * 1.5), rec_period = 2.5, max_t = 3.5)
res_low_slow_2 <- Q(run_sim_wrapper, n_per_arm=rep(1000, N_runs),const = list(s_0 = exp(log(0.9)*1.5), rec_period = 2.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_mod_slow_2 <- purrr::map_df(rep(250, 1e4), one_sim_zs, s_0 = exp(log(0.6) * 1.5), rec_period = 2.5, max_t = 3.5)
res_mod_slow_2 <- Q(run_sim_wrapper, n_per_arm=rep(250, N_runs),const = list(s_0 = exp(log(0.6)*1.5), rec_period = 2.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)
#res_high_slow_2 <- purrr::map_df(rep(150, 1e4), one_sim_zs, s_0 = exp(log(0.2) * 1.5), rec_period = 2.5, max_t = 3.5)
res_high_slow_2 <- Q(run_sim_wrapper, n_per_arm=rep(150, N_runs),const = list(s_0 = exp(log(0.2)*1.5), rec_period = 2.5, max_t = 3.5), n_jobs = 200) |> do.call(what = rbind)

##---------------------------------------------------------------------------------------

post_sim(res_low_slow_2)
post_sim(res_mod_slow_2)
post_sim(res_high_slow_2)


