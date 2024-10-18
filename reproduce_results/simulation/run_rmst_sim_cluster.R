library(survival)
source("src/sim_helper_fns.R")
source("src/survRM2_fns.R")

## Global Configuration
N_runs <- 1e4
n_jobs <- 200
random_seed <- 32624
set.seed(random_seed)

## store to capture run date and time
current_date <- format(Sys.Date(), "%Y%m%d")

global_params <- list(
  event_rates = c("Low" = 0.9, "Moderate" = 0.6, "High" = 0.2),
  recruitment_speeds = c("Instant" = 0.0001, "Fast" = 0.5, "Moderate" = 1.5, "Slow" = 2.5),
  max_t = c(3.001, 4.5),
  sample_sizes = c("Low" = 1000, "Moderate" = 250, "High" = 150)
)


## Function Definitions
post_sim <- function(res){
  pow_rmst <- mean(res$z_rmst > qnorm(0.975))
  pow_ph <- mean(res$z_cox < qnorm(0.025))
  eff_rmst_cox <- ((qnorm(0.975) + qnorm(mean(res$z_rmst > qnorm(0.975)))) / 
                     (qnorm(0.975) + qnorm(mean(res$z_cox < qnorm(0.025))))) ^ 2
  frac_post_tau <- mean(res$n_events_tau / res$n_events)
  
  return(data.frame(pow_rmst = pow_rmst,
                    pow_ph = pow_ph,
                    eff_rmst_cox = eff_rmst_cox,
                    frac_post_tau = frac_post_tau))
}

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

collect_results <- function(){
  results <- list()
  scenario_id <- 1
  
  for (event_rate in names(global_params$event_rates)) {
    for (recruitment_speed in names(global_params$recruitment_speeds)) {
      for (t in global_params$max_t) {
        cat(sprintf("Running Scenario %d: Event Rate = %s, Recruitment Speed = %s, Max t = %.3f\n", 
                    scenario_id, event_rate, recruitment_speed, t))
        
        n_per_arm <- global_params$sample_sizes[[event_rate]]
        s_0 <- exp(log(global_params$event_rates[[event_rate]]) * 1.5)
        rec_period <- global_params$recruitment_speeds[[recruitment_speed]]
        
        res <- Q(run_sim_wrapper, n_per_arm = rep(n_per_arm, N_runs),
                 const = list(s_0 = s_0, rec_period = rec_period, max_t = t), 
                 n_jobs = n_jobs) |> do.call(what = rbind)
        
        sim_results <- post_sim(res)
        sim_results$scenario_id <- scenario_id
        sim_results$event_rate <- event_rate
        sim_results$recruitment_speed <- recruitment_speed
        sim_results$max_t <- t
        
        results[[scenario_id]] <- sim_results
        scenario_id <- scenario_id + 1
      }
    }
  }
  
  final_results <- do.call(rbind, results)
  return(final_results)
}

## Run Simulations and Collect Results
final_results <- collect_results()

## Print or Save Results
print(final_results) # For printing in console

## Create Results Directory and Save with Metadata
results_folder <- "results"
if (!dir.exists(results_folder)) {
  dir.create(results_folder)
}

## save the final results 
file_name <- sprintf("%s/simulation_results_seed_%d_%s.csv", results_folder, random_seed, current_date)
write.csv(final_results, file_name, row.names = FALSE)


