library(survival)
library(clustermq)

# Load helper functions
source("src/sim_helper_fns.R")
source("src/survRM2_fns.R")

#-------------------------------------------------------------
# Global Variables
#-------------------------------------------------------------
N_runs <- 1e4  # Number of simulations
n_jobs <- 200  # Number of parallel jobs
seed_value <- 32624  # Set seed for reproducibility

#-------------------------------------------------------------
# Simulation parameters for event rates and recruitment
#-------------------------------------------------------------
# Event rates
s_0_values <- c(exp(log(0.9) * 1.5), exp(log(0.6) * 1.5), exp(log(0.2) * 1.5))  

# Recruitment periods
rec_period_values <- c(0.0001, 0.5, 1.5, 2.5)  

# Sample sizes for different event rates
n_per_arm_values <- c(1000, 250, 150)  

# Maximum follow-up time
max_t <- 3.001  

#-------------------------------------------------------------
# Post Simulation Summary Function
# Calculate evluation metrics 
#-------------------------------------------------------------
post_sim <- function(res) {
  pow_rmst <- mean(res$z_rmst > qnorm(0.975))
  
  pow_ph <- mean(res$z_cox < qnorm(0.025))
  
  eff_rmst_cox <- ((qnorm(0.975) + qnorm(mean(res$z_rmst > qnorm(0.975)))) / 
                     (qnorm(0.975) + qnorm(mean(res$z_cox < qnorm(0.025))))) ^ 2
  
  frac_post_tau <- mean(res$n_events_tau / res$n_events)
  
  return(data.frame(
    pow_rmst = round(pow_rmst, 2),
    pow_ph = round(pow_ph, 2),
    eff_rmst_cox = round(eff_rmst_cox, 2),
    frac_post_tau = round(frac_post_tau, 2)
  ))
}

#-------------------------------------------------------------
# Simulation Wrapper
#-------------------------------------------------------------
run_sim_wrapper <- function(n_per_arm, s_0, rec_period, max_t) {
  
  if (!exists('worker_setup_complete', mode = 'logical')) {
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

#-------------------------------------------------------------
# Main Simulation Function
#-------------------------------------------------------------
run_simulation <- function(s_0, rec_period, n_per_arm, max_t) {
  Q(run_sim_wrapper,
    n_per_arm = rep(n_per_arm, N_runs),
    const = list(s_0 = s_0, rec_period = rec_period, max_t = max_t),
    n_jobs = n_jobs) |> do.call(what = rbind)
}

#-------------------------------------------------------------
# Run Simulations for All Scenarios
#-------------------------------------------------------------
set.seed(seed_value)

# Function to run all scenarios for a given recruitment period
run_scenarios <- function(rec_period, max_t) {
  results <- list()
  
  # Run simulations for different event rates and sample sizes
  for (i in seq_along(s_0_values)) {
    cat("Running scenario for s_0 =", s_0_values[i], "and rec_period =", rec_period, "\n")
    res <- run_simulation(s_0_values[i], rec_period, n_per_arm_values[i], max_t)
    results[[i]] <- post_sim(res)
  }
  
  return(do.call(rbind, results))
}

# Run simulations for each recruitment period
all_results <- list()

for (rec_period in rec_period_values) {
  all_results[[as.character(rec_period)]] <- run_scenarios(rec_period, max_t)
}

#-------------------------------------------------------------
# Display Results
#-------------------------------------------------------------
print("Results for different recruitment periods:")
print(all_results)