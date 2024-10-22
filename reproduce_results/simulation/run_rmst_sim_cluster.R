#--------------------------------------------------------
# Load required packages and helper functions
library(survival)          # Survival analysis functions
library(clustermq)         # Parallel computing for simulations
library(dplyr)             # Data manipulation


# Load helper functions for simulations
source("src/sim_helper_fns.R")
source("src/survRM2_fns.R")

#--------------------------------------------------------
# Global Configuration for simulations
N_runs <- 1e4            # Number of simulation runs
n_jobs <- 200            # Number of parallel jobs for clustermq
random_seed <- 32624     # Seed for reproducibility
set.seed(random_seed)    # Set seed

#--------------------------------------------------------
# Simulation scenario metadata

global_params <- list(
  event_rates = c("Low" = 0.9, "Moderate" = 0.6, "High" = 0.2),
  recruitment_speeds = c("Instant" = 0.0001, "Fast" = 0.5, "Moderate" = 1.5, "Slow" = 2.5),
  max_t = c(3.001, 4.5),
  sample_sizes = c("Low" = 1000, "Moderate" = 250, "High" = 150)
)

#--------------------------------------------------------
# Additional simulation functions 

# Function to process simulation results and calculate metrics 
post_sim <- function(res){
  # Calculate power for RMST
  pow_rmst <- mean(res$z_rmst > qnorm(0.975))
  # Calculate power for PH
  pow_ph <- mean(res$z_cox < qnorm(0.025))
  # Calculate relative efficiency between RMST and PH
  eff_rmst_cox <- ((qnorm(0.975) + qnorm(mean(res$z_rmst > qnorm(0.975)))) / 
                     (qnorm(0.975) + qnorm(mean(res$z_cox < qnorm(0.025))))) ^ 2
  # Fraction of events occurring post tau
  frac_post_tau <- mean(res$n_events_tau / res$n_events)
  
  # Return calculated metrics as a data frame
  return(data.frame(pow_rmst = pow_rmst,
                    pow_ph = pow_ph,
                    eff_rmst_cox = eff_rmst_cox,
                    frac_post_tau = frac_post_tau))
}

# Wrapper function to run a single simulation scenario 
run_sim_wrapper <- function(n_per_arm, s_0, rec_period, max_t){
  # Load required setup in worker nodes if not already loaded
  if(!exists('worker_setup_complete', mode='logical')) {
    cat("Calling setup function:\n")
    source("src/sim_helper_fns.R")
    source("src/survRM2_fns.R")
    worker_setup_complete <<- TRUE
  } else {
    # Run garbage collection to free memory
    cat("Calling GC:\n")
    print(gc())
  }
  
  # Run one simulation scenario
  one_sim_zs(n_per_arm = n_per_arm, s_0 = s_0, rec_period = rec_period, max_t = max_t)
}

# Function to collect results across all simulation scenarios
collect_results <- function(){
  results <- list()
  scenario_id <- 1
  
  # Loop through each combination of event rate, recruitment speed, and max_t
  for (event_rate in names(global_params$event_rates)) {
    for (recruitment_speed in names(global_params$recruitment_speeds)) {
      for (t in global_params$max_t) {
        # Print information about the current scenario
        cat(sprintf("Running Scenario %d: Event Rate = %s, Recruitment Speed = %s, Max t = %.3f\n", 
                    scenario_id, event_rate, recruitment_speed, t))
        
        # Retrieve parameters for the current scenario
        n_per_arm <- global_params$sample_sizes[[event_rate]]
        s_0 <- exp(log(global_params$event_rates[[event_rate]]) * 1.5)
        rec_period <- global_params$recruitment_speeds[[recruitment_speed]]
        
        # Run the simulations in parallel using clustermq
        res <- Q(run_sim_wrapper, n_per_arm = rep(n_per_arm, N_runs),
                 const = list(s_0 = s_0, rec_period = rec_period, max_t = t), 
                 n_jobs = n_jobs) |> do.call(what = rbind)
        
        # Post-process simulation results
        sim_results <- post_sim(res)
        # Add scenario metadata
        sim_results$scenario_id <- scenario_id
        sim_results$event_rate <- event_rate
        sim_results$recruitment_speed <- recruitment_speed
        sim_results$max_t <- t
        
        # Store the results
        results[[scenario_id]] <- sim_results
        scenario_id <- scenario_id + 1
      }
    }
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, results)
  return(final_results)
}

#--------------------------------------------------------
# Run Simulations and Collect Results
final_results <- collect_results()

#--------------------------------------------------------
# Print or Save Results
print(final_results) # For printing in console

#--------------------------------------------------------
# Create Results Directory and Save with Metadata
results_folder <- "results"
if (!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Save the final results with metadata 
file_name <- sprintf("%s/simulation_results.csv", results_folder)
write_csv(final_results, file_name)
