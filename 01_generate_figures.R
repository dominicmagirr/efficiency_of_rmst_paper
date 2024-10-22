#--------------------------------------------------
# Main script to generate paper figures
#--------------------------------------------------

#-------------------------------------------------------------
# Create RMST and PH weight function figure 
#-------------------------------------------------------------
source("reproduce_results/weight_functions/plot_weight_functions.R")

#-------------------------------------------------------------
# Run simulation scenarios and generate results table (optional)
#-------------------------------------------------------------
run_simulation <- FALSE  # Set to TRUE to run the simulation code

# Run simulations (Warning: this will take time if not on a cluster)
if (run_simulation) {
  source("reproduce_results/simulation/run_rmst_sim_cluster.R")
} else {
  message("Skipping simulation as 'run_simulation' is set to FALSE.")
}

# Generate results table
source("reproduce_results/simulation/generate_results_table.R")

#-------------------------------------------------------------
# Generate case study KM plots 
#-------------------------------------------------------------
source("reproduce_results/case_studies/plot_km_case_studies.R")

#-------------------------------------------------------------
# Generate case study KM vs RMST table
#-------------------------------------------------------------
source("reproduce_results/case_studies/calc_km_rmst_case_studies.R")

#-------------------------------------------------------------
# Generate case study HR plots
#-------------------------------------------------------------
source("reproduce_results/case_studies/plot_hr_case_studies.R")

#-------------------------------------------------------------
# Generate case study standardised score plots
#-------------------------------------------------------------
source("reproduce_results/case_studies/plot_std_scores_case_studies.R")

#-------------------------------------------------------------
# Generate appendix plots - case study 
#-------------------------------------------------------------
source("reproduce_results/appendix/plot_appendix_case_studies.R")



