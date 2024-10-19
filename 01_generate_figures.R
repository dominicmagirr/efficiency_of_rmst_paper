#--------------------------------------------------
# Main script to generate paper figures
#--------------------------------------------------

#-----------------------------------------------
# Install the specific version of flexsurv (version 2.2.2) for compatibility.
# This ensures the model fitting functions match the code and results expected.
#-----------------------------------------------
local_lib <- "./local_lib"
if (!dir.exists(local_lib)) {
  dir.create(local_lib)
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", lib = local_lib, repos = "http://cran.us.r-project.org")
}

if (!requireNamespace("flexsurv", quietly = TRUE, lib.loc = local_lib)) {
  devtools::install_version("flexsurv",
                            version = "2.2.2",
                            lib = local_lib,
                            repos = "http://cran.us.r-project.org")
}


#--------------------------------------------------
# Function to check if packages are installed, 
# and install if not
#--------------------------------------------------
check_and_install_packages <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg)
  }
  invisible(lapply(pkg, library, character.only = TRUE))
}

# List of required packages
required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "ggplot2",
  "survival",
  "clustermq",
  "nphRCT",
  "survRM2",
  "ggsurvfit",
  "gt", 
  "patchwork"
)


# Check and install missing packages
check_and_install_packages(required_packages)

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
