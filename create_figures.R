## Produce paper figures

# Function to check if packages are installed, and install if not
check_and_install_packages <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg)
  }
  invisible(lapply(pkg, library, character.only = TRUE))
}

# List of required packages
required_packages <- c("dplyr",
                       "tidyr",
                       "readr",
                       "purrr",
                       "ggplot2",
                       "survival",
                       "clustermq")

# Check and install missing packages
check_and_install_packages(required_packages)

#-------------------------------------------------------------
# Create RMST and PH weight function figure 
#-------------------------------------------------------------
source("reproduce_results/weight_functions/weight_functions.R")

#-------------------------------------------------------------
# Run simulation scenarios and generate results table (optional)
#-------------------------------------------------------------
run_simulation <- FALSE  # Set to TRUE to run the simulation code

if (run_simulation) {
  # Run simulations (Warning: this will take time if not on a cluster)
  source("reproduce_results/simulation/run_rmst_sim_cluster.R")
} else {
  message("Skipping simulation as 'run_simulation' is set to FALSE.")
}

# Generate results table
source("reproduce_results/simulation/generate_results_table.R")