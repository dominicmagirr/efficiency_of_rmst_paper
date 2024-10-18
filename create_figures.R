## produce paper figures


# Function to check if packages are installed, and install if not
check_and_install_packages <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg)
  }
  invisible(lapply(pkg, library, character.only = TRUE))
}

# List of required packages
required_packages <- c("dplyr", "purrr", "ggplot2", "survival", "clustermq")

                      
# Check and install missing packages
check_and_install_packages(required_packages)



#-------------------------------------------------------------
# Create RMST and PH weight function figure 
#-------------------------------------------------------------
source("reproduce_results/weight_functions/weight_functions.R")

