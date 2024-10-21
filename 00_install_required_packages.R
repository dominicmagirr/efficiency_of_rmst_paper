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
