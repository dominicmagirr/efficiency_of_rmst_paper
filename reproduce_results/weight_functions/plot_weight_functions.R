#-------------------------------------------------------------
# Load required libraries, checking if installed first
#-------------------------------------------------------------
check_and_install_packages <- function(pkg) {
  new_pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new_pkg)) {
    install.packages(new_pkg)
  }
  invisible(lapply(pkg, library, character.only = TRUE))
}

required_packages <- c("dplyr", "purrr", "ggplot2")
check_and_install_packages(required_packages)


#-------------------------------------------------------------
# Main function for plotting weight functions for simulation scenarios
#-------------------------------------------------------------
plot_weight_f <- function(s_1 = 0.9, rec_period = 0.001, hr = 0.67, max_t = 3.001, tau_target = 3) {
  
  # Classify event rate based on s_1 value
  event_rate <- dplyr::case_when(
    s_1 == 0.9 ~ "Low",
    s_1 == 0.6 ~ "Moderate",
    s_1 == 0.2 ~ "High",
    TRUE ~ as.character(s_1)
  )
  
  # Classify recruitment rate based on rec_period value
  rec_rate <- dplyr::case_when(
    rec_period == 0.001 ~ "Instant",
    rec_period == 0.5 ~ "Fast",
    rec_period == 1.5 ~ "Moderate",
    rec_period == 2.5 ~ "Slow",
    TRUE ~ as.character(rec_period)
  )
  
  # Compute the control and experimental hazard rates
  lambda_e <- -log(s_1) / 3  # Experimental arm rate
  lambda_c <- lambda_e / hr  # Control arm rate (adjusted for hazard ratio)
  
  # Create a time sequence to evaluate the weight functions
  t_seq <- seq(0, tau_target, length.out = 100)
  
  ## Calculate survival functions
  # Censoring survival: proportion of subjects remaining in the trial at time t
  surv_c <- 1 - pmax(0, (t_seq - (max_t - rec_period)) / rec_period)
  
  # Control arm survival: exponential survival curve based on the control rate
  surv_0 <- exp(-lambda_c * t_seq)
  
  ## Calculate weights for the two tests
  # Log-Rank weight: product of censoring and control survival probabilities
  w_lr <- surv_c * surv_0
  
  # RMST weight: restricted mean survival time weight function
  w_rmst <- (exp(-lambda_c * t_seq) - exp(-lambda_c * tau_target)) / (1 - exp(-lambda_c * tau_target))
  
  # Normalize weights and prepare data for plotting
  out_dat <- data.frame(
    Time = rep(t_seq, 2),  # Time points
    Weight = c(w_lr / mean(w_lr), w_rmst / mean(w_rmst)),  # Normalized weights
    Test = rep(c("PH", "RMST"), each = length(t_seq)),  # Test types
    event_rate = event_rate,  # Event rate category
    rec_rate = rec_rate  # Recruitment rate category
  )
  
  return(out_dat)
}


#-------------------------------------------------------------
# Set up parameters for event rates and recruitment periods
#-------------------------------------------------------------

s_1_values <- c(0.9, 0.6, 0.2)  # Event rates (low, moderate, high)
rec_period_values <- c(0.001, 0.5, 1.5, 2.5)  # Recruitment periods (instant, fast, moderate, slow)

#-------------------------------------------------------------
# Generate combinations of event rates and recruitment periods using expand.grid
# Apply plot_weight_f across these combinations using purrr::map2_df
#-------------------------------------------------------------

p_res <- purrr::map2_df(
  expand.grid(s_1_values, rec_period_values)[, 1],  # Event rates
  expand.grid(s_1_values, rec_period_values)[, 2],  # Recruitment periods
  plot_weight_f
)


#-------------------------------------------------------------
# Relabel recruitment rate factors for better plot readability
#-------------------------------------------------------------
p_res$rec_rate <- factor(p_res$rec_rate,
                         levels = c("Instant", "Fast", "Moderate", "Slow"),
                         labels = c("Instant recruitment", "Fast recruitment", 
                                    "Moderate recruitment", "Slow recruitment"))

# Relabel event rate factors for better plot readability
p_res$event_rate <- factor(p_res$event_rate,
                           levels = c("Low", "Moderate", "High"),
                           labels = c("Low event rate", "Moderate event rate", 
                                      "High event rate"))

#-------------------------------------------------------------
# Generate the weight function plot with facets 
# for event rate and recruitment rate
#-------------------------------------------------------------
weight_function_plot <-
  ggplot(data = p_res, mapping = aes(x = Time, y = Weight, linetype = Test)) +
  geom_line() +
  facet_wrap(event_rate ~ rec_rate, scales = "free") +
#  theme_bw() +
  theme_minimal() +
  labs(x = "Time", y = "Weight", linetype = "Test") +
  theme(legend.position = "right")

#-------------------------------------------------------------
# Save the plot as a high-resolution PNG file
# Check if the "figs/" directory exists and create it if not
#-------------------------------------------------------------
if (!dir.exists("figs")) {
  dir.create("figs")
}

ggsave("figs/weight_functions.pdf", 
       weight_function_plot, 
       width = 10, 
       height = 9, 
       units = "in", 
       dpi = 300)