# Load necessary libraries
library(survival)
library(flexsurv)
library(nphRCT)
library(dplyr)
library(survRM2)
library(ggsurvfit)
library(readr)
library(gt)

#-------------------------------------------
# Step 1: Load Data for All Studies

# Example A - Sustain Case Study
data_sustain <- read.csv("data/sustain_ipd.csv")

# Example B - Cleopatra Study
load("data/CLEOPATRA_2A.rda")
data_cleopatra <- CLEOPATRA_2A

# Example C - Leader Study
data_leader <- read.csv("data/leader_km.csv")

# Example D - Poplar Study
data_poplar <- read_csv("data/poplar.csv")

#-------------------------------------------
# Step 2: Data Preparation for All Studies

# Sustain Case Study
data_sustain$arm <- factor(ifelse(data_sustain$arm == 1, "control", "experimental"))

# Cleopatra Study
data_cleopatra$arm <- factor(ifelse(data_cleopatra$arm == "control", "control", "experimental"))

# Leader Study
data_leader$arm <- factor(ifelse(data_leader$arm == 1, "control", "experimental"))

# Poplar Study
data_poplar$arm <- factor(ifelse(data_poplar$arm == "control", "control", "experimental"))

#-------------------------------------------
# Step 3: Define Tau Values for All Studies

tau_sustain <- 108
tau_cleopatra <- 65
tau_leader <- 48
tau_poplar <- 24

#-------------------------------------------
# Step 4: Fit Cox Proportional Hazards Model for All Studies

fit_cox_model <- function(data) {
  cox_fit <- coxph(Surv(time, event) ~ arm, data = data)
  summary_cox <- summary(cox_fit)
  return(summary_cox$coefficients[1, "z"])
}

# Calculate Cox Z-statistic for each study
cox_z_sustain <- fit_cox_model(data_sustain)
cox_z_cleopatra <- fit_cox_model(data_cleopatra)
cox_z_leader <- fit_cox_model(data_leader)
cox_z_poplar <- fit_cox_model(data_poplar)

#-------------------------------------------
# Step 5: Calculate RMST Z-Statistic for All Studies

fit_rmst_model <- function(data, tau) {
  data$arm2 <- ifelse(data$arm == "control", 0, 1)
  rmst_res <- rmst2(time = data$time, status = data$event, arm = data$arm2, tau = tau)
  rmst_est <- rmst_res$unadjusted.result[1, 1]
  rmst_se <- (rmst_res$unadjusted.result[1, 2:3] |> diff()) / (2 * qnorm(0.975))
  z_rmst <- rmst_est / rmst_se
  return(z_rmst)
}

# Calculate RMST Z-statistic for each study
rmst_z_sustain <- fit_rmst_model(data_sustain, tau_sustain)
rmst_z_cleopatra <- fit_rmst_model(data_cleopatra, tau_cleopatra)
rmst_z_leader <- fit_rmst_model(data_leader, tau_leader)
rmst_z_poplar <- fit_rmst_model(data_poplar, tau_poplar)

#-------------------------------------------
# Step 6: Calculate Percentage of Events After Tau for All Studies

calculate_percent_events_after_tau <- function(data, tau) {
  percent_events <- sum(data$time > tau & data$event == 1) / sum(data$event) * 100
  return(round(percent_events, 0))
}

# Calculate percentage of events after tau for each study
percent_events_sustain <- calculate_percent_events_after_tau(data_sustain, tau_sustain)
percent_events_cleopatra <- calculate_percent_events_after_tau(data_cleopatra, tau_cleopatra)
percent_events_leader <- calculate_percent_events_after_tau(data_leader, tau_leader)
percent_events_poplar <- calculate_percent_events_after_tau(data_poplar, tau_poplar)

#-------------------------------------------
# Step 7: Create Results Table for All Studies

results_table <- data.frame(
  Example = c("A - Sustain", "B - Cleopatra", "C - Leader", "D - Poplar"),
  Cox_Z_Statistic = round(c(cox_z_sustain, cox_z_cleopatra, cox_z_leader, cox_z_poplar), 2),
  RMST_Z_Statistic = round(c(rmst_z_sustain, rmst_z_cleopatra, rmst_z_leader, rmst_z_poplar), 2),
  Tau = c(tau_sustain, tau_cleopatra, tau_leader, tau_poplar),
  Percent_Events_After_Tau = c(percent_events_sustain, percent_events_cleopatra, percent_events_leader, percent_events_poplar)
)

# Print the final results table
print(results_table)

#-------------------------------------------
# Step 8: Create GT Table and Save to RTF

gt_results_table <- results_table |> 
  gt() |> 
  tab_header(
    title = "Summary of Cox PH and RMST Analyses for All Studies"
  )


gt::gtsave(gt_results_table, "figs/km_rmst_case_studies.rtf")
gt::gtsave(gt_results_table, "figs/km_rmst_case_studies.tex")

