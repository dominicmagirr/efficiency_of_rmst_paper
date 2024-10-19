#-----------------------------------------------
# Load required libraries
#-----------------------------------------------
library(survival)
library(flexsurv, lib.loc = local_lib)
library(ggsurvfit)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

#-----------------------------------------------
# Source custom functions
#-----------------------------------------------
source("src/survRM2_fns.R")

#-----------------------------------------------
# Define a function to fit the spline model and 
# generate hazard ratio confidence intervals
#-----------------------------------------------
fit_and_plot <- function(data,
                         time_points,
                         plot_title,
                         ylim_range,
                         xlim_range,
                         abline_value) {
  # Fit spline model
  fit_spline <- flexsurvspline(
    Surv(time, event) ~ arm +
      gamma1(arm) +
      gamma2(arm),
    data = data,
    k = 4,
    scale = "hazard"
  )
  
  # Simulate from multivariate normal distribution
  sims <- normboot.flexsurvreg(fit_spline, B = 1e5, raw = TRUE)
  
  # Generate hazard ratio confidence intervals
  res_mat <- matrix(0, nrow = length(time_points), ncol = 3)
  for (i in time_points) {
    h_c <- hsurvspline(i, gamma = sims[, c("gamma0", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5")], knots = fit_spline$knots)
    h_e <- hsurvspline(
      i,
      gamma = cbind(sims[, "gamma0"], sims[, "gamma1"] + sims[, "gamma1(armexperimental)"], sims[, "gamma2"] + sims[, "gamma2(armexperimental)"], sims[, "gamma3"], sims[, "gamma4"], sims[, "gamma5"]),
      offset = sims[, "armexperimental"],
      knots = fit_spline$knots
    )
    hr <- h_e / h_c
    res_mat[i, ] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    Time = time_points,
    Median = res_mat[, 2],
    Lower = res_mat[, 1],
    Upper = res_mat[, 3]
  )
  
  # Plot using ggplot2 with confidence interval as a band
  ggplot(plot_data, aes(x = Time)) +
    geom_line(aes(y = Median), color = "black", alpha = 0.9) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper),
                alpha = 0.2,
                fill = "red") +
    geom_hline(
      yintercept = 1,
      linetype = "solid",
      color = "grey",
      alpha = 0.6
    ) +
    labs(title = plot_title, x = "Time", y = "Hazard Ratio") +
    theme_minimal() +
    coord_cartesian(ylim = ylim_range, xlim = xlim_range) +
    scale_y_log10()
}

#-----------------------------------------------
# Load datasets and apply the function to generate plots
#-----------------------------------------------

# Dataset 1: sustain_ipd.csv
dat1 <- read.csv("data/sustain_ipd.csv")
dat1$arm <- factor(ifelse(dat1$arm == 1, "control", "experimental"))
plot1 <- fit_and_plot(dat1, 1:109, "(A)", c(0.3, 1.5), c(5, 109), 0.7373)

# Dataset 2: CLEOPATRA_2A.rda
load("data/CLEOPATRA_2A.rda")
dat2 <- CLEOPATRA_2A
dat2$arm <- factor(ifelse(dat2$arm == "control", "control", "experimental"))
plot2 <- fit_and_plot(dat2, 1:70, "(B)", c(0.3, 1.5), c(5, 70), 0.6767)

# Dataset 3: leader_km.csv
dat3 <- read.csv("data/leader_km.csv")
dat3$arm <- factor(ifelse(dat3$arm == 1, "control", "experimental"))
plot3 <- fit_and_plot(dat3, 1:50, "(C)", c(0.5, 1.2), c(1, 50), 0.869)

# Dataset 4: External Excel data
dat4 <- read_csv("data/poplar.csv")
plot4 <- fit_and_plot(dat4, 1:30, "(D)", c(0.1, 1.2), c(1, 30), 0.6752)

#-----------------------------------------------
# Combine all four plots using patchwork
#-----------------------------------------------

combined_plot <- plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 2)
print(combined_plot)

#-----------------------------------------------
# Save the combined plot
#-----------------------------------------------
ggsave("figs/hr_case_studies.pdf",
       combined_plot,
       width = 12,
       height = 10)
