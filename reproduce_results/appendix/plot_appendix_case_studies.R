# Function to generate survival curve, cumulative hazards and hazard ratio function plots for case studies

fit_and_plot <- function(data, file_name, t, xlim_surv, ylim_surv, ylim_hr, xlim_hr) {
  dat$arm <- factor(ifelse(dat$arm %in% c(1, "control"), "control", "experimental"))
  
  # Fit a flexible parametric survival model with spline terms for 'arm'
  fit_spline <- flexsurvspline(Surv(time, event) ~ arm + gamma1(arm) + gamma2(arm), data = dat, k = 4, scale = "hazard")
  
  # Calculate the overall hazard ratio using a Cox proportional hazards model
  hr_overall <- exp(coxph(Surv(time, event) ~ arm, data = dat)$coef)
  
  # Create a PDF file to save the plots
  pdf(file = file_name, width = 12, height = 10)
  par(mfrow = c(2, 2))
  
  # Plot the survival curves and fitted spline without confidence intervals
  
  plot(fit_spline, type = "survival", ci = FALSE, col = c(2, 4), ylim = ylim_surv, xlab = "Time", ylab = "Survival", main = "(1)")
  legend("bottomleft", c("Control", "Experimental"), col = c(2, 4), lty = c(1, 1))
  
  # Plot the cumulative hazard curves with confidence intervals
  plot(fit_spline, type = "cumhaz", ci = TRUE, col = c(2, 4), xlab = "Time", ylab = "Cumulative Hazard", main = "(2)")
  
  # Plot the cumulative hazard curves on a log-log scale
  plot(fit_spline, type = "cumhaz", log = "xy", col = c(2, 4), xlim = xlim_surv, xlab = "Time", ylab = "Cumulative Hazard", main = "(3)")
  
  sims <- normboot.flexsurvreg(fit_spline, B = 1e5, raw = TRUE)
  res_mat <- matrix(0, nrow = t, ncol = 3)
  
  for (i in 1:t) {
    h_c <- hsurvspline(i, gamma = sims[, c("gamma0", "gamma1", "gamma2", "gamma3", "gamma4", "gamma5")], knots = fit_spline$knots)
    h_e <- hsurvspline(i, gamma = cbind(sims[, "gamma0"], sims[, "gamma1"] + sims[, "gamma1(armexperimental)"], sims[, "gamma2"] + sims[, "gamma2(armexperimental)"], sims[, "gamma3"], sims[, "gamma4"], sims[, "gamma5"]), offset = sims[, "armexperimental"], knots = fit_spline$knots)
    hr <- h_e / h_c
    res_mat[i, ] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  }
  
  # Plot the hazard ratio over time with confidence intervals
  plot(1:t, res_mat[, 2], type = "l", ylim = ylim_hr, xlim = xlim_hr, xlab = "Time", ylab = "Hazard Ratio", main = "(4)")
  points(1:t, res_mat[, 1], type = "l", lty = 2)
  points(1:t, res_mat[, 3], type = "l", lty = 2)
  
  # Add a horizontal line representing the overall hazard ratio
  abline(h = hr_overall, col = 1, lty = 3)
  
  # Close the PDF file
  dev.off()
}

# Sustain study
dat <- read.csv("data/sustain_ipd.csv")
fit_and_plot(dat, "figs/appendix_sustain.pdf", t = 109, c(10, 109), c(0.9, 1), c(0.3, 1.5), c(5, 109))

# Cleopatra study
dat <- read.csv("data/CLEOPATRA_2A.csv")
fit_and_plot(dat, "figs/appendix_cleopatra.pdf", t = 70, c(5, 70), c(0, 1), c(0.3, 1.5), c(5, 70))

# Leader study
dat <- read.csv("data/leader_km.csv")
fit_and_plot(dat, "figs/appendix_leader.pdf", t = 50, c(5, 50), c(0.8, 1), c(0.5, 1.2), c(1, 50))

# Poplar study
dat <- read_csv("data/poplar.csv")
fit_and_plot(dat, "figs/appendix_poplar.pdf", t = 30, c(1, 30), c(0, 1), c(0.1, 1.2), c(1, 30))
