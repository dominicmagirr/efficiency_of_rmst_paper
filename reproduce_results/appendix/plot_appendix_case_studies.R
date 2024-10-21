library(survival)
library(flexsurv, lib.loc = "local_lib") # Load the specific version of flexsurv
library(readr)

# A - sustain study 

dat <- read.csv("data/sustain_ipd.csv")
dat$arm <- factor(ifelse(dat$arm == 1, "control", "experimental"))


## assess ph
fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")

# Step 1: Call the pdf command to start the plot
pdf(file = "figs/appendix_sustain.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plots with R code

par(mfrow = c(2,2))
plot(fit_spline, type = "survival", ci = FALSE, col = c(2,4), ylim = c(0.9, 1), xlab = "Time", ylab = "Survival", main = "(A)")
legend("bottomleft", c("Control","Experimental"), col = c(2,4), lty = c(1,1))
plot(fit_spline, type = "cumhaz", ci = TRUE, col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(B)")
plot(fit_spline, type = "cumhaz", log = "xy", xlim = c(10, 109), col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(C)")


## conf int for hazard ratio
## simulate from MVN distribution
sims <- normboot.flexsurvreg(fit_spline, 
                             B = 1e5, 
                             raw = TRUE)

res_mat <- matrix(0, nrow = 109, ncol = 3)
for (i in 1:109){
  ## hazard
  h_c <-  hsurvspline(i, gamma = sims[,c("gamma0", 
                                         "gamma1", 
                                         "gamma2", 
                                         "gamma3", 
                                         "gamma4", 
                                         "gamma5")], 
                      knots = fit_spline$knots)
  
  ## hazard
  h_e <-  hsurvspline(i, 
                      gamma = cbind(sims[,"gamma0"], 
                                    sims[,"gamma1"] + sims[,"gamma1(armexperimental)"], 
                                    sims[,"gamma2"] + sims[,"gamma2(armexperimental)"], 
                                    sims[,"gamma3"], 
                                    sims[,"gamma4"], 
                                    sims[,"gamma5"]),
                      offset = sims[,"armexperimental"],
                      knots = fit_spline$knots)
  
  hr <- h_e / h_c
  res_mat[i,] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

plot(1:109, res_mat[,2], type = "l", ylim = c(0.3,1.5), xlim = c(5, 109),
     xlab = "Time", ylab = "Hazard Ratio", main = "(D)")
points(1:109, res_mat[,1],type = "l", lty = 2)
points(1:109, res_mat[,3],type = "l", lty = 2)
abline(h = 0.7373, col = 1, lty = 3)

# Step 3: Run dev.off() to create the file!
dev.off()

# 2 cleopatra study 

load("data/CLEOPATRA_2A.rda")
dat <- CLEOPATRA_2A
survfit(Surv(time, event) ~ arm, data = dat) |> plot(col=1:2)
dat$arm <- factor(ifelse(dat$arm == "control", "control", "experimental"))


## assess ph

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")

# Step 1: Call the pdf command to start the plot
pdf(file = "figs/appendix_cleopatra.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plots with R code


par(mfrow = c(2,2))
plot(fit_spline, type = "survival", ci = FALSE, col = c(2,4), ylim = c(0, 1), xlab = "Time", ylab = "Survival", main = "(A)")
legend("bottomleft", c("Control", "Experimental"), col = c(2,4), lty = c(1,1))
plot(fit_spline, type = "cumhaz", ci = TRUE, col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(B)")
plot(fit_spline, type = "cumhaz", log = "xy", col = c(2,4), xlim = c(5, 70), xlab = "Time", ylab = "Cumulative Hazard", main = "(C)")

## conf int for hazard ratio
## simulate from MVN distribution
sims <- normboot.flexsurvreg(fit_spline, 
                             B = 1e5, 
                             raw = TRUE)

res_mat <- matrix(0, nrow = 70, ncol = 3)
for (i in 1:70){
  ## hazard
  h_c <-  hsurvspline(i, gamma = sims[,c("gamma0", 
                                         "gamma1", 
                                         "gamma2", 
                                         "gamma3", 
                                         "gamma4", 
                                         "gamma5")], 
                      knots = fit_spline$knots)
  
  ## hazard
  h_e <-  hsurvspline(i, 
                      gamma = cbind(sims[,"gamma0"], 
                                    sims[,"gamma1"] + sims[,"gamma1(armexperimental)"], 
                                    sims[,"gamma2"] + sims[,"gamma2(armexperimental)"], 
                                    sims[,"gamma3"], 
                                    sims[,"gamma4"], 
                                    sims[,"gamma5"]),
                      offset = sims[,"armexperimental"],
                      knots = fit_spline$knots)
  
  hr <- h_e / h_c
  res_mat[i,] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

plot(1:70, res_mat[,2], type = "l", ylim = c(0.3,1.5), xlim = c(5, 70),
     xlab = "Time", ylab = "Hazard Ratio", main = "(D)")
points(1:70, res_mat[,1],type = "l", lty = 2)
points(1:70, res_mat[,3],type = "l", lty = 2)

# Step 3: Run dev.off() to create the file!
dev.off()


#3 - leader study 
dat <- read.csv("data/leader_km.csv")
dat$arm <- factor(ifelse(dat$arm == 1, "control", "experimental"))


## assess ph

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")


# Step 1: Call the pdf command to start the plot
pdf(file = "figs/appendix_leader.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plots with R code

par(mfrow = c(2,2))
plot(fit_spline, type = "survival", ci = FALSE, col = c(2,4), ylim = c(0.8, 1), xlab = "Time", ylab = "Survival", main = "(A)")
legend("bottomleft", c("Control", "Experimental"), col = c(2,4), lty = c(1,1))
plot(fit_spline, type = "cumhaz", ci = TRUE, col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(B)")
plot(fit_spline, type = "cumhaz", log = "xy", xlim = c(5, 50), col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(C)")

## conf int for hazard ratio
## simulate from MVN distribution
sims <- normboot.flexsurvreg(fit_spline, 
                             B = 1e5, 
                             raw = TRUE)

res_mat <- matrix(0, nrow = 50, ncol = 3)
for (i in 1:50){
  ## hazard
  h_c <-  hsurvspline(i, gamma = sims[,c("gamma0", 
                                         "gamma1", 
                                         "gamma2", 
                                         "gamma3", 
                                         "gamma4", 
                                         "gamma5")], 
                      knots = fit_spline$knots)
  
  ## hazard
  h_e <-  hsurvspline(i, 
                      gamma = cbind(sims[,"gamma0"], 
                                    sims[,"gamma1"] + sims[,"gamma1(armexperimental)"], 
                                    sims[,"gamma2"] + sims[,"gamma2(armexperimental)"], 
                                    sims[,"gamma3"], 
                                    sims[,"gamma4"], 
                                    sims[,"gamma5"]),
                      offset = sims[,"armexperimental"],
                      knots = fit_spline$knots)
  
  hr <- h_e / h_c
  res_mat[i,] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}



plot(1:50, res_mat[,2], type = "l", ylim = c(0.5,1.2), xlim = c(1, 50),
     xlab = "Time", ylab = "Hazard Ratio", main = "(D)")
points(1:50, res_mat[,1],type = "l", lty = 2)
points(1:50, res_mat[,3],type = "l", lty = 2)
abline(h = 0.869, col = 1, lty = 3)

# Step 3: Run dev.off() to create the file!
dev.off()



### 4 - poplar study

# load data
dat <- read_csv("data/poplar.csv")


## assess ph

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")


# Step 1: Call the pdf command to start the plot
pdf(file = "figs/appendix_poplar.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 10) # The height of the plot in inches

# Step 2: Create the plots with R code


par(mfrow = c(2,2))
plot(fit_spline, type = "survival", ci = FALSE, col = c(2,4), ylim = c(0, 1), xlab = "Time", ylab = "Survival", main = "(A)")
legend("bottomleft", c("Control", "Experimental"), col = c(2,4), lty = c(1,1))
plot(fit_spline, type = "cumhaz", ci = TRUE, col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(B)")
plot(fit_spline, type = "cumhaz", log = "xy", xlim = c(1, 30), col = c(2,4), xlab = "Time", ylab = "Cumulative Hazard", main = "(C)")


## conf int for hazard ratio
## simulate from MVN distribution
sims <- normboot.flexsurvreg(fit_spline, 
                             B = 1e5, 
                             raw = TRUE)

res_mat <- matrix(0, nrow = 30, ncol = 3)
for (i in 1:30){
  ## hazard
  h_c <-  hsurvspline(i, gamma = sims[,c("gamma0", 
                                         "gamma1", 
                                         "gamma2", 
                                         "gamma3", 
                                         "gamma4", 
                                         "gamma5")], 
                      knots = fit_spline$knots)
  
  ## hazard
  h_e <-  hsurvspline(i, 
                      gamma = cbind(sims[,"gamma0"], 
                                    sims[,"gamma1"] + sims[,"gamma1(armexperimental)"], 
                                    sims[,"gamma2"] + sims[,"gamma2(armexperimental)"], 
                                    sims[,"gamma3"], 
                                    sims[,"gamma4"], 
                                    sims[,"gamma5"]),
                      offset = sims[,"armexperimental"],
                      knots = fit_spline$knots)
  
  hr <- h_e / h_c
  res_mat[i,] <- quantile(hr, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
}

plot(1:30, res_mat[,2], type = "l", ylim = c(0.1,1.2), xlim = c(1, 30),
     xlab = "Time", ylab = "Hazard Ratio", main = "(D)")
points(1:30, res_mat[,1],type = "l", lty = 2)
points(1:30, res_mat[,3],type = "l", lty = 2)
abline(h = 0.6752, col = 1, lty = 3)


# Step 3: Run dev.off() to create the file!
dev.off()


