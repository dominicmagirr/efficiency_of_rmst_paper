library(survival)
library(flexsurv)
library(ggsurvfit)
source("src/survRM2_fns.R")

dat <- read.csv("data/sustain_ipd.csv")
dat$arm <- factor(ifelse(dat$arm == 1, "control", "experimental"))

## plot

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")

par(mfrow = c(2,2))

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
     xlab = "Time", ylab = "Hazard Ratio", main = "(A)")
points(1:109, res_mat[,1],type = "l", lty = 2)
points(1:109, res_mat[,3],type = "l", lty = 2)
abline(h = 0.7373, col = 1, lty = 3)

###################################################################

###########################
load("data/CLEOPATRA_2A.rda")
dat <- CLEOPATRA_2A



dat$arm <- factor(ifelse(dat$arm == "control", "control", "experimental"))

## plot

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")

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
     xlab = "Time", ylab = "Hazard Ratio", main = "(B)")
points(1:70, res_mat[,1],type = "l", lty = 2)
points(1:70, res_mat[,3],type = "l", lty = 2)


abline(h = 0.6767, col = 1, lty = 3)

#############################################
library(survival)
source("src/survRM2_fns.R")
dat <- read.csv("data/leader_km.csv")

dat$arm <- factor(ifelse(dat$arm == 1, "control", "experimental"))


fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")
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
     xlab = "Time", ylab = "Hazard Ratio", main = "(C)")
points(1:50, res_mat[,1],type = "l", lty = 2)
points(1:50, res_mat[,3],type = "l", lty = 2)
abline(h = 0.869, col = 1, lty = 3)

#####################################

### get and plot data
library(dplyr)
template <- tempfile(fileext = ".xlsx")

httr::GET(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-018-0134-3/MediaObjects/41591_2018_134_MOESM3_ESM.xlsx", 
          httr::write_disk(template))

dat <- readxl::read_excel(template,sheet = 2) %>% 
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) %>%
  mutate(event = -1 * (OS.CNSR - 1),
         time = OS,
         arm = factor(ifelse(TRT01P == "Docetaxel", "control", "experimental"))) %>% 
  select(time, event, arm) %>%
  as.data.frame()


fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")
# conf int for hazard ratio
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


## plot



