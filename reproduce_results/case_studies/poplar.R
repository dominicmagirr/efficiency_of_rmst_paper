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


## plot

p_4 <- survfit2(Surv(time, event) ~ arm, data = dat) |> 
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_quantile() +
  add_risktable()

patchwork::wrap_plots(ggsurvfit_build(p_1 + theme(legend.position = c(0.2, 0.2)) + ggtitle("(A)")),
                      ggsurvfit_build(p_2 + theme(legend.position = "none")  + ggtitle("(B)")),
                      ggsurvfit_build(p_3 + theme(legend.position = "none")  + ggtitle("(C)")),
                      ggsurvfit_build(p_4 + theme(legend.position = "none")  + ggtitle("(D)")))


## assess ph

fit_spline <- flexsurvspline(Surv(time, event) ~ arm +
                               gamma1(arm) +
                               gamma2(arm),
                             data = dat, k = 4, scale = "hazard")

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

## Z statistics
coxph(Surv(time, event) ~ arm, data = dat) |> summary()
dat$arm2 <- ifelse(dat$arm == "control", 0, 1)


res_rmst <- rmst2(time = dat$time,
                  status = dat$event,
                  arm = dat$arm2,
                  tau = 24)


rmst_est <- res_rmst$unadjusted.result[1,1]
rmst_se <- (res_rmst$unadjusted.result[1,2:3] |> diff()) / (2 * qnorm(0.975))
z_rmst <- rmst_est / rmst_se

rmst_est
z_rmst
sum(dat$time > 24 & dat$event == 1) / sum(dat$event)

####################################################
## score plots

library(nphRCT)
df_lr <- find_scores(formula=Surv(time, event) ~ arm, data=dat, method = "lr")

df_lr <- df_lr$df %>% 
  mutate(label = "Log-rank", time = t_j) %>% 
  dplyr::select(c(time, event, group, standardized_score, label))

df_rmst <- find_scores(formula=Surv(time, event) ~ arm, tau=24, data=dat, method = "rmst")

df_rmst <- df_rmst$df %>% 
  mutate(label = "KM RMST", time = t_j) %>% 
  dplyr::select(c(time, event, group, standardized_score, label))


# stack data frames
rmst_comp = rbind(df_lr,
                  df_rmst) %>%
  mutate(type = ifelse(event == 1, "Event", "Censored"),
         example = "(D)")


# re-order
rmst_comp$label = factor(rmst_comp$label,
                         levels = unique(rmst_comp$label))

#plot
p_rmst_comp = ggplot() +
  geom_point(data = rmst_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = rmst_comp %>% group_by(group, label, example) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(label ~ example, scales = "free_x", nrow = 2)

p_rmst_comp


rmst_comp_D <- rmst_comp

rmst_comp <- rbind(rmst_comp_A, rmst_comp_B, rmst_comp_C, rmst_comp_D)


# plot together
p_rmst_comp = ggplot() +
  geom_point(data = rmst_comp,  aes(x = time, y = standardized_score, color = group, alpha = type)) + 
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() + 
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  geom_hline(data = rmst_comp %>% group_by(group, label, example) %>% dplyr::summarize(mean_score = mean(standardized_score)), 
             aes(yintercept = mean_score, colour = group), linetype = 2) +
  facet_wrap(label ~ example, scales = "free_x", nrow = 2)

p_rmst_comp

