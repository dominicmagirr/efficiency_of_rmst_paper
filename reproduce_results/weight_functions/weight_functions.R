## plot weight functions for simulation scenarios...
plot_weight_f <- function(s_1 = 0.9,
                          rec_period = 0.001,
                          hr = 0.67,
                          max_t = 3.001,
                          tau_target = 3){
  
  
  event_rate <- dplyr::case_when(s_1 == 0.9 ~ "Low",
                                 s_1 == 0.6 ~ "Moderate",
                                 s_1 == 0.2 ~ "High",
                                 .default = as.character(s_1))
  
  rec_rate <- dplyr::case_when(rec_period == 0.001 ~ "Instant",
                               rec_period == 0.5 ~ "Fast",
                               rec_period == 1.5 ~ "Moderate",
                               rec_period == 2.5 ~ "Slow",
                               .default = as.character(s_1))
  
  
  lambda_e = -log(s_1) / 3
  lambda_c = lambda_e / hr
  
  t_seq <- seq(0, tau_target, length.out = 100)
  
  ## censoring distribution
  surv_c <- 1 - pmax(0, (t_seq - (max_t - rec_period)) / rec_period)
  
  ## control survival
  surv_0 <- exp(-lambda_c*t_seq)
  
  ## L-R weight
  w_lr <- surv_c * surv_0
  
  ## RMST weight
  
  w_rmst <- (exp(-lambda_c*t_seq) - exp(-lambda_c * tau_target)) / (1 - exp(-lambda_c * tau_target))
  
  
  out_dat <- data.frame(Time = rep(t_seq, 2),
                        Weight = c(w_lr / mean(w_lr), w_rmst / mean(w_rmst)),
                        Test = rep(c("PH", "RMST"), each = length(t_seq)),
                        event_rate = event_rate,
                        rec_rate = rec_rate)
  
  return(out_dat)
  
}
#########################################################
s_1 <- c(0.9, 0.6, 0.2)
rec_period <- c(0.001, 0.5, 1.5, 2.5)
p_res <- purrr::map2_df(expand.grid(s_1, rec_period)[,1],
                        expand.grid(s_1, rec_period)[,2],
                        plot_weight_f)
#########################################################
library(ggplot2)
p_res$rec_rate <- factor(p_res$rec_rate,
                           levels = c("Instant", "Fast", "Moderate", "Slow"),
                         labels = c("Instant recruitment", "Fast recruitment", "Moderate recruitment", "Slow recruitment"))

p_res$event_rate <- factor(p_res$event_rate,
                         levels = c("Low", "Moderate", "High"),
                         labels = c("Low event rate", "Moderate event rate", "High event rate"))


ggplot(data = p_res,
       mapping = aes(x = Time, y = Weight, linetype = Test)) +
  geom_line() +
  facet_wrap(event_rate ~ rec_rate, scales = "free") +
  theme_bw()












