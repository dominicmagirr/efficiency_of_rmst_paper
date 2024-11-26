# Load required libraries
library(survival)
library(readr)
library(dplyr)
library(ggsurvfit)
library(patchwork)

#-----------------------------------------------------------
# SUSTAIN dataset
#-----------------------------------------------------------
sustain_data <- read_csv("data/sustain_ipd.csv") |> 
  mutate(arm = factor(ifelse(arm == 1, "control", "experimental")))

# Plot SUSTAIN data
p_1 <- survfit2(Surv(time, event) ~ arm, data = sustain_data) |> 
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.position = c(0.2, 0.2)) +
  ggtitle(expression("(A) SUSTAIN-6:"~Z^"Cox"==2.38*"; "~Z^"RMST"==2.22)) +
  geom_vline(xintercept = 108, linetype = 2)

  

#-----------------------------------------------------------
# CLEOPATRA dataset
#-----------------------------------------------------------
cleopatra_data <- read.csv("data/CLEOPATRA_2A.csv") |> 
  mutate(arm = factor(ifelse(arm == "control", "control", "experimental")))

# Plot CLEOPATRA data
p_2 <- survfit2(Surv(time, event) ~ arm, data = cleopatra_data) |> 
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_risktable() + 
  theme(legend.position = "none")  + 
  ggtitle(expression("(B) CLEOPATRA:"~Z^"Cox"==3.77*"; "~Z^"RMST"==3.75)) +
  geom_vline(xintercept = 65, linetype = 2)

#-----------------------------------------------------------
# LEADER dataset
#-----------------------------------------------------------
leader_data <- read_csv("data/leader_km.csv")
leader_data <- leader_data %>% mutate(arm = factor(ifelse(arm == 1, "control", "experimental")))

# Plot LEADER data
p_3 <- survfit2(Surv(time, event) ~ arm, data = leader_data) |> 
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_risktable() +
  theme(legend.position = "none")  +
  ggtitle(expression("(C) LEADER:"~Z^"Cox"==2.53*"; "~Z^"RMST"==2.63)) +
  geom_vline(xintercept = 48, linetype = 2)


#-----------------------------------------------------------
# POPLAR dataset
#-----------------------------------------------------------
poplar_data <- read_csv("data/poplar.csv")

# Plot POPLAR data
p_4 <- survfit2(Surv(time, event) ~ arm, data = poplar_data) |> 
  ggsurvfit() +
  add_censor_mark() +
  add_confidence_interval() +
  add_risktable() + 
  theme(legend.position = "none")  + 
  ggtitle(expression("(D) POPLAR:"~Z^"Cox"==2.75*"; "~Z^"RMST"==2.24))+
  geom_vline(xintercept = 24, linetype = 2)


#-----------------------------------------------------------
# Combine plots into a single figure
#-----------------------------------------------------------
final_plot <- wrap_plots(ggsurvfit_build(p_1),
                      ggsurvfit_build(p_2),
                      ggsurvfit_build(p_3),
                      ggsurvfit_build(p_4)
                      )

#-----------------------------------------------------------
# save final_plot
#-----------------------------------------------------------
ggsave(
  "figs/km_plots_case_studies.pdf",
  final_plot,
  width = 12,
  height = 9,
  units = "in",
  dpi = 300
)
