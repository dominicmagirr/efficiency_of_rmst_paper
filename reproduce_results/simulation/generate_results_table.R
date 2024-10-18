## Create Publication-Ready Table
library(dplyr)
library(tidyr)
library(readr)
library(gt)

# read in simulations results csv file and store as final_results
# change if simulation results change
final_results <- read_csv("results/simulation_results.csv")

final_results_wide <- final_results |>
  pivot_wider(
    id_cols = c(event_rate, recruitment_speed),
    names_from = max_t,
    values_from = c(pow_rmst, pow_ph, eff_rmst_cox, frac_post_tau),
    names_glue = "{.value}{ifelse(max_t == 4.5, '_plus', '')}"
  )


publication_table <- final_results_wide |>
  select(-frac_post_tau) |>
  gt(groupname_col = 'event_rate') |>
  tab_header(
    title = "Simulation Results",
    subtitle = "Power and Relative Efficiency Across Different Scenarios"
  ) |>
  cols_label(
    event_rate = "Event Rate",
    recruitment_speed = "Recruitment",
    pow_rmst = "Power (RMST, t_H)",
    pow_rmst_plus = "Power (RMST, t_H+)",
    pow_ph = "Power (PH, t_H)",
    pow_ph_plus = "Power (PH, t_H+)",
    eff_rmst_cox = "Relative Efficiency (t_H)",
    eff_rmst_cox_plus = "Relative Efficiency (t_H+)",
    frac_post_tau_plus = "% Events > τ (t_H+)"
  ) |>
  fmt_number(
    columns = c(pow_rmst, pow_rmst_plus, pow_ph, pow_ph_plus, eff_rmst_cox, eff_rmst_cox_plus),
    decimals = 2
  ) |>
  fmt_percent(
    columns = frac_post_tau_plus,
    decimals = 1
  ) |>
  tab_spanner(
    label = "Follow-up Time t_H",
    columns = c(pow_rmst, pow_ph, eff_rmst_cox)
  ) |>
  tab_spanner(
    label = "Follow-up Time t_H+",
    columns = c(pow_rmst_plus, pow_ph_plus, eff_rmst_cox_plus, frac_post_tau_plus)
  )

# Print the publication-ready table
print(publication_table)
