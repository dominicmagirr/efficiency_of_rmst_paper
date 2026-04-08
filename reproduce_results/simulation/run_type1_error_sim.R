## Type I error simulation under H0 (HR = 1)
## Runs all 12 scenarios from Table 1 with equal survival distributions
## to verify both tests maintain nominal significance level.
##
## This addresses Reviewer 2 comment on Type I error control
## and the Associate Editor's request to reassure that tests are level.

library(survival)
source("src/sim_helper_fns.R")
source("src/survRM2_fns.R")

## ---- Post-processing function for Type I error ----
post_sim_t1e <- function(res, label = "") {

  # Two-sided rejection at alpha = 0.05
  rej_rmst_2sided <- mean(abs(res$z_rmst) > qnorm(0.975))

  # One-sided rejection at alpha = 0.025 (matching power analysis direction)
  rej_rmst_1sided <- mean(res$z_rmst > qnorm(0.975))
  rej_ph_1sided   <- mean(res$z_cox < qnorm(0.025))

  # Two-sided for PH
  rej_ph_2sided <- mean(abs(res$z_cox) > qnorm(0.975))

  n_sims <- nrow(res)
  mc_se <- sqrt(0.05 * 0.95 / n_sims)

  return(data.frame(
    scenario    = label,
    n_sims      = n_sims,
    rej_rmst_2s = round(rej_rmst_2sided, 4),
    rej_ph_2s   = round(rej_ph_2sided, 4),
    rej_rmst_1s = round(rej_rmst_1sided, 4),
    rej_ph_1s   = round(rej_ph_1sided, 4),
    mc_se       = round(mc_se, 4)
  ))
}

## ---- Scenario definitions ----
## Same scenarios as Table 1 but with HR = 1 (null hypothesis)
## s_0 values are the control arm survival at 3 years
## Under H0, both arms have the same survival distribution

scenarios <- data.frame(
  label      = c("1: Low/Instant",  "2: Low/Fast",  "3: Low/Moderate",  "4: Low/Slow",
                 "5: Mod/Instant",  "6: Mod/Fast",  "7: Mod/Moderate",  "8: Mod/Slow",
                 "9: High/Instant", "10: High/Fast", "11: High/Moderate", "12: High/Slow"),
  n_per_arm  = c(rep(1000, 4), rep(250, 4), rep(150, 4)),
  s1_3       = c(rep(0.9, 4), rep(0.6, 4), rep(0.2, 4)),
  rec_period = rep(c(0.0001, 0.5, 1.5, 2.5), 3)
)

# Under H0 with HR=1, s_0 = s_1 = s1_3 (no treatment effect)
# In the original code, s_0 = exp(log(s1_3) * 1.5) was the CONTROL arm survival
# because HR = 1/1.5 = 0.67, so control is worse: s_0 = s1^(1/HR) = s1^1.5
# Under H0 (HR=1), both arms have the same distribution.
# We use s_0 = s1_3 directly (both arms identical).

## ---- Run simulations ----
N_runs <- 10000
set.seed(42026)

cat("Running Type I error simulations (HR = 1) across", nrow(scenarios), "scenarios\n")
cat("N_runs per scenario:", N_runs, "\n")
cat("Expected MC SE for alpha=0.05:", round(sqrt(0.05 * 0.95 / N_runs), 4), "\n\n")

results <- list()

for (i in seq_len(nrow(scenarios))) {

  sc <- scenarios[i, ]
  cat("Scenario", sc$label, "... ")

  t_start <- Sys.time()

  res_i <- purrr::map_df(
    seq_len(N_runs),
    ~ one_sim_zs(
      n_per_arm  = sc$n_per_arm,
      s_0        = sc$s1_3,          # Under H0, control survival = experimental survival
      hr         = 1,                # NULL HYPOTHESIS
      rec_period = sc$rec_period,
      max_t      = 3.000
    )
  )

  t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
  cat(sprintf("done (%.1f min)\n", as.numeric(t_elapsed)))

  results[[i]] <- post_sim_t1e(res_i, label = sc$label)
}

## ---- Compile and display results ----
results_df <- do.call(rbind, results)

cat("\n========================================\n")
cat("Type I Error Results (H0: HR = 1)\n")
cat("========================================\n")
cat(sprintf("Monte Carlo SE for alpha=0.05: %.4f\n", results_df$mc_se[1]))
cat(sprintf("Expected range: [%.4f, %.4f]\n\n",
            0.05 - 2*results_df$mc_se[1],
            0.05 + 2*results_df$mc_se[1]))

print(results_df[, c("scenario", "rej_rmst_2s", "rej_ph_2s")], row.names = FALSE)

## ---- Save results ----
write.csv(results_df,
          file = "reproduce_results/simulation/type1_error_results.csv",
          row.names = FALSE)

cat("\nResults saved to reproduce_results/simulation/type1_error_results.csv\n")
