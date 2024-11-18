#------------------------------------------------
# load required packages
#------------------------------------------------
library(survival)
library(nphRCT)
library(dplyr)
library(ggplot2)
library(readr)

#------------------------------------------------
# Utility function to process data and generate scores
#------------------------------------------------

# This function takes the input dataset, calculates Log-rank and KM RMST scores,
# and returns a combined dataframe with the standardized scores along with event types.
process_data <- function(data, tau, example_label) {
  # Calculate Log-rank scores
  df_lr <- find_scores(formula = Surv(time, event) ~ arm, 
                       data = data, 
                       method = "lr")$df |>
    mutate(label = "Log-rank", time = t_j) |>
    select(time, event, group, standardized_score, label)
  
  # Calculate Restricted Mean Survival Time (RMST) scores
  df_rmst <- 
    find_scores(formula = Surv(time, event) ~ arm, tau = tau, 
                data = data, 
                method = "rmst")$df |>
    mutate(label = "KM RMST", time = t_j) |>
    select(time, event, group, standardized_score, label)
  
  # Stack data frames for Log-rank and RMST scores
  # Add example label to identify the data source in the final plot
  rbind(df_lr, df_rmst) |>
    mutate(type = ifelse(event == 1, "Event", "Censored"),
           example = example_label)
}

#------------------------------------------------
# Load and process datasets
#------------------------------------------------

# List of data files with corresponding tau values and labels to be used in the analysis
data_files <- list(
  list(file = "data/sustain_ipd.csv", tau = 108, label = "(A)"),
  list(file = "data/CLEOPATRA_2A.csv", tau = 65, label = "(B)"),
  list(file = "data/leader_km.csv", tau = 48, label = "(C)"),
  list(file = "data/poplar.csv", tau = 24, label = "(D)")
)

# Initialize combined dataset
rmst_comp <- NULL

#------------------------------------------------
# Loop through datasets and process each
#------------------------------------------------
# This loop loads each dataset, applies necessary transformations, and calculates scores
# The resulting processed datasets are combined into a single dataset
datasets <- lapply(data_files, function(df_info) {
  # Load data from file (either .csv or .rda)
  if (grepl(".rda$", df_info$file)) {
    load(df_info$file)
    data <- get(ls()[ls() != "df_info"])
  } else {
    data <- read_csv(df_info$file)
  }
  
  # Convert arm variable to a factor with levels "control" and "experimental"
  data$arm <- factor(ifelse(data$arm == 1 | data$arm == "control", "control", "experimental"))
  
  # Process data to generate scores for both Log-rank and RMST
  process_data(data, df_info$tau, df_info$label)
})

# Combine all processed datasets into a single dataframe
rmst_comp <- do.call(rbind, datasets)

# Re-order factor levels for labels to ensure consistent ordering in the plot
rmst_comp$label <- factor(rmst_comp$label, levels = unique(rmst_comp$label))

#------------------------------------------------
# Generate the plot
#------------------------------------------------

# Create a ggplot object to visualize the standardized scores for each dataset
p_rmst_comp <- ggplot() +
  geom_point(data = rmst_comp,
             aes(
               x = time,
               y = standardized_score,
               color = group,
               alpha = type
             )) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  ylab("Score") +
  labs(color = "Arm", alpha = "Observation type") +
  # Add horizontal lines indicating the mean score for each group and label
  geom_hline(
    data = rmst_comp |> 
      group_by(group, label, example) |> 
      summarize(mean_score = mean(standardized_score)),
    aes(yintercept = mean_score, colour = group),
    linetype = 2
  ) +
  facet_wrap(label ~ example, scales = "free_x", nrow = 2)

#------------------------------------------------
# Save the final figure to a file
#------------------------------------------------
ggsave(
  "figs/std_scores_case_studies.pdf",
  plot = p_rmst_comp,
  width = 12,
  height = 9,
  units = "in",
  dpi = 300
)
