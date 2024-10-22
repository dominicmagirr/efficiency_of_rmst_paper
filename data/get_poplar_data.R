# Load necessary libraries
library(dplyr)
library(httr)
library(readxl)
library(readr)

# Define the temporary file path for the downloaded Excel file
template <- tempfile(fileext = ".xlsx")

# Download the Poplar trial data from the specified URL
GET(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-018-0134-3/MediaObjects/41591_2018_134_MOESM3_ESM.xlsx", 
    write_disk(template, overwrite = TRUE))

# Read the second sheet of the Excel file and select relevant columns
dat <- read_excel(template, sheet = 2) |>
  select(PtID, ECOGGR, OS, OS.CNSR, TRT01P) |?
  # Create event and arm variables and select required columns
  mutate(
    event = -1 * (OS.CNSR - 1),
    time = OS,
    arm = factor(ifelse(TRT01P == "Docetaxel", "control", "experimental"))
  ) |>
  select(time, event, arm) |>
  as.data.frame()

# Write the cleaned data to a CSV file
write_csv(dat, "data/poplar.csv")
