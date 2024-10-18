# PH & RMST Efficiency

## Description

This repository contains code to reproduce all outputs for the paper titled **"Efficiency of nonparametric two-sample superiority tests based on restricted mean survival time under proportional hazards."** The code includes simulations, case studies, and necessary functions for data analysis.

## Authors

- Dominic Magirr (maintainer)
- Craig Wang
- Xinlei Deng
- Mark Baillie

## Repository Overview

The repository is organized into the following main directories:

### 1. `data/`
- **CLEOPATRA_2a.rda**: Data from the CLEOPATRA trial.
- **leader_km.csv**: Data for the leader case study.
- **sustain_idp.csv**: Data for the sustain case study.

This directory contains all the data sets used for the analyses in the paper.

### 2. `reproduce_results/`
This directory includes code to reproduce all the results presented in the paper. It is subdivided into the following folders:

- **case_studies/**: Code to reproduce the case studies, including analyses of clinical trial data.
- **simulation/**: Code to run the simulation studies discussed in the paper.
- **weight_functions/**: Code related to different weight functions used in the paper's analysis.

### 3. `src/`
This folder contains helper functions that are used throughout the analysis:
- **sim_helper_fns.R**: Helper functions for simulations.
- **survRM2_fns.R**: Functions for performing RMST calculations and related survival analysis.

## Installation

To run the analysis, make sure you have R installed with the necessary packages. You can install the required packages by running the following script:

```R
source("install_packages.R")
```

## How to Reproduce Results

1. Clone the repository to your local machine.
   
   ```bash
   git clone https://github.com/yourusername/ph-rmst-efficiency.git
   cd ph-rmst-efficiency
   ```

2. Open R or RStudio and set your working directory to the repository:

   ```R
   setwd("path_to_repository")
   ```

3. Load the required data:

   ```R
   load("data/CLEOPATRA_2a.rda")
   ```

4. Run the following scripts to reproduce specific sections:
   - **Case Studies**: 
     ```R
     source("reproduce_results/case_studies/run_case_studies.R")
     ```
   - **Simulation**: 
     ```R
     source("reproduce_results/simulation/run_simulations.R")
     ```
   - **Weight Functions**: 
     ```R
     source("reproduce_results/weight_functions/run_weight_functions.R")
     ```

## Citation

Please cite the following paper when using this code:

> Magirr, D., Wang, C., Deng, X., & Baillie, M. (2024). "Efficiency of nonparametric two-sample superiority tests based on restricted mean survival time under proportional hazards." *TBC*.

