# PH & RMST Efficiency

## Overview

This repository contains the code needed to reproduce all outputs presented in the paper titled **"Efficiency of Nonparametric Two-Sample Superiority Tests Based on Restricted Mean Survival Time versus The Log-Rank Test Under Proportional Hazards."** The repository includes simulations, case studies, and supporting functions for data analysis.

## Authors

- Dominic Magirr (Maintainer)
- Craig Wang
- Xinlei Deng
- Mark Baillie

## Repository Structure

The repository is organized into the following main directories:

### 1. `data/`
This directory contains all datasets utilized in the analyses presented in the paper:

- **CLEOPATRA_2a.rda**: Data from the CLEOPATRA trial.
- **leader_km.csv**: Data for the LEADER case study.
- **sustain_idp.csv**: Data for the SUSTAIN case study.
- **poplar.csv**: Data for the POPLAR case study.

### 2. `reproduce_results/`
This directory includes all the scripts required to reproduce the results in the paper. It contains the following subfolders:

- **case_studies/**: Code for analyzing clinical trial data presented as case studies.
- **censoring/**: Code for illustrating the recruitment rate and censoring distributions.
- **simulation/**: Scripts to perform simulation studies discussed in the paper.
- **weight_functions/**: Code for applying the weight functions used in the analyses.

### 3. `src/`
This folder includes helper functions used across various parts of the analysis:

- **sim_helper_fns.R**: Contains helper functions for simulations.
- **survRM2_fns.R**: Functions for calculating Restricted Mean Survival Time (RMST) and related survival analyses.

## Installation and Setup

To run the analysis, you need to have R installed, along with the necessary R packages. You can install the required packages and run the main analysis by executing the provided script:

```r
source("00_install_required_packages.R")
```

This script will install all the dependencies automatically.

## Reproducing the Results

To reproduce the results presented in the paper, follow these steps:

1. Clone the repository to your local machine:
   
   ```bash
   git clone https://github.com/yourusername/ph-rmst-efficiency.git
   cd ph-rmst-efficiency
   ```

2. Install the required R packages by running the following command within your R session:

   ```r
   source("00_install_required_packages.R")
   ```

3. Run the main script to generate all figures and results. Note that for full simulation studies, you will need to toggle the relevant simulation flag:

   ```r
   source("01_generate_figures.R")
   ```

   The script includes both figure generation and case study analyses. Make sure to configure the settings as needed before running simulations, as these can be computationally intensive.

## Citation

If you use this code or data in your research, please cite the following paper:

> Magirr, D., Wang, C., Deng, X., & Baillie, M. (2024). *Efficiency of Nonparametric Two-Sample Superiority Tests Based on Restricted Mean Survival Time versus The Log-Rank Test Under Proportional Hazards.* (Journal To Be Confirmed).

## Notes and Considerations

- **Computational Requirements**: Some scripts, particularly the simulation scripts, may require significant computational resources. We recommend running these scripts on a machine with adequate memory and processing power.


## Contact

For questions or issues related to this repository, please contact the maintainer.

