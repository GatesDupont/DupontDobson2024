# Dupont and Dobson 2024

This repository contains the code and data for the analysis presented in Dupont and Dobson 2024.

## Project Overview

This project investigates the population dynamics of North American bird species across different ecological regions.

## Repository Structure

- **data/**: Contains all the data files used in the analysis.
  - `AVONET2_eBird.xlsx`: Dataset used only for aligning taxonomic naming conventions.
  - `BCR-results/`: Results from the Bird Conservation Region (BCR) analyses.
  - `BCR-species/`: Contains species data for each of the six BCRs.
  - `Stage2_Hackett_MCC_no_neg.tre`: Phylogenetic tree data.
  - `detectability-factors/`: Data on detectability adjustments and factors.
- **scripts/**: Contains all the R scripts used for data processing and analysis.
  - `00_Functions.R`: Utility functions for analysis.
  - `01_01_Select_species_Atlantic_Northern_Forest_.R` to `01_06_Select_species_Shortgrass_Prairie.R`: Scripts for selecting species in different BCRs.
  - `02_Model_density_v062724.R`: Script for density modeling.
  - `03_drivers_phylo_v2.R`: Script analyzing phylogenetic drivers.
  - `fit_model_with_fallback.R`: Script for fitting neg-bin and Pois models with fallback options.
  - `species_abundance_estimates.R`: Utility script for estimating species abundances.
  - `theory/`: Contains script for theoretical analysis.
    - `habitat-loss.R`: Script for theoretical analysis.
- **README.md**: This file.
- **DupontDobson2024.Rproj**: R project file.
- **.gitignore**: Git ignore file.

## Installation and Setup

- Setup: R and several R packages are required to run these scripts.
- Installation: git clone https://github.com/GatesDupont/DupontDobson2024.git
