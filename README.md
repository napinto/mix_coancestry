# Navigating Population-level Coancestry in Forensic Mixture Analysis - Simulation and Likelihood Ratio Analysis

## Overview

This repository contains R scripts developed to simulate DNA mixture profiles and calculate likelihood ratios (LRs) used in forensic genetics relatedness studies. The work supports the methodology and analysis described in the manuscript "Navigating Population-level Coancestry in Forensic Mixture Analysis", Costa C et al. submitted to Scientif Reports, September 2025

The scripts simulate allelic profiles based on population allele frequencies, compute LRs per genetic marker for mixtures with two or three contributors, and summarize the results for manuscript figures and tables. Population substructure is modeled through a theta correction.

---

## Repository Structure

- `2cont_auxiliary_funcs.R`: Functions for LR calculations in two-contributor mixture scenarios.
- `2cont_main.R`: Main driver script to run two-contributor simulations and LR calculations.
- `3cont_auxiliary_funcs.R`: Functions extending LR calculations to three or more contributors.
- `3cont_main.R`: Main script for simulations and LR calculations on three-contributor mixtures.
- `FREQ_ALELICAS.csv`: Allele frequency data by marker (must be provided).
- `Navigating Population-level Coancestry in Forensic Mixture Analysis.pdf`: Manuscript reporting the methodology and findings, once published.

---

## Getting Started

### Prerequisites

- R (version ≥ 4.0 recommended)
- Packages: `dplyr`, `readr`, `purrr`, `openxlsx`

Install required packages with:

```install.packages(c("dplyr","readr","purrr","openxlsx"))```


### Allele Frequency Data

Prepare the `FREQ_ALELICAS.csv` file containing allele frequencies:

- First column: Allele identifiers
- Other columns: Markers with allele frequencies
- Frequencies must sum to 1 per marker

---

## Running the Simulations

### For Two Contributors:

```
source("2cont_auxiliary_funcs.R")
source("2cont_main.R")
```


### For Three or More Contributors:

```
source("3cont_auxiliary_funcs.R")
source("3cont_main.R")
```


---

## Output and Results

The scripts produce:

- Simulated allele profiles per marker
- Likelihood ratios computed with theta adjustment (0 and 0.01)
- Summary statistics per marker and mixture profile
- Excel files with LR results and summaries for manuscript use

---

## Notes

- Adjust the simulation size (variable `nsims`) in main scripts for balancing performance and precision.
- Seed initialization ensures reproducibility. If not set, it will use one derived from the marker name.
- Modify `cuts` parameters in main scripts to refine likelihood ratio grouping.

---

## Citation

Please cite the accompanying manuscript when using this data or methodology.

---

## Contact

For questions or contributions, please contact the repository maintainer at npinto@i3s.up.pt.

---

*Thank you for your interest in this project!*
