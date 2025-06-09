# Hybrid Health Outcomes Analysis

This repository contains the analysis code for investigating the impact of genetic hybridization on health outcomes in wild house mice (*Mus musculus*) using immune-based predictive models.

## Repository Structure

- `data/`: Raw and processed data files
- `R/`: Core function files
- `analysis/`: Analysis pipeline scripts
- `outputs/`: Generated figures, tables, and model objects
- `manuscript/`: Publication-ready outputs
- `tests/`: Function validation tests

## Quick Start

1. Load required packages and functions:
```r
source("R/01_data_preparation.R")
source("R/02_hybrid_analysis_functions.R") 
source("R/03_visualization_functions.R")
```

2. Run main analysis:
```r
source("analysis/01_main_hybrid_analysis.R")
```

## Dependencies

- R (>= 4.0.0)
- Alice Balard's parasiteLoad package
- tidyverse, ggplot2, patchwork

## Citation

[Your citation information here]
