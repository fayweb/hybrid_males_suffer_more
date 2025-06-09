# Hybrid Males Suffer More: Chapter 2 Analysis

**Author:** Fay Webster  
**Institution:** Humboldt-University Berlin & Leibniz Institute for Zoo and Wildlife Research  

## Project Overview

This repository contains the analysis code for **Chapter 2** of a cumulative PhD thesis, which investigates how genetic admixture affects immune responses and tolerance to *Eimeria* infections in wild mice, with a particular focus on **sex-specific differences** in hybrid populations.

**Chapter Integration:** This project builds on **Chapter 1** (immune-based health prediction) by applying the validated Random Forest model to examine hybrid and sex effects in natural mouse populations.

## Research Questions

1. **Do hybrid mice exhibit elevated inflammatory responses** compared to parental subspecies?
2. **Are there sex-specific differences in infection tolerance** among hybrid mice?
3. **How does genetic admixture affect the balance** between resistance and tolerance?
4. **Do hybrid males suffer disproportionately** from parasitic infections?

## Study System

- **Species:** Wild *Mus musculus* from the European hybrid zone
- **Subspecies:** *M. m. domesticus* (Western Europe) × *M. m. musculus* (Eastern Europe)
- **Pathogen:** *Eimeria* spp. (intracellular coccidia)
- **Sample Size:** 336 wild-caught mice with complete data
- **Approach:** Field-collected mice with hybrid indices, immune gene expression, infection loads, and predicted health outcomes

## Data Integration Strategy

### **Source:** Chapter 1 Repository (`Hybrid_health_outcomes`)
This project uses **processed, analysis-ready data** from the completed Chapter 1 pipeline:

#### **Integrated Files:**
1. **`field_mice_complete.csv`** - Primary dataset (336 wild mice)
   - Hybrid index scores (genetic admixture)
   - Sex information (male/female)
   - 19 immune gene expression profiles
   - *Eimeria* infection status and intensity
   - **Predicted weight loss** (from Chapter 1 Random Forest model)

2. **`chapter1_rf_model.rds`** - Validated Random Forest model
   - Trained on 136 laboratory mice
   - 47% variance explained in weight loss
   - CXCL9 as top predictor
   - Applied to predict health outcomes in wild mice

3. **`utility_functions.R`** - Custom analysis functions
   - Plotting helpers and statistical utilities
   - Visualization functions for publication figures

### **Why This Integration Approach:**
- ✅ **Leverages validated Chapter 1 model** rather than rebuilding from scratch
- ✅ **Focuses on natural populations** relevant to hybrid/sex questions
- ✅ **Uses proven, clean datasets** with quality control already applied
- ✅ **Maintains reproducibility** through saved model objects
- ✅ **Enables direct comparison** with Chapter 1 results

## Repository Structure

```
hybrid_males_suffer_more/
├── 00_master_script.R              # Main control script
├── README.md                       # This documentation
├── 
├── data/
│   └── processed/                  # Analysis-ready datasets
│       ├── field_mice_complete.csv     # Primary dataset (n=336)
│       ├── chapter1_rf_model.rds       # Trained RF model
│       └── [analysis outputs]
├── 
├── scripts/
│   ├── 01_data_preparation/        # Data loading and setup
│   │   └── utility_functions.R        # Custom functions from Chapter 1
│   ├── 02_exploratory_analysis/    # Descriptive statistics
│   ├── 03_statistical_models/      # PCA, linear models, predictions
│   ├── 04_figure_generation/       # Publication-ready figures
│   └── 05_supplementary/           # Additional analyses
├── 
├── results/
│   ├── figures/                    # Main publication figures
│   ├── tables/                     # Statistical output tables
│   └── supplementary/              # Supplementary materials
└── 
└── manuscript/
    ├── drafts/                     # Manuscript versions
    └── figures/                    # Final manuscript figures
```

## Getting Started

### Prerequisites
- R version 4.0 or higher
- RStudio (recommended)
- Required R packages (automatically installed by master script)

### Quick Start
1. **Clone this repository**
2. **Open `hybrid_males_suffer_more.Rproj` in RStudio**
3. **Run the master script:**
   ```r
   source("00_master_script.R")
   ```

The master script will:
- Load all required packages
- Import the processed datasets
- Run all analysis modules
- Generate figures and tables
- Save results to appropriate folders

## Key Analyses

### **1. Data Integration & Setup**
- Load processed field data with predictions
- Import Chapter 1 Random Forest model
- Quality control and data validation

### **2. Exploratory Analysis**
- Hybrid index distributions and sex ratios
- Infection prevalence by hybrid status and sex
- Immune gene expression patterns
- Predicted health outcome distributions

### **3. Statistical Modeling**
- **PCA:** Immune gene expression patterns by hybrid status
- **Linear models:** Hybrid and sex effects on inflammation
- **Interaction analysis:** Hybrid × sex effects on tolerance
- **Validation:** Model predictions vs. observed patterns

### **4. Key Hypotheses Testing**
- **H1:** Hybrid mice show elevated inflammatory responses
- **H2:** Males are more susceptible than females in hybrid populations
- **H3:** Genetic admixture disrupts optimal immune regulation
- **H4:** Sex-specific costs of hybridization affect infection tolerance

## Expected Outputs

### **Main Figures**
- **Figure 1:** Hybrid index distributions and study population
- **Figure 2:** PCA of immune gene expression by hybrid status and sex
- **Figure 3:** Predicted weight loss by hybrid index and sex
- **Figure 4:** Sex-specific inflammatory responses in hybrids
- **Figure 5:** Tolerance vs. resistance patterns across the hybrid zone

### **Statistical Tables**
- Population characteristics and sample sizes
- PCA loadings and variance explained
- Linear model results for hybrid and sex effects
- Post-hoc comparisons and effect sizes

## Dependencies

### Core Packages
- `tidyverse` - Data manipulation and visualization
- `FactoMineR`, `factoextra` - PCA analysis
- `randomForest` - Model application and validation
- `broom`, `car`, `emmeans` - Statistical modeling

See `00_master_script.R` for complete package list.

## Integration with Chapter 1

This analysis directly builds on:
- **Immune gene selection:** 19 genes validated in Chapter 1
- **Health prediction model:** Random Forest trained on lab infections
- **Methodological framework:** Cross-population validation approach
- **Biological interpretation:** CXCL9 and inflammatory pathway focus

## Citation

*[To be updated with publication information]*

Webster, F., Heitlinger, E., et al. (2024). Hybrid males suffer more: Sex-specific tolerance to parasitic infections in European house mouse populations. *[Journal TBD]*.

**Related Work:**
Webster, F., et al. (2024). Bridging controlled and natural systems: Immune markers allow assessment of parasite infection health effects in the wild. [Chapter 1 - submitted]

## Contact

**Fay Webster**  
Email: [fay.webster@hu-berlin.de]  
GitHub: [@fayweb](https://github.com/fayweb)

## Acknowledgments

- Research Training Group 2046 "Parasite Infections: From Experimental Models to Natural Systems" (DFG)
- Chapter 1 collaborators and data collection teams
- Wild mouse sampling and laboratory analysis teams

---

*Last updated: June 9, 2025*
