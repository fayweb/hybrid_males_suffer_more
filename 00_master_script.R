# ==============================================================================
# MASTER SCRIPT: Hybrid Males Suffer More - Chapter 2 Analysis
# ==============================================================================
# Project: Tolerance of hybrid hosts against infections
# Chapter: 2 - Hybrid inflammation and sex-specific infection tolerance
# Author: Fay Webster
# Institution: Humboldt-University Berlin & Leibniz Institute for Zoo and Wildlife Research
#
# Description: Master script for Chapter 2 analysis examining how genetic
# admixture affects immune responses and tolerance to Eimeria infections in
# wild mice, with focus on sex-specific differences. Builds on Chapter 1
# validated Random Forest model and processed field data.
#
# Data Integration: Uses analysis-ready datasets from Chapter 1 pipeline
# - field_mice_complete.csv: 336 wild mice with predictions
# - chapter1_rf_model.rds: Validated Random Forest model
# - utility_functions.R: Custom analysis functions
# ==============================================================================

# Clear workspace
rm(list = ls())
gc()

# ==============================================================================
# PROJECT SETUP & CONFIGURATION
# ==============================================================================

# Set project root (automatically detects if using RStudio project)
if (rstudioapi::isAvailable()) {
  project_root <- dirname(rstudioapi::getActiveDocumentContext()$path)
} else {
  project_root <- getwd()
}

# Create directory structure if it doesn't exist
required_dirs <- c(
  "data/processed",
  "scripts/01_data_preparation",
  "scripts/02_exploratory_analysis",
  "scripts/03_statistical_models",
  "scripts/04_figure_generation",
  "scripts/05_supplementary",
  "results/figures",
  "results/tables",
  "results/supplementary"
)

for (dir in required_dirs) {
  if (!dir.exists(file.path(project_root, dir))) {
    dir.create(file.path(project_root, dir), recursive = TRUE)
  }
}

# ==============================================================================
# PACKAGE MANAGEMENT & LOADING
# ==============================================================================
required_packages <- c(
  # Data manipulation
  "tidyverse",      # ggplot2, dplyr, tidyr, readr, etc.
  "data.table",     # Fast data manipulation
  "janitor",        # Data cleaning

  # Statistical analysis
  "broom",          # Tidy statistical output
  "car",            # ANOVA, regression diagnostics
  "emmeans",        # Estimated marginal means
  "multcomp",       # Multiple comparisons
  "lme4",           # Mixed-effects models
  "performance",    # Model diagnostics
  "emmeans",

  # Multivariate analysis
  "FactoMineR",     # PCA analysis
  "factoextra",     # PCA visualization
  "corrplot",       # Correlation plots
  "vegan",          # Ecological statistics

  # Machine Learning (for applying Chapter 1 model)
  "randomForest",   # Random forest models
  "caret",          # Model training and validation

  # Visualization
  "ggplot2",        # Grammar of graphics
  "ggpubr",         # Publication-ready plots
  "patchwork",      # Combine plots
  "RColorBrewer",   # Color palettes
  "viridis",        # Color scales
  "scales",         # Scale functions
  "cowplot",        # Plot arrangements
  "ggridges",       # Ridge plots
  "ggbeeswarm",     # Bee swarm plots
  "ggsignif",       # Significance brackets
  "ggeffects",      # Effect plots
  "gt",
  "glue",
  "forcats",
  "ggrepel",
  "ggtext",
  "gridExtra",


  # File I/O and utilities
  "readr",          # Fast file reading
  "here",           # Path management
  "stringr",        # String manipulation
  "forcats",        # Factor manipulation

  # Hybrid analysis framework
  "parasiteLoad",    # Alice Balard's hybrid analysis framework (NO COMMA!)

  # Maps
  "sf",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggspatial",
  "viridis",
  "ggmap",
  "rosm"

)

# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Install and load all required packages
cat("Loading required packages...\n")

# Special handling for parasiteLoad dependencies
if (!require("optimx", character.only = TRUE)) {
  # Install specific version of optimx required for parasiteLoad
  install.packages("optimx", version = "2021-10.12", dependencies = TRUE)
  library("optimx", character.only = TRUE)
}

if (!require("parasiteLoad", character.only = TRUE)) {
  # Install from Alice Balard's GitHub if not available on CRAN
  if (!require("devtools", character.only = TRUE)) {
    install.packages("devtools")
    library("devtools")
  }
  devtools::install_github("alicebalard/parasiteLoad@v2.0", force = TRUE)
  library("parasiteLoad", character.only = TRUE)
}

install_and_load(required_packages)
cat("All packages loaded successfully!\n\n")

# ==============================================================================
# GLOBAL SETTINGS & PARAMETERS
# ==============================================================================

# Set global options
options(
  stringsAsFactors = FALSE,
  scipen = 999,  # Avoid scientific notation
  digits = 4
)

# Set random seed for reproducibility
set.seed(42)

# ggplot2 theme settings
theme_set(theme_classic() +
            theme(
              text = element_text(size = 12),
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 12),
              strip.text = element_text(size = 12)
            ))

# Color palettes for consistent visualization
hybrid_colors <- c(
  "M.m.domesticus" = "firebrick1",      # Red
  "Hybrid" = "purple",              # Blue
  "M.m.musculus" = "blue"         # Green
)

sex_colors <- c(
  "Female" = "#4daf4a",   # Green
  "Male"   = "#ff7f00"    # Orange
)


infection_colors <- c(
  "Uninfected"       = "#A6CEE3",   # Light blue
  "E. ferrisi"        = "#FF8EE0",   # Soft pink
  "E. falciformis"    = "#FF197A"    # Darker pink/purple
)


infection_status_colors <- c(
  "Uninfected" = "#A6CEE3",  # Soft light blue
  "Infected with Eimeria spp." = "#FF7094"  # Medium pink
)




infection_factors <- c("Uninfected",
                     "E. ferrisi",
                     "E. falciformis")

infection_presence_factors <- c("FALSE", "TRUE")

# Gene names (19 immune genes from Chapter 1)
immune_genes <- c(
  "IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1",
  "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB",
  "SOCS1", "TICAM1", "TNF"
)

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

# Define analysis parameters
ANALYSIS_PARAMS <- list(
  # PCA settings
  pca_scale = TRUE,
  pca_center = TRUE,
  n_components = 5,

  # Statistical significance threshold
  alpha = 0.05,

  # Multiple testing correction
  p_adjust_method = "fdr",

  # Sample size thresholds
  min_group_size = 5,

  # Figure dimensions (inches)
  fig_width = 8,
  fig_height = 6,
  fig_dpi = 300
)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

# Load custom functions from Chapter 1
source(file.path("scripts", "01_data_preparation", "utility_functions.R"))



# ==============================================================================
# ANALYSIS WORKFLOW CONTROL
# ==============================================================================

# Control which analyses to run (set to TRUE/FALSE as needed)
RUN_ANALYSIS <- list(
  data_loading = TRUE,
  exploratory_analysis = TRUE,
  statistical_models = TRUE,
  figure_generation = TRUE,
  supplementary_analysis = TRUE
)

# ==============================================================================
# MAIN ANALYSIS WORKFLOW
# ==============================================================================

cat("Starting Hybrid Males Suffer More Analysis Pipeline\n")
cat("Project Root:", project_root, "\n")
cat("Timestamp:", Sys.time(), "\n\n")

# 1. DATA LOADING & SETUP
if (RUN_ANALYSIS$data_loading) {
  print_section("DATA LOADING & SETUP")

  # Load primary dataset
  field_mice <- load_field_data()

  # Load Chapter 1 model
  rf_model <- load_chapter1_model()

  # Basic data summary
  cat("\nDataset Summary:\n")
  cat("- Sample size:", nrow(field_mice), "mice\n")
  cat("- Predicted weight loss range:",
      round(min(field_mice$predicted_weight_loss, na.rm = TRUE), 2), "to",
      round(max(field_mice$predicted_weight_loss, na.rm = TRUE), 2), "%\n")

  # Check for key variables
  key_vars <- c("HI", "Sex", "predicted_weight_loss")
  missing_vars <- key_vars[!key_vars %in% names(field_mice)]
  if (length(missing_vars) > 0) {
    warning("⚠ Missing key variables:", paste(missing_vars, collapse = ", "))
  } else {
    cat("✓ All key variables present for hybrid/sex analysis\n")
  }
}

# 2. EXPLORATORY ANALYSIS
if (RUN_ANALYSIS$exploratory_analysis) {
  print_section("EXPLORATORY ANALYSIS")

  # Run Figure 1 data exploration and overview
  cat("Running data exploration and Figure 1 generation...\n")
  source(file.path("scripts", "02_exploratory_analysis", "01_data_exploration.R"))
  cat("✓ Data exploration completed\n")
  cat("✓ Figure 1 panels created and saved\n")
  cat("✓ Statistical models for sex/infection effects completed\n")
  cat("✓ Summary tables generated\n\n")
}


# Add this after the exploratory analysis section

# Run distribution analysis
cat("Running distribution analysis...\n")
source(file.path("scripts", "02_exploratory_analysis", "02_distribution_analysis.R"))
cat("✓ Distribution analysis completed\n")
cat("✓ Supplementary Figure 1 panels created and saved\n\n")



# ==============================================================================
# ANALYSIS COMPLETION
# ==============================================================================

print_section("SETUP COMPLETE")
cat("Chapter 2 analysis environment initialized successfully!\n")
cat("Datasets loaded and ready for hybrid/sex analysis.\n\n")

cat("Next steps:\n")
cat("1. Run exploratory analysis to examine hybrid patterns\n")
cat("2. Analyze sex-specific differences in infection tolerance\n")
cat("3. Generate publication figures\n")
cat("4. Create manuscript tables\n\n")

cat("Key objects in environment:\n")
cat("- field_mice: Primary dataset (n =", nrow(field_mice), ")\n")
if (!is.null(rf_model)) {
  cat("- rf_model: Chapter 1 Random Forest model\n")
}
cat("- immune_genes: 19 gene names for analysis\n")
cat("- Color palettes: hybrid_colors, sex_colors, infection_colors\n\n")

# Save workspace for future reference
save.image(file.path("results", "chapter2_workspace.RData"))
cat("Workspace saved to results/chapter2_workspace.RData\n")

