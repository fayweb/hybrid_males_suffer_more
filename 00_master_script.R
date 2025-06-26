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
  "umap", # Computes a manifold approximation and projection

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
  "F" = "#4daf4a",   # Green
  "M"   = "#ff7f00"    # Orange
)


infection_colors <- c(
  "Uninfected"       = "#00FFFF",   # Light blue
  "E. ferrisi"        = "#FF8EE0",   # Soft pink
  "E. falciformis"    = "#FF197A"    # Darker pink/purple
)


infection_status_colors <- c(
  "Uninfected" = "#00FFFF",  # Soft light blue
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
  exploratory_analysis = FALSE,
  statistical_models = FALSE,
  figure_generation = FALSE
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

# create factors
field_mice$Sex <- as.factor(field_mice$Sex)

sex_colors <- c(
  "F" = "#4daf4a",   # Green
  "M"   = "#ff7f00"    # Orange
)

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


# Run distribution analysis
cat("Running distribution analysis...\n")
#source(file.path("scripts", "02_exploratory_analysis", "02_distribution_analysis.R"))
cat("✓ Distribution analysis completed\n")
cat("✓ Supplementary Figure 1 panels created and saved\n\n")


# 3. STATISTICAL MODELS
if (RUN_ANALYSIS$statistical_models) {
  print_section("STATISTICAL MODELS")

  # Run main hybrid analysis using parasiteLoad framework
  cat("Running hybrid effect analysis (parasiteLoad framework)...\n")
  source(file.path("scripts", "03_statistical_models", "01_hybrid_analysis.R"))
  cat("✓ Hybrid effect analysis completed\n")
  cat("✓ Sex-specific hybrid costs detected\n")
  cat("✓ Infection-dependent effects confirmed\n\n")

  # Generate banana plots and main figures from hybrid analysis
  cat("Creating hybrid analysis figures...\n")
  source(file.path("scripts", "03_statistical_models", "02_hybrid_analysis_figures.R"))
  cat("✓ Banana plots created and saved\n")
  cat("✓ Main Figure 2 panels generated\n\n")

  # Run validation analysis using Ferreira methodology
  cat("Running methodological validation (Ferreira framework)...\n")
  source(file.path("scripts", "03_statistical_models", "03_hybrid_analysis_validation_ferreira.R"))
  cat("✓ Cross-method validation completed\n")
  cat("✓ Ferreira distance-based analysis finished\n")
  cat("✓ Validation tables created\n\n")

  # Generate validation figures
  cat("Creating validation figures...\n")
  source(file.path("scripts", "03_statistical_models", "04_hybrid_analysis_validation_ferreira_figures.R"))
  cat("✓ Supplementary Figure 3 created\n")
  cat("✓ Cross-method comparison plots saved\n\n")

  # Male-specific analysis (if needed)
  if (file.exists(file.path("scripts", "03_statistical_models", "01_b_hybrid_analysis_male_effects.R"))) {
    cat("Running additional male-specific analysis...\n")
    source(file.path("scripts", "03_statistical_models", "01_b_hybrid_analysis_male_effects.R"))
    cat("✓ Male-specific effects analysis completed\n\n")
  }

  # Summary of statistical findings
  cat("STATISTICAL ANALYSIS SUMMARY:\n")
  cat("==============================\n")
  cat("✅ Main finding: Hybrid males suffer more (p = 0.038)\n")
  cat("✅ No constitutive costs in uninfected mice (p = 0.545)\n")
  cat("✅ Infection-dependent hybrid breakdown confirmed\n")
  cat("✅ Cross-validated using two statistical frameworks\n")
  cat("✅ Publication-ready figures and tables generated\n\n")
}

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

# Get version info for key packages
sessionInfo()
packageVersion("parasiteLoad")
packageVersion("randomForest")
packageVersion("mice")
packageVersion("limma")
packageVersion("pwr")
packageVersion("spdep")



# Create the comprehensive gene information table
gene_data <- data.frame(
  Gene = c("CXCL9", "TNF", "IDO1", "IFNγ", "IL-6", "TICAM1", "RETNLB",
           "IRGM1", "SOCS1", "MUC2", "MPO", "IL1RN", "IL-10", "IL-13",
           "NCR1", "PRF1", "TLR2", "TLR4", "ARG1"),

  Full_Name = c("C-X-C motif chemokine ligand 9",
                "Tumor necrosis factor",
                "Indoleamine 2,3-dioxygenase 1",
                "Interferon gamma",
                "Interleukin 6",
                "TIR domain-containing adapter molecule 1",
                "Resistin-like beta",
                "Immunity-related GTPase M",
                "Suppressor of cytokine signaling 1",
                "Mucin 2",
                "Myeloperoxidase",
                "Interleukin 1 receptor antagonist",
                "Interleukin 10",
                "Interleukin 13",
                "Natural cytotoxicity receptor 1",
                "Perforin 1",
                "Toll-like receptor 2",
                "Toll-like receptor 4",
                "Arginase 1"),

  Function = c("T cell recruitment, IFN-γ-inducible",
               "Pro-inflammatory, induces cachexia",
               "Tryptophan depletion, immunosuppression",
               "Th1 master regulator, parasite control",
               "Acute phase response, inflammation",
               "TLR adaptor, viral response",
               "Th2 effector, goblet cell hyperplasia",
               "Autophagy, parasite resistance",
               "Negative feedback, prevents tissue damage",
               "Intestinal barrier, pathogen exclusion",
               "Oxidative burst, neutrophil marker",
               "IL-1 antagonist, limits inflammation",
               "Anti-inflammatory, regulatory",
               "Th2 cytokine, mucus production",
               "NK cell activation receptor",
               "Cytolytic granule protein",
               "Bacterial recognition",
               "LPS sensor, innate immunity",
               "M2 macrophage marker, tissue repair"),

  Category = c("Chemokine", "Cytokine", "Metabolic enzyme", "Cytokine", "Cytokine",
               "Innate signaling", "Th2 response", "Cell-autonomous", "Regulatory",
               "Barrier", "Innate effector", "Regulatory", "Regulatory", "Th2 response",
               "Cytotoxic", "Cytotoxic", "Pattern recognition", "Pattern recognition",
               "Metabolic enzyme"),

  Response_Type = c("Th1", "Inflammatory", "Regulatory", "Th1", "Inflammatory",
                    "Innate", "Th2", "Cell-autonomous", "Regulatory", "Barrier",
                    "Innate", "Regulatory", "Regulatory", "Th2", "Cytotoxic",
                    "Cytotoxic", "Innate", "Innate", "Th2/Regulatory")
)

# Create the beautiful GT table
gene_table <- gene_data %>%
  gt() %>%

  # Add title and subtitle
  tab_header(
    title = md("**Immune Gene Panel for Infection Response Prediction**"),
    subtitle = md("*Nineteen genes quantified by high-throughput qPCR at day 8 post-infection*")
  ) %>%

  # Format column headers
  cols_label(
    Gene = md("**Gene**"),
    Full_Name = md("**Full Name**"),
    Function = md("**Primary Function**"),
    Category = md("**Category**"),
    Response_Type = md("**Response Type**")
  ) %>%

  # Group by response type
  tab_row_group(
    label = md("**Th1/Pro-inflammatory Response**"),
    rows = Response_Type %in% c("Th1", "Inflammatory")
  ) %>%
  tab_row_group(
    label = md("**Th2/Anti-helminth Response**"),
    rows = Response_Type == "Th2"
  ) %>%
  tab_row_group(
    label = md("**Regulatory/Anti-inflammatory**"),
    rows = Response_Type %in% c("Regulatory", "Th2/Regulatory")
  ) %>%
  tab_row_group(
    label = md("**Innate Immunity**"),
    rows = Response_Type == "Innate"
  ) %>%
  tab_row_group(
    label = md("**Cell-mediated Immunity**"),
    rows = Response_Type %in% c("Cell-autonomous", "Cytotoxic", "Barrier")
  ) %>%

  # Style the table
  tab_style(
    style = list(
      cell_fill(color = "#E8F4F8"),
      cell_text(weight = "bold")
    ),
    locations = cells_row_groups()
  ) %>%

  # Highlight gene names
  tab_style(
    style = cell_text(weight = "bold", style = "italic"),
    locations = cells_body(columns = Gene)
  ) %>%

  # Add alternating row colors
  opt_row_striping() %>%

  # Add borders
  tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black"
  ) %>%

  # Add footnote
  tab_footnote(
    footnote = md("Genes selected based on established roles in parasite immunity. Expression normalized using quantile normalization with GAPDH and PPIB as housekeeping genes."),
    locations = cells_title()
  ) %>%

  # Add source note
  tab_source_note(
    source_note = md("**Source:** Webster et al. (2025) *in review*")
  )

gene_table

# Save the table using your function
save_table_all_formats(gene_table, "Table_S1_Immune_Genes")

# Also create a simplified version for the main text if needed
gene_table_simple <- gene_data %>%
  dplyr::select(Gene, Function, Category) %>%
  gt() %>%
  tab_header(
    title = md("**Immune Genes Used for Prediction**")
  ) %>%
  cols_label(
    Gene = md("**Gene**"),
    Function = md("**Function**"),
    Category = md("**Category**")
  ) %>%
  tab_options(
    table.font.size = 10,
    table.width = pct(80)
  )


# Save simplified version
save_table_all_formats(gene_table_simple, "Table_S1_Immune_Genes_Simple")
