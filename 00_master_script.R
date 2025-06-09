# ==============================================================================
# MASTER SCRIPT: Hybrid Males Suffer More - Chapter 2 Analysis
# ==============================================================================
# Project: Tolerance of hybrid hosts against infections
# Author: Fay Webster
# Institution: Humboldt-University Berlin & Leibniz Institute for Zoo and Wildlife Research

# Clear workspace
rm(list = ls())
gc()

# Required packages
required_packages <- c(
  "tidyverse", "FactoMineR", "factoextra", "randomForest",
  "ggplot2", "patchwork", "broom", "car"
)

# Install and load packages
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Global settings
set.seed(42)
theme_set(theme_classic())

# Color palettes
hybrid_colors <- c("M.m.domesticus" = "#E31A1C", "Hybrid" = "#1F78B4", "M.m.musculus" = "#33A02C")
sex_colors <- c("Male" = "#FF7F00", "Female" = "#6A3D9A")

# Workflow control
RUN_ANALYSIS <- list(
  data_preparation = TRUE,
  exploratory_analysis = TRUE,
  statistical_models = TRUE,
  figure_generation = TRUE
)

cat("Master script loaded successfully!\n")
cat("Project ready for analysis.\n")
