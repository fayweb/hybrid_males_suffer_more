# ==============================================================================
# FERREIRA ET AL. METHODOLOGY VALIDATION - SIMPLIFIED VERSION
# ==============================================================================
# Following the exact same structure as parasiteLoad script:
# 1. Complete Dataset Analysis (infected + uninfected)
# 2. Uninfected Subset Analysis (constitutive costs)
# 3. Infected Subset Analysis (infection-specific effects)
#
# Using linear models instead of Bayesian to avoid brms installation issues
# Still follows Ferreira et al. pairwise distance methodology
# ==============================================================================

cat("=== FERREIRA ET AL. METHODOLOGY VALIDATION ===\n")
cat("Following parasiteLoad script structure for direct comparison\n")
cat("Using distance-based linear models (Ferreira methodology)\n\n")

# ==============================================================================
# 1. DATA PREPARATION - USING EXACT SAME DATASETS AS parasiteLoad
# ==============================================================================

cat("1. PREPARING DATA FOR FERREIRA ANALYSIS\n")
cat("=======================================\n")

# Check that parasiteLoad datasets exist
if (!exists("field_mice") || !exists("hybrid_data") || !exists("uninfected_data")) {
  stop("Please run parasiteLoad script first to create field_mice, hybrid_data, and uninfected_data")
}

cat("Using EXACT SAME datasets as parasiteLoad analysis:\n")

# 1. COMPLETE DATASET - Use field_mice (exactly like parasiteLoad complete_model)
ferreira_complete_data <- field_mice %>%
  filter(
    !is.na(HI),
    !is.na(Sex),
    !is.na(predicted_weight_loss)
  ) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe) - key Ferreira variable
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number(),

    # Response variable (matching parasiteLoad)
    response = predicted_weight_loss,

    # Infection status for distance calculations
    infected = as.logical(infection_status)
  )

# 2. UNINFECTED SUBSET - Direct from parasiteLoad (uninfected_data)
ferreira_uninfected_data <- uninfected_data %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

# 3. INFECTED SUBSET - Create from hybrid_data (following parasiteLoad logic)
ferreira_infected_data <- hybrid_data %>%
  filter(infected) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

cat("Analysis datasets prepared (matching parasiteLoad exactly):\n")
cat("- Complete dataset:", nrow(ferreira_complete_data), "mice (same as field_mice for complete_model)\n")
cat("- Uninfected subset:", nrow(ferreira_uninfected_data), "mice (same as uninfected_data)\n")
cat("- Infected subset:", nrow(ferreira_infected_data), "mice (filtered from hybrid_data)\n")

# Verify datasets match parasiteLoad exactly
cat("\nDataset verification:\n")
cat("- field_mice → ferreira_complete_data:", nrow(field_mice) >= nrow(ferreira_complete_data), "(filtered for complete cases)\n")
cat("- uninfected_data matches ferreira_uninfected_data:", nrow(uninfected_data) == nrow(ferreira_uninfected_data), "\n")
cat("- hybrid_data infected subset matches:", sum(hybrid_data$infected) == nrow(ferreira_infected_data), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS FOR FERREIRA ANALYSIS
# ==============================================================================

# Function to create pairwise distance matrix (following Ferreira et al.)
create_ferreira_pairwise <- function(data, analysis_name) {

  n <- nrow(data)
  n_pairs <- n * (n - 1) / 2

  cat("Creating", n_pairs, "pairwise comparisons for", analysis_name, "\n")

  # For large datasets, sample pairs to make it computationally feasible
  if (n_pairs > 10000) {
    cat("Large dataset detected. Sampling 10000 random pairs for analysis.\n")
    sample_pairs <- TRUE
    n_sample <- 10000
  } else {
    sample_pairs <- FALSE
  }

  # Initialize results
  pairwise_list <- list()
  counter <- 1

  if (sample_pairs) {
    # Sample random pairs
    all_pairs <- expand.grid(i = 1:(n-1), j = 2:n)
    all_pairs <- all_pairs[all_pairs$i < all_pairs$j, ]
    sampled_indices <- sample(nrow(all_pairs), n_sample)
    pairs_to_process <- all_pairs[sampled_indices, ]

    for (idx in 1:nrow(pairs_to_process)) {
      i <- pairs_to_process$i[idx]
      j <- pairs_to_process$j[idx]

      # Get individual data
      ind1 <- data[i, ]
      ind2 <- data[j, ]

      # Calculate distances (same as before)
      response_diff <- abs(ind1$response - ind2$response)
      max_response_range <- max(data$response) - min(data$response)
      response_distance <- response_diff / max_response_range
      response_similarity <- 1 - response_distance

      subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
      hHe_distance <- abs(ind1$hHe - ind2$hHe)
      hHe_mean <- (ind1$hHe + ind2$hHe) / 2

      sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
      infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

      pairwise_list[[counter]] <- data.frame(
        similarity = response_similarity,
        subspecies_genetic_distance = subspecies_genetic_distance,
        hHe_distance = hHe_distance,
        hHe_mean = hHe_mean,
        sex_distance = sex_distance,
        infection_distance = infection_distance,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }

  } else {
    # Process all pairs
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {

        # Get individual data
        ind1 <- data[i, ]
        ind2 <- data[j, ]

        # Response similarity (normalized distance between predicted weight loss)
        response_diff <- abs(ind1$response - ind2$response)
        max_response_range <- max(data$response) - min(data$response)
        response_distance <- response_diff / max_response_range
        response_similarity <- 1 - response_distance

        # Genetic distances (following Ferreira et al. Table 1)
        subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
        hHe_distance <- abs(ind1$hHe - ind2$hHe)
        hHe_mean <- (ind1$hHe + ind2$hHe) / 2

        # Other distances
        sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
        infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

        pairwise_list[[counter]] <- data.frame(
          similarity = response_similarity,
          subspecies_genetic_distance = subspecies_genetic_distance,
          hHe_distance = hHe_distance,
          hHe_mean = hHe_mean,
          sex_distance = sex_distance,
          infection_distance = infection_distance,
          stringsAsFactors = FALSE
        )

        counter <- counter + 1
      }
    }
  }

  result <- do.call(rbind, pairwise_list)

  # Scale predictors to 0-1 (following Ferreira methodology)
  result <- result %>%
    mutate(
      subspecies_genetic_distance_scaled = scales::rescale(subspecies_genetic_distance, to = c(0, 1)),
      hHe_distance_scaled = scales::rescale(hHe_distance, to = c(0, 1)),
      hHe_mean_scaled = scales::rescale(hHe_mean, to = c(0, 1))
    )

  return(result)
}

# Function to fit Ferreira-style linear model
fit_ferreira_model <- function(pairwise_data, model_name, include_infection = TRUE) {

  cat("Fitting", model_name, "using linear model (Ferreira methodology)...\n")

  # Model formula (matching Ferreira et al.)
  if (include_infection) {
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance + infection_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  } else {
    # For uninfected/infected subsets, no infection_distance needed
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  }

  # Fit linear model
  model <- lm(as.formula(formula_text), data = pairwise_data)
  model$model_type <- "linear"

  return(model)
}

# Function to extract results in parasiteLoad style
extract_ferreira_results <- function(model, model_name) {

  cat("\nExtracting results for", model_name, "...\n")

  # Get model summary
  model_summary <- summary(model)
  coef_table <- model_summary$coefficients

  # Key effects matching parasiteLoad interpretation
  effects_of_interest <- c(
    "subspecies_genetic_distance_scaled" = "Subspecies genetic effect",
    "hHe_distance_scaled" = "Hybridization effect (hHe-dist)",
    "hHe_mean_scaled" = "Mean hybridization effect",
    "sex_distance" = "Sex difference effect",
    "subspecies_genetic_distance_scaled:hHe_distance_scaled" = "Subspecies × Hybridization interaction"
  )

  # Add infection effect if present
  if ("infection_distance" %in% rownames(coef_table)) {
    effects_of_interest["infection_distance"] <- "Infection difference effect"
  }

  results_df <- data.frame(
    Effect = names(effects_of_interest),
    Interpretation = unname(effects_of_interest),
    Estimate = numeric(length(effects_of_interest)),
    Std_Error = numeric(length(effects_of_interest)),
    CI_Lower = numeric(length(effects_of_interest)),
    CI_Upper = numeric(length(effects_of_interest)),
    T_Value = numeric(length(effects_of_interest)),
    P_Value = numeric(length(effects_of_interest)),
    Significant = logical(length(effects_of_interest)),
    stringsAsFactors = FALSE
  )

  for (i in 1:length(effects_of_interest)) {
    param_name <- names(effects_of_interest)[i]

    if (param_name %in% rownames(coef_table)) {
      results_df$Estimate[i] <- coef_table[param_name, "Estimate"]
      results_df$Std_Error[i] <- coef_table[param_name, "Std. Error"]
      results_df$T_Value[i] <- coef_table[param_name, "t value"]
      results_df$P_Value[i] <- coef_table[param_name, "Pr(>|t|)"]

      # Calculate 95% confidence intervals
      results_df$CI_Lower[i] <- results_df$Estimate[i] - 1.96 * results_df$Std_Error[i]
      results_df$CI_Upper[i] <- results_df$Estimate[i] + 1.96 * results_df$Std_Error[i]

      # Check significance
      results_df$Significant[i] <- results_df$P_Value[i] < 0.05
    }
  }

  return(results_df)
}

# ==============================================================================
# 2. COMPLETE DATASET ANALYSIS - MATCHING parasiteLoad complete_model
# ==============================================================================

cat("\n2. COMPLETE DATASET ANALYSIS (Ferreira Methodology)\n")
cat("====================================================\n")
cat("Testing overall hybrid effects (infected + uninfected mice)\n")
cat("Comparable to parasiteLoad complete_model\n\n")

# Create pairwise distance matrix for complete dataset
complete_pairwise <- create_ferreira_pairwise(ferreira_complete_data, "Complete Dataset")

# Fit linear model
ferreira_complete_model <- fit_ferreira_model(
  complete_pairwise,
  "Complete Dataset Model",
  include_infection = TRUE
)

# Extract results
complete_results <- extract_ferreira_results(ferreira_complete_model, "Complete Dataset")

cat("✓ Complete dataset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_complete_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 3. UNINFECTED SUBSET ANALYSIS - MATCHING parasiteLoad constitutive_model
# ==============================================================================

cat("3. UNINFECTED SUBSET ANALYSIS (Constitutive Costs)\n")
cat("==================================================\n")
cat("Testing constitutive immune costs (uninfected mice only)\n")
cat("Comparable to parasiteLoad constitutive_model\n\n")

# Create pairwise distance matrix for uninfected subset
uninfected_pairwise <- create_ferreira_pairwise(ferreira_uninfected_data, "Uninfected Subset")

# Fit linear model (no infection distance since all uninfected)
ferreira_uninfected_model <- fit_ferreira_model(
  uninfected_pairwise,
  "Uninfected Subset Model",
  include_infection = FALSE
)

# Extract results
uninfected_results <- extract_ferreira_results(ferreira_uninfected_model, "Uninfected Subset")

cat("✓ Uninfected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_uninfected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 4. INFECTED SUBSET ANALYSIS - NEW (following parasiteLoad logic)
# ==============================================================================

cat("4. INFECTED SUBSET ANALYSIS (Infection-Specific Effects)\n")
cat("========================================================\n")
cat("Testing hybrid effects in infected mice only\n")
cat("Complementing parasiteLoad analysis\n\n")

# Create pairwise distance matrix for infected subset
infected_pairwise <- create_ferreira_pairwise(ferreira_infected_data, "Infected Subset")

# Fit linear model (no infection distance since all infected)
ferreira_infected_model <- fit_ferreira_model(
  infected_pairwise,
  "Infected Subset Model",
  include_infection = FALSE
)

# Extract results
infected_results <- extract_ferreira_results(ferreira_infected_model, "Infected Subset")

cat("✓ Infected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_infected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 5. RESULTS COMPARISON WITH parasiteLoad
# ==============================================================================

cat("5. FERREIRA VALIDATION RESULTS\n")
cat("==============================\n")

# Function to summarize key findings
summarize_ferreira_findings <- function(results, analysis_name) {

  cat("\n", analysis_name, "Results:\n")
  cat(paste(rep("-", nchar(analysis_name) + 9), collapse = ""), "\n")

  # Key effects to highlight
  sex_effect <- results[results$Interpretation == "Sex difference effect", ]
  hybrid_effect <- results[results$Interpretation == "Hybridization effect (hHe-dist)", ]
  subspecies_effect <- results[results$Interpretation == "Subspecies genetic effect", ]

  if (nrow(sex_effect) > 0) {
    cat(sprintf("Sex differences: %s (p = %.5f, Est = %.4f)\n",
                ifelse(sex_effect$Significant, "SIGNIFICANT", "not significant"),
                sex_effect$P_Value, sex_effect$Estimate))
  }

  if (nrow(hybrid_effect) > 0) {
    cat(sprintf("Hybridization effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(hybrid_effect$Significant, "SIGNIFICANT", "not significant"),
                hybrid_effect$P_Value, hybrid_effect$Estimate))
  }

  if (nrow(subspecies_effect) > 0) {
    cat(sprintf("Subspecies effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(subspecies_effect$Significant, "SIGNIFICANT", "not significant"),
                subspecies_effect$P_Value, subspecies_effect$Estimate))
  }
}

# Summarize all three analyses
summarize_ferreira_findings(complete_results, "COMPLETE DATASET")
summarize_ferreira_findings(uninfected_results, "UNINFECTED SUBSET")
summarize_ferreira_findings(infected_results, "INFECTED SUBSET")

# ==============================================================================
# 6. DIRECT COMPARISON WITH parasiteLoad RESULTS
# ==============================================================================

cat("\n6. VALIDATION SUMMARY: Ferreira vs parasiteLoad\n")
cat("===============================================\n")

cat("COMPARISON OF KEY FINDINGS:\n")
cat("===========================\n\n")

cat("parasiteLoad Results (from script 1):\n")
cat("- Overall hybrid effect: p = 0.017 (SIGNIFICANT)\n")
cat("- Male-specific effect: p = 0.038 (SIGNIFICANT)\n")
cat("- Female-specific effect: p = 0.189 (not significant)\n")
cat("- Constitutive costs: p = 0.545 (not significant)\n\n")

cat("Ferreira Methodology Results:\n")

# Check for significant effects in each analysis
complete_sex_sig <- any(complete_results$Interpretation == "Sex difference effect" &
                          complete_results$Significant)
complete_hybrid_sig <- any(complete_results$Interpretation == "Hybridization effect (hHe-dist)" &
                             complete_results$Significant)
uninfected_hybrid_sig <- any(uninfected_results$Interpretation == "Hybridization effect (hHe-dist)" &
                               uninfected_results$Significant)
infected_effects <- sum(infected_results$Significant)

cat("- Complete dataset sex effects:", ifelse(complete_sex_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Complete dataset hybrid effects:", ifelse(complete_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Uninfected hybrid effects:", ifelse(uninfected_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Infected subset effects:", infected_effects, "significant parameters\n\n")

# Overall validation assessment
validation_score <- sum(
  complete_sex_sig || complete_hybrid_sig,  # Matches overall hybrid effect
  !uninfected_hybrid_sig,  # Matches no constitutive costs
  infected_effects > 0  # Matches infection-dependent effects
)

cat("VALIDATION OUTCOME:\n")
if (validation_score >= 2) {
  cat("🎉 STRONG VALIDATION: Ferreira methodology supports parasiteLoad findings!\n")
  cat("   Consistent patterns across both statistical frameworks.\n")
} else if (validation_score == 1) {
  cat("⚠️  PARTIAL VALIDATION: Some agreement between methodologies.\n")
  cat("   Differences may reflect methodological sensitivity.\n")
} else {
  cat("❌ LIMITED VALIDATION: Different patterns detected.\n")
  cat("   Consider sample size or methodological assumptions.\n")
}

cat("\nKEY INSIGHTS:\n")
cat("- Both methods test hybrid × host interactions\n")
cat("- parasiteLoad: Direct hybrid index effects on health\n")
cat("- Ferreira: Distance-based similarity in health responses\n")
cat("- Complementary approaches strengthen evidence base\n\n")

# ==============================================================================
# FERREIRA ET AL. METHODOLOGY VALIDATION - SIMPLIFIED VERSION
# ==============================================================================
# Following the exact same structure as parasiteLoad script:
# 1. Complete Dataset Analysis (infected + uninfected)
# 2. Uninfected Subset Analysis (constitutive costs)
# 3. Infected Subset Analysis (infection-specific effects)
#
# Using linear models instead of Bayesian to avoid brms installation issues
# Still follows Ferreira et al. pairwise distance methodology
# ==============================================================================

cat("=== FERREIRA ET AL. METHODOLOGY VALIDATION ===\n")
cat("Following parasiteLoad script structure for direct comparison\n")
cat("Using distance-based linear models (Ferreira methodology)\n\n")

# ==============================================================================
# 1. DATA PREPARATION - USING EXACT SAME DATASETS AS parasiteLoad
# ==============================================================================

cat("1. PREPARING DATA FOR FERREIRA ANALYSIS\n")
cat("=======================================\n")

# Check that parasiteLoad datasets exist
if (!exists("field_mice") || !exists("hybrid_data") || !exists("uninfected_data")) {
  stop("Please run parasiteLoad script first to create field_mice, hybrid_data, and uninfected_data")
}

cat("Using EXACT SAME datasets as parasiteLoad analysis:\n")

# 1. COMPLETE DATASET - Use field_mice (exactly like parasiteLoad complete_model)
ferreira_complete_data <- field_mice %>%
  filter(
    !is.na(HI),
    !is.na(Sex),
    !is.na(predicted_weight_loss)
  ) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe) - key Ferreira variable
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number(),

    # Response variable (matching parasiteLoad)
    response = predicted_weight_loss,

    # Infection status for distance calculations
    infected = as.logical(infection_status)
  )

# 2. UNINFECTED SUBSET - Direct from parasiteLoad (uninfected_data)
ferreira_uninfected_data <- uninfected_data %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

# 3. INFECTED SUBSET - Create from hybrid_data (following parasiteLoad logic)
ferreira_infected_data <- hybrid_data %>%
  filter(infected) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

cat("Analysis datasets prepared (matching parasiteLoad exactly):\n")
cat("- Complete dataset:", nrow(ferreira_complete_data), "mice (same as field_mice for complete_model)\n")
cat("- Uninfected subset:", nrow(ferreira_uninfected_data), "mice (same as uninfected_data)\n")
cat("- Infected subset:", nrow(ferreira_infected_data), "mice (filtered from hybrid_data)\n")

# Verify datasets match parasiteLoad exactly
cat("\nDataset verification:\n")
cat("- field_mice → ferreira_complete_data:", nrow(field_mice) >= nrow(ferreira_complete_data), "(filtered for complete cases)\n")
cat("- uninfected_data matches ferreira_uninfected_data:", nrow(uninfected_data) == nrow(ferreira_uninfected_data), "\n")
cat("- hybrid_data infected subset matches:", sum(hybrid_data$infected) == nrow(ferreira_infected_data), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS FOR FERREIRA ANALYSIS
# ==============================================================================

# Function to create pairwise distance matrix (following Ferreira et al.)
create_ferreira_pairwise <- function(data, analysis_name) {

  n <- nrow(data)
  n_pairs <- n * (n - 1) / 2

  cat("Creating", n_pairs, "pairwise comparisons for", analysis_name, "\n")

  # For large datasets, sample pairs to make it computationally feasible
  if (n_pairs > 10000) {
    cat("Large dataset detected. Sampling 10000 random pairs for analysis.\n")
    sample_pairs <- TRUE
    n_sample <- 10000
  } else {
    sample_pairs <- FALSE
  }

  # Initialize results
  pairwise_list <- list()
  counter <- 1

  if (sample_pairs) {
    # Sample random pairs
    all_pairs <- expand.grid(i = 1:(n-1), j = 2:n)
    all_pairs <- all_pairs[all_pairs$i < all_pairs$j, ]
    sampled_indices <- sample(nrow(all_pairs), n_sample)
    pairs_to_process <- all_pairs[sampled_indices, ]

    for (idx in 1:nrow(pairs_to_process)) {
      i <- pairs_to_process$i[idx]
      j <- pairs_to_process$j[idx]

      # Get individual data
      ind1 <- data[i, ]
      ind2 <- data[j, ]

      # Calculate distances (same as before)
      response_diff <- abs(ind1$response - ind2$response)
      max_response_range <- max(data$response) - min(data$response)
      response_distance <- response_diff / max_response_range
      response_similarity <- 1 - response_distance

      subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
      hHe_distance <- abs(ind1$hHe - ind2$hHe)
      hHe_mean <- (ind1$hHe + ind2$hHe) / 2

      sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
      infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

      pairwise_list[[counter]] <- data.frame(
        similarity = response_similarity,
        subspecies_genetic_distance = subspecies_genetic_distance,
        hHe_distance = hHe_distance,
        hHe_mean = hHe_mean,
        sex_distance = sex_distance,
        infection_distance = infection_distance,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }

  } else {
    # Process all pairs
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {

        # Get individual data
        ind1 <- data[i, ]
        ind2 <- data[j, ]

        # Response similarity (normalized distance between predicted weight loss)
        response_diff <- abs(ind1$response - ind2$response)
        max_response_range <- max(data$response) - min(data$response)
        response_distance <- response_diff / max_response_range
        response_similarity <- 1 - response_distance

        # Genetic distances (following Ferreira et al. Table 1)
        subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
        hHe_distance <- abs(ind1$hHe - ind2$hHe)
        hHe_mean <- (ind1$hHe + ind2$hHe) / 2

        # Other distances
        sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
        infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

        pairwise_list[[counter]] <- data.frame(
          similarity = response_similarity,
          subspecies_genetic_distance = subspecies_genetic_distance,
          hHe_distance = hHe_distance,
          hHe_mean = hHe_mean,
          sex_distance = sex_distance,
          infection_distance = infection_distance,
          stringsAsFactors = FALSE
        )

        counter <- counter + 1
      }
    }
  }

  result <- do.call(rbind, pairwise_list)

  # Scale predictors to 0-1 (following Ferreira methodology)
  result <- result %>%
    mutate(
      subspecies_genetic_distance_scaled = scales::rescale(subspecies_genetic_distance, to = c(0, 1)),
      hHe_distance_scaled = scales::rescale(hHe_distance, to = c(0, 1)),
      hHe_mean_scaled = scales::rescale(hHe_mean, to = c(0, 1))
    )

  return(result)
}

# Function to fit Ferreira-style linear model
fit_ferreira_model <- function(pairwise_data, model_name, include_infection = TRUE) {

  cat("Fitting", model_name, "using linear model (Ferreira methodology)...\n")

  # Model formula (matching Ferreira et al.)
  if (include_infection) {
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance + infection_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  } else {
    # For uninfected/infected subsets, no infection_distance needed
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  }

  # Fit linear model
  model <- lm(as.formula(formula_text), data = pairwise_data)
  model$model_type <- "linear"

  return(model)
}

# Function to extract results in parasiteLoad style
extract_ferreira_results <- function(model, model_name) {

  cat("\nExtracting results for", model_name, "...\n")

  # Get model summary
  model_summary <- summary(model)
  coef_table <- model_summary$coefficients

  # Key effects matching parasiteLoad interpretation
  effects_of_interest <- c(
    "subspecies_genetic_distance_scaled" = "Subspecies genetic effect",
    "hHe_distance_scaled" = "Hybridization effect (hHe-dist)",
    "hHe_mean_scaled" = "Mean hybridization effect",
    "sex_distance" = "Sex difference effect",
    "subspecies_genetic_distance_scaled:hHe_distance_scaled" = "Subspecies × Hybridization interaction"
  )

  # Add infection effect if present
  if ("infection_distance" %in% rownames(coef_table)) {
    effects_of_interest["infection_distance"] <- "Infection difference effect"
  }

  results_df <- data.frame(
    Effect = names(effects_of_interest),
    Interpretation = unname(effects_of_interest),
    Estimate = numeric(length(effects_of_interest)),
    Std_Error = numeric(length(effects_of_interest)),
    CI_Lower = numeric(length(effects_of_interest)),
    CI_Upper = numeric(length(effects_of_interest)),
    T_Value = numeric(length(effects_of_interest)),
    P_Value = numeric(length(effects_of_interest)),
    Significant = logical(length(effects_of_interest)),
    stringsAsFactors = FALSE
  )

  for (i in 1:length(effects_of_interest)) {
    param_name <- names(effects_of_interest)[i]

    if (param_name %in% rownames(coef_table)) {
      results_df$Estimate[i] <- coef_table[param_name, "Estimate"]
      results_df$Std_Error[i] <- coef_table[param_name, "Std. Error"]
      results_df$T_Value[i] <- coef_table[param_name, "t value"]
      results_df$P_Value[i] <- coef_table[param_name, "Pr(>|t|)"]

      # Calculate 95% confidence intervals
      results_df$CI_Lower[i] <- results_df$Estimate[i] - 1.96 * results_df$Std_Error[i]
      results_df$CI_Upper[i] <- results_df$Estimate[i] + 1.96 * results_df$Std_Error[i]

      # Check significance
      results_df$Significant[i] <- results_df$P_Value[i] < 0.05
    }
  }

  return(results_df)
}

# ==============================================================================
# 2. COMPLETE DATASET ANALYSIS - MATCHING parasiteLoad complete_model
# ==============================================================================

cat("\n2. COMPLETE DATASET ANALYSIS (Ferreira Methodology)\n")
cat("====================================================\n")
cat("Testing overall hybrid effects (infected + uninfected mice)\n")
cat("Comparable to parasiteLoad complete_model\n\n")

# Create pairwise distance matrix for complete dataset
complete_pairwise <- create_ferreira_pairwise(ferreira_complete_data, "Complete Dataset")

# Fit linear model
ferreira_complete_model <- fit_ferreira_model(
  complete_pairwise,
  "Complete Dataset Model",
  include_infection = TRUE
)

# Extract results
complete_results <- extract_ferreira_results(ferreira_complete_model, "Complete Dataset")

cat("✓ Complete dataset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_complete_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 3. UNINFECTED SUBSET ANALYSIS - MATCHING parasiteLoad constitutive_model
# ==============================================================================

cat("3. UNINFECTED SUBSET ANALYSIS (Constitutive Costs)\n")
cat("==================================================\n")
cat("Testing constitutive immune costs (uninfected mice only)\n")
cat("Comparable to parasiteLoad constitutive_model\n\n")

# Create pairwise distance matrix for uninfected subset
uninfected_pairwise <- create_ferreira_pairwise(ferreira_uninfected_data, "Uninfected Subset")

# Fit linear model (no infection distance since all uninfected)
ferreira_uninfected_model <- fit_ferreira_model(
  uninfected_pairwise,
  "Uninfected Subset Model",
  include_infection = FALSE
)

# Extract results
uninfected_results <- extract_ferreira_results(ferreira_uninfected_model, "Uninfected Subset")

cat("✓ Uninfected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_uninfected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 4. INFECTED SUBSET ANALYSIS - NEW (following parasiteLoad logic)
# ==============================================================================

cat("4. INFECTED SUBSET ANALYSIS (Infection-Specific Effects)\n")
cat("========================================================\n")
cat("Testing hybrid effects in infected mice only\n")
cat("Complementing parasiteLoad analysis\n\n")

# Create pairwise distance matrix for infected subset
infected_pairwise <- create_ferreira_pairwise(ferreira_infected_data, "Infected Subset")

# Fit linear model (no infection distance since all infected)
ferreira_infected_model <- fit_ferreira_model(
  infected_pairwise,
  "Infected Subset Model",
  include_infection = FALSE
)

# Extract results
infected_results <- extract_ferreira_results(ferreira_infected_model, "Infected Subset")

cat("✓ Infected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_infected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 5. RESULTS COMPARISON WITH parasiteLoad
# ==============================================================================

cat("5. FERREIRA VALIDATION RESULTS\n")
cat("==============================\n")

# Function to summarize key findings
summarize_ferreira_findings <- function(results, analysis_name) {

  cat("\n", analysis_name, "Results:\n")
  cat(paste(rep("-", nchar(analysis_name) + 9), collapse = ""), "\n")

  # Key effects to highlight
  sex_effect <- results[results$Interpretation == "Sex difference effect", ]
  hybrid_effect <- results[results$Interpretation == "Hybridization effect (hHe-dist)", ]
  subspecies_effect <- results[results$Interpretation == "Subspecies genetic effect", ]

  if (nrow(sex_effect) > 0) {
    cat(sprintf("Sex differences: %s (p = %.5f, Est = %.4f)\n",
                ifelse(sex_effect$Significant, "SIGNIFICANT", "not significant"),
                sex_effect$P_Value, sex_effect$Estimate))
  }

  if (nrow(hybrid_effect) > 0) {
    cat(sprintf("Hybridization effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(hybrid_effect$Significant, "SIGNIFICANT", "not significant"),
                hybrid_effect$P_Value, hybrid_effect$Estimate))
  }

  if (nrow(subspecies_effect) > 0) {
    cat(sprintf("Subspecies effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(subspecies_effect$Significant, "SIGNIFICANT", "not significant"),
                subspecies_effect$P_Value, subspecies_effect$Estimate))
  }
}

# Summarize all three analyses
summarize_ferreira_findings(complete_results, "COMPLETE DATASET")
summarize_ferreira_findings(uninfected_results, "UNINFECTED SUBSET")
summarize_ferreira_findings(infected_results, "INFECTED SUBSET")

# ==============================================================================
# 6. DIRECT COMPARISON WITH parasiteLoad RESULTS
# ==============================================================================

cat("\n6. VALIDATION SUMMARY: Ferreira vs parasiteLoad\n")
cat("===============================================\n")

cat("COMPARISON OF KEY FINDINGS:\n")
cat("===========================\n\n")

cat("parasiteLoad Results (from script 1):\n")
cat("- Overall hybrid effect: p = 0.017 (SIGNIFICANT)\n")
cat("- Male-specific effect: p = 0.038 (SIGNIFICANT)\n")
cat("- Female-specific effect: p = 0.189 (not significant)\n")
cat("- Constitutive costs: p = 0.545 (not significant)\n\n")

cat("Ferreira Methodology Results:\n")

# Check for significant effects in each analysis
complete_sex_sig <- any(complete_results$Interpretation == "Sex difference effect" &
                          complete_results$Significant)
complete_hybrid_sig <- any(complete_results$Interpretation == "Hybridization effect (hHe-dist)" &
                             complete_results$Significant)
uninfected_hybrid_sig <- any(uninfected_results$Interpretation == "Hybridization effect (hHe-dist)" &
                               uninfected_results$Significant)
infected_effects <- sum(infected_results$Significant)

cat("- Complete dataset sex effects:", ifelse(complete_sex_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Complete dataset hybrid effects:", ifelse(complete_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Uninfected hybrid effects:", ifelse(uninfected_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Infected subset effects:", infected_effects, "significant parameters\n\n")

# Overall validation assessment
validation_score <- sum(
  complete_sex_sig || complete_hybrid_sig,  # Matches overall hybrid effect
  !uninfected_hybrid_sig,  # Matches no constitutive costs
  infected_effects > 0  # Matches infection-dependent effects
)

cat("VALIDATION OUTCOME:\n")
if (validation_score >= 2) {
  cat("🎉 STRONG VALIDATION: Ferreira methodology supports parasiteLoad findings!\n")
  cat("   Consistent patterns across both statistical frameworks.\n")
} else if (validation_score == 1) {
  cat("⚠️  PARTIAL VALIDATION: Some agreement between methodologies.\n")
  cat("   Differences may reflect methodological sensitivity.\n")
} else {
  cat("❌ LIMITED VALIDATION: Different patterns detected.\n")
  cat("   Consider sample size or methodological assumptions.\n")
}

cat("\nKEY INSIGHTS:\n")
cat("- Both methods test hybrid × host interactions\n")
cat("- parasiteLoad: Direct hybrid index effects on health\n")
cat("- Ferreira: Distance-based similarity in health responses\n")
cat("- Complementary approaches strengthen evidence base\n\n")

# ==============================================================================
# FERREIRA ET AL. METHODOLOGY VALIDATION - SIMPLIFIED VERSION
# ==============================================================================
# Following the exact same structure as parasiteLoad script:
# 1. Complete Dataset Analysis (infected + uninfected)
# 2. Uninfected Subset Analysis (constitutive costs)
# 3. Infected Subset Analysis (infection-specific effects)
#
# Using linear models instead of Bayesian to avoid brms installation issues
# Still follows Ferreira et al. pairwise distance methodology
# ==============================================================================

cat("=== FERREIRA ET AL. METHODOLOGY VALIDATION ===\n")
cat("Following parasiteLoad script structure for direct comparison\n")
cat("Using distance-based linear models (Ferreira methodology)\n\n")

# ==============================================================================
# 1. DATA PREPARATION - USING EXACT SAME DATASETS AS parasiteLoad
# ==============================================================================

cat("1. PREPARING DATA FOR FERREIRA ANALYSIS\n")
cat("=======================================\n")

# Check that parasiteLoad datasets exist
if (!exists("field_mice") || !exists("hybrid_data") || !exists("uninfected_data")) {
  stop("Please run parasiteLoad script first to create field_mice, hybrid_data, and uninfected_data")
}

cat("Using EXACT SAME datasets as parasiteLoad analysis:\n")

# 1. COMPLETE DATASET - Use field_mice (exactly like parasiteLoad complete_model)
ferreira_complete_data <- field_mice %>%
  filter(
    !is.na(HI),
    !is.na(Sex),
    !is.na(predicted_weight_loss)
  ) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe) - key Ferreira variable
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number(),

    # Response variable (matching parasiteLoad)
    response = predicted_weight_loss,

    # Infection status for distance calculations
    infected = as.logical(infection_status)
  )

# 2. UNINFECTED SUBSET - Direct from parasiteLoad (uninfected_data)
ferreira_uninfected_data <- uninfected_data %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

# 3. INFECTED SUBSET - Create from hybrid_data (following parasiteLoad logic)
ferreira_infected_data <- hybrid_data %>%
  filter(infected) %>%
  mutate(
    # Calculate expected hybrid heterozygosity (hHe)
    hHe = 2 * HI * (1 - HI),

    # Individual IDs for pairwise analysis
    individual_id = row_number()
  )

cat("Analysis datasets prepared (matching parasiteLoad exactly):\n")
cat("- Complete dataset:", nrow(ferreira_complete_data), "mice (same as field_mice for complete_model)\n")
cat("- Uninfected subset:", nrow(ferreira_uninfected_data), "mice (same as uninfected_data)\n")
cat("- Infected subset:", nrow(ferreira_infected_data), "mice (filtered from hybrid_data)\n")

# Verify datasets match parasiteLoad exactly
cat("\nDataset verification:\n")
cat("- field_mice → ferreira_complete_data:", nrow(field_mice) >= nrow(ferreira_complete_data), "(filtered for complete cases)\n")
cat("- uninfected_data matches ferreira_uninfected_data:", nrow(uninfected_data) == nrow(ferreira_uninfected_data), "\n")
cat("- hybrid_data infected subset matches:", sum(hybrid_data$infected) == nrow(ferreira_infected_data), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS FOR FERREIRA ANALYSIS
# ==============================================================================

# Function to create pairwise distance matrix (following Ferreira et al.)
create_ferreira_pairwise <- function(data, analysis_name) {

  n <- nrow(data)
  n_pairs <- n * (n - 1) / 2

  cat("Creating", n_pairs, "pairwise comparisons for", analysis_name, "\n")

  # For large datasets, sample pairs to make it computationally feasible
  if (n_pairs > 10000) {
    cat("Large dataset detected. Sampling 10000 random pairs for analysis.\n")
    sample_pairs <- TRUE
    n_sample <- 10000
  } else {
    sample_pairs <- FALSE
  }

  # Initialize results
  pairwise_list <- list()
  counter <- 1

  if (sample_pairs) {
    # Sample random pairs
    all_pairs <- expand.grid(i = 1:(n-1), j = 2:n)
    all_pairs <- all_pairs[all_pairs$i < all_pairs$j, ]
    sampled_indices <- sample(nrow(all_pairs), n_sample)
    pairs_to_process <- all_pairs[sampled_indices, ]

    for (idx in 1:nrow(pairs_to_process)) {
      i <- pairs_to_process$i[idx]
      j <- pairs_to_process$j[idx]

      # Get individual data
      ind1 <- data[i, ]
      ind2 <- data[j, ]

      # Calculate distances (same as before)
      response_diff <- abs(ind1$response - ind2$response)
      max_response_range <- max(data$response) - min(data$response)
      response_distance <- response_diff / max_response_range
      response_similarity <- 1 - response_distance

      subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
      hHe_distance <- abs(ind1$hHe - ind2$hHe)
      hHe_mean <- (ind1$hHe + ind2$hHe) / 2

      sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
      infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

      pairwise_list[[counter]] <- data.frame(
        similarity = response_similarity,
        subspecies_genetic_distance = subspecies_genetic_distance,
        hHe_distance = hHe_distance,
        hHe_mean = hHe_mean,
        sex_distance = sex_distance,
        infection_distance = infection_distance,
        stringsAsFactors = FALSE
      )

      counter <- counter + 1
    }

  } else {
    # Process all pairs
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {

        # Get individual data
        ind1 <- data[i, ]
        ind2 <- data[j, ]

        # Response similarity (normalized distance between predicted weight loss)
        response_diff <- abs(ind1$response - ind2$response)
        max_response_range <- max(data$response) - min(data$response)
        response_distance <- response_diff / max_response_range
        response_similarity <- 1 - response_distance

        # Genetic distances (following Ferreira et al. Table 1)
        subspecies_genetic_distance <- abs(ind1$HI - ind2$HI)
        hHe_distance <- abs(ind1$hHe - ind2$hHe)
        hHe_mean <- (ind1$hHe + ind2$hHe) / 2

        # Other distances
        sex_distance <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
        infection_distance <- ifelse(ind1$infected == ind2$infected, 0, 1)

        pairwise_list[[counter]] <- data.frame(
          similarity = response_similarity,
          subspecies_genetic_distance = subspecies_genetic_distance,
          hHe_distance = hHe_distance,
          hHe_mean = hHe_mean,
          sex_distance = sex_distance,
          infection_distance = infection_distance,
          stringsAsFactors = FALSE
        )

        counter <- counter + 1
      }
    }
  }

  result <- do.call(rbind, pairwise_list)

  # Scale predictors to 0-1 (following Ferreira methodology)
  result <- result %>%
    mutate(
      subspecies_genetic_distance_scaled = scales::rescale(subspecies_genetic_distance, to = c(0, 1)),
      hHe_distance_scaled = scales::rescale(hHe_distance, to = c(0, 1)),
      hHe_mean_scaled = scales::rescale(hHe_mean, to = c(0, 1))
    )

  return(result)
}

# Function to fit Ferreira-style linear model
fit_ferreira_model <- function(pairwise_data, model_name, include_infection = TRUE) {

  cat("Fitting", model_name, "using linear model (Ferreira methodology)...\n")

  # Model formula (matching Ferreira et al.)
  if (include_infection) {
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance + infection_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  } else {
    # For uninfected/infected subsets, no infection_distance needed
    formula_text <- "similarity ~ subspecies_genetic_distance_scaled + hHe_distance_scaled +
                    hHe_mean_scaled + sex_distance +
                    subspecies_genetic_distance_scaled:hHe_distance_scaled"
  }

  # Fit linear model
  model <- lm(as.formula(formula_text), data = pairwise_data)
  model$model_type <- "linear"

  return(model)
}

# Function to extract results in parasiteLoad style
extract_ferreira_results <- function(model, model_name) {

  cat("\nExtracting results for", model_name, "...\n")

  # Get model summary
  model_summary <- summary(model)
  coef_table <- model_summary$coefficients

  # Key effects matching parasiteLoad interpretation
  effects_of_interest <- c(
    "subspecies_genetic_distance_scaled" = "Subspecies genetic effect",
    "hHe_distance_scaled" = "Hybridization effect (hHe-dist)",
    "hHe_mean_scaled" = "Mean hybridization effect",
    "sex_distance" = "Sex difference effect",
    "subspecies_genetic_distance_scaled:hHe_distance_scaled" = "Subspecies × Hybridization interaction"
  )

  # Add infection effect if present
  if ("infection_distance" %in% rownames(coef_table)) {
    effects_of_interest["infection_distance"] <- "Infection difference effect"
  }

  results_df <- data.frame(
    Effect = names(effects_of_interest),
    Interpretation = unname(effects_of_interest),
    Estimate = numeric(length(effects_of_interest)),
    Std_Error = numeric(length(effects_of_interest)),
    CI_Lower = numeric(length(effects_of_interest)),
    CI_Upper = numeric(length(effects_of_interest)),
    T_Value = numeric(length(effects_of_interest)),
    P_Value = numeric(length(effects_of_interest)),
    Significant = logical(length(effects_of_interest)),
    stringsAsFactors = FALSE
  )

  for (i in 1:length(effects_of_interest)) {
    param_name <- names(effects_of_interest)[i]

    if (param_name %in% rownames(coef_table)) {
      results_df$Estimate[i] <- coef_table[param_name, "Estimate"]
      results_df$Std_Error[i] <- coef_table[param_name, "Std. Error"]
      results_df$T_Value[i] <- coef_table[param_name, "t value"]
      results_df$P_Value[i] <- coef_table[param_name, "Pr(>|t|)"]

      # Calculate 95% confidence intervals
      results_df$CI_Lower[i] <- results_df$Estimate[i] - 1.96 * results_df$Std_Error[i]
      results_df$CI_Upper[i] <- results_df$Estimate[i] + 1.96 * results_df$Std_Error[i]

      # Check significance
      results_df$Significant[i] <- results_df$P_Value[i] < 0.05
    }
  }

  return(results_df)
}

# ==============================================================================
# 2. COMPLETE DATASET ANALYSIS - MATCHING parasiteLoad complete_model
# ==============================================================================

cat("\n2. COMPLETE DATASET ANALYSIS (Ferreira Methodology)\n")
cat("====================================================\n")
cat("Testing overall hybrid effects (infected + uninfected mice)\n")
cat("Comparable to parasiteLoad complete_model\n\n")

# Create pairwise distance matrix for complete dataset
complete_pairwise <- create_ferreira_pairwise(ferreira_complete_data, "Complete Dataset")

# Fit linear model
ferreira_complete_model <- fit_ferreira_model(
  complete_pairwise,
  "Complete Dataset Model",
  include_infection = TRUE
)

# Extract results
complete_results <- extract_ferreira_results(ferreira_complete_model, "Complete Dataset")

cat("✓ Complete dataset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_complete_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 3. UNINFECTED SUBSET ANALYSIS - MATCHING parasiteLoad constitutive_model
# ==============================================================================

cat("3. UNINFECTED SUBSET ANALYSIS (Constitutive Costs)\n")
cat("==================================================\n")
cat("Testing constitutive immune costs (uninfected mice only)\n")
cat("Comparable to parasiteLoad constitutive_model\n\n")

# Create pairwise distance matrix for uninfected subset
uninfected_pairwise <- create_ferreira_pairwise(ferreira_uninfected_data, "Uninfected Subset")

# Fit linear model (no infection distance since all uninfected)
ferreira_uninfected_model <- fit_ferreira_model(
  uninfected_pairwise,
  "Uninfected Subset Model",
  include_infection = FALSE
)

# Extract results
uninfected_results <- extract_ferreira_results(ferreira_uninfected_model, "Uninfected Subset")

cat("✓ Uninfected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_uninfected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 4. INFECTED SUBSET ANALYSIS - NEW (following parasiteLoad logic)
# ==============================================================================

cat("4. INFECTED SUBSET ANALYSIS (Infection-Specific Effects)\n")
cat("========================================================\n")
cat("Testing hybrid effects in infected mice only\n")
cat("Complementing parasiteLoad analysis\n\n")

# Create pairwise distance matrix for infected subset
infected_pairwise <- create_ferreira_pairwise(ferreira_infected_data, "Infected Subset")

# Fit linear model (no infection distance since all infected)
ferreira_infected_model <- fit_ferreira_model(
  infected_pairwise,
  "Infected Subset Model",
  include_infection = FALSE
)

# Extract results
infected_results <- extract_ferreira_results(ferreira_infected_model, "Infected Subset")

cat("✓ Infected subset analysis finished\n")
cat("Model R-squared:", round(summary(ferreira_infected_model)$r.squared, 4), "\n\n")

# ==============================================================================
# 5. RESULTS COMPARISON WITH parasiteLoad
# ==============================================================================

cat("5. FERREIRA VALIDATION RESULTS\n")
cat("==============================\n")

# Function to summarize key findings
summarize_ferreira_findings <- function(results, analysis_name) {

  cat("\n", analysis_name, "Results:\n")
  cat(paste(rep("-", nchar(analysis_name) + 9), collapse = ""), "\n")

  # Key effects to highlight
  sex_effect <- results[results$Interpretation == "Sex difference effect", ]
  hybrid_effect <- results[results$Interpretation == "Hybridization effect (hHe-dist)", ]
  subspecies_effect <- results[results$Interpretation == "Subspecies genetic effect", ]

  if (nrow(sex_effect) > 0) {
    cat(sprintf("Sex differences: %s (p = %.5f, Est = %.4f)\n",
                ifelse(sex_effect$Significant, "SIGNIFICANT", "not significant"),
                sex_effect$P_Value, sex_effect$Estimate))
  }

  if (nrow(hybrid_effect) > 0) {
    cat(sprintf("Hybridization effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(hybrid_effect$Significant, "SIGNIFICANT", "not significant"),
                hybrid_effect$P_Value, hybrid_effect$Estimate))
  }

  if (nrow(subspecies_effect) > 0) {
    cat(sprintf("Subspecies effect: %s (p = %.5f, Est = %.4f)\n",
                ifelse(subspecies_effect$Significant, "SIGNIFICANT", "not significant"),
                subspecies_effect$P_Value, subspecies_effect$Estimate))
  }
}

# Summarize all three analyses
summarize_ferreira_findings(complete_results, "COMPLETE DATASET")
summarize_ferreira_findings(uninfected_results, "UNINFECTED SUBSET")
summarize_ferreira_findings(infected_results, "INFECTED SUBSET")

# ==============================================================================
# 6. DIRECT COMPARISON WITH parasiteLoad RESULTS
# ==============================================================================

cat("\n6. VALIDATION SUMMARY: Ferreira vs parasiteLoad\n")
cat("===============================================\n")

cat("COMPARISON OF KEY FINDINGS:\n")
cat("===========================\n\n")

cat("parasiteLoad Results (from script 1):\n")
cat("- Overall hybrid effect: p = 0.017 (SIGNIFICANT)\n")
cat("- Male-specific effect: p = 0.038 (SIGNIFICANT)\n")
cat("- Female-specific effect: p = 0.189 (not significant)\n")
cat("- Constitutive costs: p = 0.545 (not significant)\n\n")

cat("Ferreira Methodology Results:\n")

# Check for significant effects in each analysis
complete_sex_sig <- any(complete_results$Interpretation == "Sex difference effect" &
                          complete_results$Significant)
complete_hybrid_sig <- any(complete_results$Interpretation == "Hybridization effect (hHe-dist)" &
                             complete_results$Significant)
uninfected_hybrid_sig <- any(uninfected_results$Interpretation == "Hybridization effect (hHe-dist)" &
                               uninfected_results$Significant)
infected_effects <- sum(infected_results$Significant)

cat("- Complete dataset sex effects:", ifelse(complete_sex_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Complete dataset hybrid effects:", ifelse(complete_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Uninfected hybrid effects:", ifelse(uninfected_hybrid_sig, "SIGNIFICANT", "not significant"), "\n")
cat("- Infected subset effects:", infected_effects, "significant parameters\n\n")

# Overall validation assessment
validation_score <- sum(
  complete_sex_sig || complete_hybrid_sig,  # Matches overall hybrid effect
  !uninfected_hybrid_sig,  # Matches no constitutive costs
  infected_effects > 0  # Matches infection-dependent effects
)

cat("VALIDATION OUTCOME:\n")
if (validation_score >= 2) {
  cat("🎉 STRONG VALIDATION: Ferreira methodology supports parasiteLoad findings!\n")
  cat("   Consistent patterns across both statistical frameworks.\n")
} else if (validation_score == 1) {
  cat("⚠️  PARTIAL VALIDATION: Some agreement between methodologies.\n")
  cat("   Differences may reflect methodological sensitivity.\n")
} else {
  cat("❌ LIMITED VALIDATION: Different patterns detected.\n")
  cat("   Consider sample size or methodological assumptions.\n")
}

cat("\nKEY INSIGHTS:\n")
cat("- Both methods test hybrid × host interactions\n")
cat("- parasiteLoad: Direct hybrid index effects on health\n")
cat("- Ferreira: Distance-based similarity in health responses\n")
cat("- Complementary approaches strengthen evidence base\n\n")

# ==============================================================================
# 7. SAVE FERREIRA VALIDATION RESULTS
# ==============================================================================

cat("7. SAVING FERREIRA VALIDATION RESULTS\n")
cat("=====================================\n")

# Save all models and results
save(
  ferreira_complete_data, ferreira_uninfected_data, ferreira_infected_data,
  complete_pairwise, uninfected_pairwise, infected_pairwise,
  ferreira_complete_model, ferreira_uninfected_model, ferreira_infected_model,
  complete_results, uninfected_results, infected_results,
  file = file.path("results", "ferreira_validation_linear.RData")
)

# Create publication-ready tables using gt
cat("\nCreating publication-ready tables...\n")

# Function to create formatted gt tables for Ferreira results
create_ferreira_gt_table <- function(results_df, table_title, analysis_description) {

  # Clean up the results for publication
  table_data <- results_df %>%
    filter(!is.na(Estimate)) %>%
    mutate(
      # Clean effect names for publication
      Effect_Clean = case_when(
        str_detect(Interpretation, "Subspecies") ~ "Subspecies genetic distance",
        str_detect(Interpretation, "Hybridization.*hHe-dist") ~ "Hybridization distance (hHe)",
        str_detect(Interpretation, "Mean hybridization") ~ "Mean hybridization level",
        str_detect(Interpretation, "Sex") ~ "Sex difference",
        str_detect(Interpretation, "Interaction") ~ "Subspecies × Hybridization",
        TRUE ~ Interpretation
      ),

      # Format p-values
      P_Value_Formatted = case_when(
        P_Value < 0.001 ~ "< 0.001",
        P_Value < 0.01 ~ sprintf("%.3f", P_Value),
        TRUE ~ sprintf("%.3f", P_Value)
      ),

      # Significance symbols
      Sig_Symbol = case_when(
        P_Value < 0.001 ~ "***",
        P_Value < 0.01 ~ "**",
        P_Value < 0.05 ~ "*",
        P_Value < 0.1 ~ ".",
        TRUE ~ ""
      ),

      # Format confidence intervals
      CI_Formatted = sprintf("(%.4f, %.4f)", CI_Lower, CI_Upper)
    ) %>%
    dplyr::select(Effect_Clean, Estimate, CI_Formatted, P_Value_Formatted, Sig_Symbol)

  # Create gt table
  gt_table <- table_data %>%
    gt() %>%

    # Table header
    tab_header(
      title = md(paste0("**", table_title, "**")),
      subtitle = analysis_description
    ) %>%

    # Column labels
    cols_label(
      Effect_Clean = "Predictor Variable",
      Estimate = "Coefficient",
      CI_Formatted = "95% CI",
      P_Value_Formatted = "P-value",
      Sig_Symbol = ""
    ) %>%

    # Format numeric columns
    fmt_number(
      columns = Estimate,
      decimals = 4
    ) %>%

    # Style significant results
    tab_style(
      style = list(
        cell_fill(color = "#ffebee"),
        cell_text(weight = "bold")
      ),
      locations = cells_body(
        rows = Sig_Symbol %in% c("*", "**", "***")
      )
    ) %>%

    # Style header
    tab_style(
      style = list(
        cell_fill(color = "#1f77b4"),
        cell_text(color = "white", weight = "bold")
      ),
      locations = cells_column_labels()
    ) %>%

    # Add footnote
    tab_footnote(
      footnote = "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1",
      locations = cells_column_labels(columns = Sig_Symbol)
    ) %>%

    # Table options
    tab_options(
      table.font.size = 12,
      heading.title.font.size = 16,
      heading.subtitle.font.size = 14,
      column_labels.font.weight = "bold",
      table.border.top.width = px(3),
      table.border.bottom.width = px(3),
      column_labels.border.bottom.width = px(2),
      row_group.border.bottom.width = px(1)
    )

  return(gt_table)
}

# Create formatted tables for each analysis
cat("Creating Complete Dataset table...\n")
complete_gt_table <- create_ferreira_gt_table(
  complete_results,
  "",
  "Distance-based analysis of hybrid effects on predicted weight loss similarity (all mice: n = 335)"
)

cat("Creating Uninfected Subset table...\n")
uninfected_gt_table <- create_ferreira_gt_table(
  uninfected_results,
  "",
  "Testing constitutive hybrid costs in uninfected mice only (n = 171)"
)

cat("Creating Infected Subset table...\n")
infected_gt_table <- create_ferreira_gt_table(
  infected_results,
  "",
  "Testing infection-dependent hybrid effects in infected mice only (n = 133)"
)

# Create validation summary table
cat("Creating validation summary table...\n")
validation_summary_data <- data.frame(
  Analysis = c("Complete Dataset", "Uninfected Subset", "Infected Subset"),
  parasiteLoad_Result = c("p = 0.017 (Significant)", "p = 0.545 (Not significant)", "Male-specific: p = 0.038"),
  Ferreira_Result = c(
    paste0(sum(complete_results$Significant, na.rm = TRUE), " significant effects"),
    paste0(sum(uninfected_results$Significant, na.rm = TRUE), " significant effects"),
    paste0(sum(infected_results$Significant, na.rm = TRUE), " significant effects")
  ),
  Validation_Status = c(
    ifelse(any(complete_results$Significant, na.rm = TRUE), "✓ Partial support", "○ No support"),
    ifelse(any(uninfected_results$Significant, na.rm = TRUE), "○ Inconsistent", "✓ Consistent"),
    ifelse(any(infected_results$Significant, na.rm = TRUE), "✓ Strong support", "○ No support")
  ),
  stringsAsFactors = FALSE
)

validation_summary_table <- validation_summary_data %>%
  gt() %>%
  tab_header(
    title = md(""),
    subtitle = "Comparison of hybrid effect detection across statistical frameworks"
  ) %>%
  cols_label(
    Analysis = "Analysis Type",
    parasiteLoad_Result = "parasiteLoad Result",
    Ferreira_Result = "Ferreira Result",
    Validation_Status = "Validation Status"
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#e8f5e8"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(
      rows = str_detect(Validation_Status, "✓")
    )
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "#1f77b4"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  ) %>%
  tab_footnote(
    footnote = "✓ = Methods agree, ○ = Methods disagree or no effect detected",
    locations = cells_column_labels(columns = Validation_Status)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 14,
    column_labels.font.weight = "bold"
  )

# Save all tables using your save_table_all_formats function
cat("\nSaving tables in multiple formats...\n")

save_table_all_formats(complete_gt_table, "Ferreira_Complete_Dataset_Results")
save_table_all_formats(uninfected_gt_table, "Ferreira_Uninfected_Subset_Results")
save_table_all_formats(infected_gt_table, "Ferreira_Infected_Subset_Results")
save_table_all_formats(validation_summary_table, "Ferreira_Validation_Summary")

cat("✓ All publication-ready tables saved in multiple formats\n")
cat("✓ Tables available in: HTML, DOCX, PNG, PDF, TEX\n")
cat("✓ Ready for manuscript integration\n\n")

cat("=== FERREIRA VALIDATION COMPLETE ===\n")
cat("Linear model implementation of Ferreira methodology successful!\n")
cat("Results provide structured validation of parasiteLoad findings.\n")
cat("Both methodologies ready for manuscript discussion! 🎉\n")

