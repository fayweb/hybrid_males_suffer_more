# FERREIRA METHODOLOGY VALIDATION WITH SEX-STRATIFIED ANALYSIS
# Cross-validate parasiteLoad findings using distance-based framework
# Testing sex-specific hybridization effects

# Load required packages
library(dplyr)
library(scales)
library(gt)

# VALIDATION: Ensure prerequisite datasets exist
required_datasets <- c("field_mice", "hybrid_data", "uninfected_data")
missing_datasets <- required_datasets[!sapply(required_datasets, exists)]
if (length(missing_datasets) > 0) {
  stop("Missing required datasets: ", paste(missing_datasets, collapse = ", "))
}

# DATA PREPARATION
prepare_ferreira_data <- function(data, analysis_name) {
  data %>%
    filter(!is.na(HI), !is.na(Sex), !is.na(predicted_weight_loss)) %>%
    mutate(
      hHe = 2 * HI * (1 - HI),
      response = predicted_weight_loss,
      infected = as.logical(infection_status)
    )
}

# Prepare datasets
complete_data <- prepare_ferreira_data(field_mice, "Complete")
uninfected_data_clean <- prepare_ferreira_data(uninfected_data, "Uninfected")
infected_data_clean <- hybrid_data %>% filter(infected) %>% prepare_ferreira_data("Infected")

# PAIRWISE DISTANCE CALCULATION
create_pairwise_distances <- function(data, analysis_name) {
  n <- nrow(data)
  n_pairs <- n * (n - 1) / 2

  cat("Creating", n_pairs, "pairwise comparisons for", analysis_name, "\n")

  # Sample large datasets for computational efficiency
  if (n_pairs > 10000) {
    cat("Large dataset: sampling 10,000 pairs\n")
    all_pairs <- expand.grid(i = 1:(n-1), j = 2:n)
    all_pairs <- all_pairs[all_pairs$i < all_pairs$j, ]
    pairs_to_process <- all_pairs[sample(nrow(all_pairs), 10000), ]
  } else {
    pairs_to_process <- expand.grid(i = 1:(n-1), j = 2:n)
    pairs_to_process <- pairs_to_process[pairs_to_process$i < pairs_to_process$j, ]
  }

  # Calculate pairwise distances
  result_list <- list()
  for (idx in 1:nrow(pairs_to_process)) {
    i <- pairs_to_process$i[idx]
    j <- pairs_to_process$j[idx]

    ind1 <- data[i, ]
    ind2 <- data[j, ]

    # Response similarity
    response_diff <- abs(ind1$response - ind2$response)
    max_range <- max(data$response) - min(data$response)
    similarity <- 1 - (response_diff / max_range)

    # Genetic distances
    subspecies_dist <- abs(ind1$HI - ind2$HI)
    hHe_dist <- abs(ind1$hHe - ind2$hHe)
    hHe_mean <- (ind1$hHe + ind2$hHe) / 2

    # Categorical distances
    sex_dist <- ifelse(ind1$Sex == ind2$Sex, 0, 1)
    infection_dist <- ifelse(ind1$infected == ind2$infected, 0, 1)

    result_list[[idx]] <- data.frame(
      similarity = similarity,
      subspecies_distance = subspecies_dist,
      hHe_distance = hHe_dist,
      hHe_mean = hHe_mean,
      sex_distance = sex_dist,
      infection_distance = infection_dist
    )
  }

  result <- do.call(rbind, result_list)

  # Scale predictors to 0-1
  result %>%
    mutate(
      subspecies_distance = rescale(subspecies_distance, to = c(0, 1)),
      hHe_distance = rescale(hHe_distance, to = c(0, 1)),
      hHe_mean = rescale(hHe_mean, to = c(0, 1))
    )
}

# FERREIRA MODEL FITTING
fit_ferreira_model <- function(pairwise_data, include_infection = TRUE) {
  if (include_infection) {
    formula <- "similarity ~ subspecies_distance + hHe_distance + hHe_mean +
                sex_distance + infection_distance + subspecies_distance:hHe_distance"
  } else {
    formula <- "similarity ~ subspecies_distance + hHe_distance + hHe_mean +
                sex_distance + subspecies_distance:hHe_distance"
  }

  lm(as.formula(formula), data = pairwise_data)
}

# RESULTS EXTRACTION
extract_results <- function(model, model_name) {
  coef_table <- summary(model)$coefficients

  # Define effects of interest
  effects <- c(
    "subspecies_distance" = "Subspecies distance",
    "hHe_distance" = "Hybridization distance",
    "hHe_mean" = "Mean hybridization",
    "sex_distance" = "Sex difference",
    "infection_distance" = "Infection difference",
    "subspecies_distance:hHe_distance" = "Subspecies × Hybridization"
  )

  # Extract results for each effect
  results <- data.frame(
    Effect = names(effects),
    Interpretation = unname(effects),
    Estimate = NA,
    Std_Error = NA,
    P_Value = NA,
    Significant = FALSE
  )

  for (i in 1:length(effects)) {
    param <- names(effects)[i]
    if (param %in% rownames(coef_table)) {
      results$Estimate[i] <- coef_table[param, "Estimate"]
      results$Std_Error[i] <- coef_table[param, "Std. Error"]
      results$P_Value[i] <- coef_table[param, "Pr(>|t|)"]
      results$Significant[i] <- results$P_Value[i] < 0.05
    }
  }

  results$Analysis <- model_name
  results$R_squared <- summary(model)$r.squared

  return(results)
}

# MAIN ANALYSIS: Run all analyses including sex-stratified

cat("Running Ferreira validation analyses...\n\n")

# 1. Complete Dataset Analysis
cat("1. Complete Dataset Analysis\n")
complete_pairwise <- create_pairwise_distances(complete_data, "Complete Dataset")
complete_model <- fit_ferreira_model(complete_pairwise, include_infection = TRUE)
complete_results <- extract_results(complete_model, "Complete Dataset")

# 2. Uninfected Subset Analysis
cat("\n2. Uninfected Subset Analysis\n")
uninfected_pairwise <- create_pairwise_distances(uninfected_data_clean, "Uninfected Subset")
uninfected_model <- fit_ferreira_model(uninfected_pairwise, include_infection = FALSE)
uninfected_results <- extract_results(uninfected_model, "Uninfected Subset")

# 3. Infected Subset Analysis
cat("\n3. Infected Subset Analysis\n")
infected_pairwise <- create_pairwise_distances(infected_data_clean, "Infected Subset")
infected_model <- fit_ferreira_model(infected_pairwise, include_infection = FALSE)
infected_results <- extract_results(infected_model, "Infected Subset")

# 4. SEX-STRATIFIED ANALYSES
cat("\n4. Sex-Stratified Analyses\n")

# Males only - Complete dataset
cat("4a. Males Only - Complete Dataset\n")
male_complete_data <- complete_data %>% filter(Sex == "M")
male_complete_pairwise <- create_pairwise_distances(male_complete_data, "Males Complete")
male_complete_model <- fit_ferreira_model(male_complete_pairwise, include_infection = TRUE)
male_complete_results <- extract_results(male_complete_model, "Males Complete")
male_complete_results

# Females only - Complete dataset
cat("4b. Females Only - Complete Dataset\n")
female_complete_data <- complete_data %>% filter(Sex == "F")
female_complete_pairwise <- create_pairwise_distances(female_complete_data, "Females Complete")
female_complete_model <- fit_ferreira_model(female_complete_pairwise, include_infection = TRUE)
female_complete_results <- extract_results(female_complete_model, "Females Complete")
female_complete_results

# Males only - Infected subset
cat("4c. Males Only - Infected Subset\n")
male_infected_data <- infected_data_clean %>% filter(Sex == "Male")
male_infected_pairwise <- create_pairwise_distances(male_infected_data, "Males Infected")
male_infected_model <- fit_ferreira_model(male_infected_pairwise, include_infection = FALSE)
male_infected_results <- extract_results(male_infected_model, "Males Infected")
male_infected_results


# Females only - Infected subset
cat("4d. Females Only - Infected Subset\n")
female_infected_data <- infected_data_clean %>% filter(Sex == "Female")
female_infected_pairwise <- create_pairwise_distances(female_infected_data, "Females Infected")
female_infected_model <- fit_ferreira_model(female_infected_pairwise, include_infection = FALSE)
female_infected_results <- extract_results(female_infected_model, "Females Infected")
female_infected_results

# COMBINE ALL RESULTS
all_results <- bind_rows(
  complete_results,
  uninfected_results,
  infected_results,
  male_complete_results,
  female_complete_results,
  male_infected_results,
  female_infected_results
)

# SUMMARY FUNCTION
summarize_findings <- function(results_df, analysis_name) {
  cat("\n", analysis_name, " (R² = ", sprintf("%.4f", unique(results_df$R_squared)), ")\n", sep = "")
  cat(rep("-", nchar(analysis_name) + 15), "\n", sep = "")

  sig_effects <- results_df %>% filter(Significant == TRUE)
  if (nrow(sig_effects) > 0) {
    for (i in 1:nrow(sig_effects)) {
      cat(sprintf("• %s: β = %.4f, p = %.4f\n",
                  sig_effects$Interpretation[i],
                  sig_effects$Estimate[i],
                  sig_effects$P_Value[i]))
    }
  } else {
    cat("• No significant effects detected\n")
  }
}

# PRINT RESULTS SUMMARY
cat("\n\nRESULTS SUMMARY\n")
cat("===============\n")

summarize_findings(complete_results, "Complete Dataset")
summarize_findings(uninfected_results, "Uninfected Subset")
summarize_findings(infected_results, "Infected Subset")
summarize_findings(male_complete_results, "Males Complete")
summarize_findings(female_complete_results, "Females Complete")
summarize_findings(male_infected_results, "Males Infected")
summarize_findings(female_infected_results, "Females Infected")

# SEX-SPECIFIC COMPARISON
cat("\n\nSEX-SPECIFIC HYBRID EFFECTS COMPARISON\n")
cat("======================================\n")

# Compare hybridization effects between sexes
male_hybrid_effect <- male_infected_results %>% filter(Effect == "hHe_distance")
female_hybrid_effect <- female_infected_results %>% filter(Effect == "hHe_distance")

if (nrow(male_hybrid_effect) > 0 && nrow(female_hybrid_effect) > 0) {
  cat("Hybridization Distance Effects in Infected Mice:\n")
  cat(sprintf("• Males: β = %.4f, p = %.4f (%s)\n",
              male_hybrid_effect$Estimate,
              male_hybrid_effect$P_Value,
              ifelse(male_hybrid_effect$Significant, "SIGNIFICANT", "not significant")))
  cat(sprintf("• Females: β = %.4f, p = %.4f (%s)\n",
              female_hybrid_effect$Estimate,
              female_hybrid_effect$P_Value,
              ifelse(female_hybrid_effect$Significant, "SIGNIFICANT", "not significant")))
}

# Save results
save(all_results, complete_model, uninfected_model, infected_model,
     male_complete_model, female_complete_model,
     male_infected_model, female_infected_model,
     file = "ferreira_sex_stratified_results.RData")

cat("\n✓ Results saved to: ferreira_sex_stratified_results.RData\n")
cat("✓ Analysis complete!\n")
