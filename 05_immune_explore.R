# ==============================================================================
# IMMUNE ANALYSIS WITH HYBRID INDEX × HETEROZYGOSITY INTERACTIONS
# ==============================================================================
# Purpose: Test sex-specific hybrid effects on immune responses with infection
# Focus: Complete dataset analysis with infection as interaction term
# Author: Fay Webster
# Date: December 2024
# ==============================================================================

# Check if running from master script
if (!exists("field_mice") || !exists("immune_genes")) {
  stop("Please run 00_master_script.R first to load data and packages")
}

cat("\n=== IMMUNE ANALYSIS: HI × He INTERACTIONS WITH INFECTION ===\n")
cat("Complete dataset approach for maximum statistical power\n\n")

# Load required packages
library(tidyverse)
library(vegan)
library(ggeffects)
library(gt)
library(car)
library(umap)

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("1. PREPARING COMPLETE DATASET\n")
cat("=============================\n")

# Filter to available genes
available_genes <- immune_genes[immune_genes %in% names(field_mice)]
cat("Available immune genes:", length(available_genes), "\n")

# Create complete dataset with all interactions
immune_data <- field_mice %>%
  dplyr::select(Mouse_ID, HI, Sex, infection_status, predicted_weight_loss,
                all_of(available_genes)) %>%
  filter(!is.na(HI) & !is.na(Sex) & !is.na(infection_status) &
           !is.na(predicted_weight_loss)) %>%
  mutate(
    # Core variables
    Sex = factor(Sex, levels = c("F", "M")),
    He = 2 * HI * (1 - HI),  # Heterozygosity - KEY for non-linear effects

    # Infection status
    infection = factor(infection_status,
                       levels = c(FALSE, TRUE),
                       labels = c("Uninfected", "Infected")),

    # Hybrid categories for visualization
    hybrid_category = case_when(
      HI < 0.2 ~ "Parental",
      HI > 0.8 ~ "Parental",
      TRUE ~ "Hybrid"
    )
  )

# Handle missing values in immune genes
immune_data <- immune_data %>%
  group_by(Sex, infection) %>%
  mutate(across(all_of(available_genes),
                ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  ungroup()

cat("✓ Total mice:", nrow(immune_data), "\n")
cat("✓ By infection status:\n")
print(table(immune_data$infection, immune_data$Sex))
cat("\n")

# ==============================================================================
# 2. CREATE BIOLOGICAL PATHWAY SCORES
# ==============================================================================

cat("2. CREATING IMMUNE PATHWAY SCORES\n")
cat("=================================\n")

# Calculate pathway scores using complete dataset
calculate_pathways <- function(data, genes) {
  pathways <- list()

  # Th1 response (anti-intracellular pathogens)
  th1_genes <- intersect(c("IFNy", "TNF", "CXCL9", "IDO1"), genes)
  if (length(th1_genes) >= 1) {
    pathways$Th1_response <- rowMeans(data[, th1_genes, drop = FALSE], na.rm = TRUE)
  }

  # Th2/mucus response (anti-Eimeria)
  th2_genes <- intersect(c("IL.13", "IL.10", "MUC2", "MUC5AC"), genes)
  if (length(th2_genes) >= 1) {
    pathways$Th2_mucus <- rowMeans(data[, th2_genes, drop = FALSE], na.rm = TRUE)
  }

  # Inflammation/tissue damage
  inflam_genes <- intersect(c("IL.6", "TNF", "IL.17", "CASP1", "MPO"), genes)
  if (length(inflam_genes) >= 1) {
    pathways$Inflammation <- rowMeans(data[, inflam_genes, drop = FALSE], na.rm = TRUE)
  }

  # Regulatory/tolerance
  reg_genes <- intersect(c("IL1RN", "IL.10", "SOCS1", "IDO1"), genes)
  if (length(reg_genes) >= 1) {
    pathways$Regulation <- rowMeans(data[, reg_genes, drop = FALSE], na.rm = TRUE)
  }

  # Cytotoxicity
  cyto_genes <- intersect(c("PRF1", "NCR1"), genes)
  if (length(cyto_genes) >= 1) {
    pathways$Cytotoxicity <- rowMeans(data[, cyto_genes, drop = FALSE], na.rm = TRUE)
  }

  # Innate immunity
  innate_genes <- intersect(c("MYD88", "TICAM1", "IRGM1"), genes)
  if (length(innate_genes) >= 1) {
    pathways$Innate_immunity <- rowMeans(data[, innate_genes, drop = FALSE], na.rm = TRUE)
  }

  # Th1/Th2 balance (critical for parasite response)
  if ("Th1_response" %in% names(pathways) && "Th2_mucus" %in% names(pathways)) {
    pathways$Th1_Th2_balance <- pathways$Th1_response - pathways$Th2_mucus
  }

  return(pathways)
}

# Add pathway scores to dataset
pathway_scores <- calculate_pathways(immune_data, available_genes)
for (pathway in names(pathway_scores)) {
  immune_data[[pathway]] <- pathway_scores[[pathway]]
}
pathway_names <- names(pathway_scores)

cat("✓ Created", length(pathway_names), "pathway scores\n\n")

# ==============================================================================
# 3. TEST MULTIPLE MODEL STRUCTURES FOR PREDICTED WEIGHT LOSS
# ==============================================================================

cat("3. TESTING MODEL STRUCTURES FOR PREDICTED WEIGHT LOSS\n")
cat("====================================================\n")

# Model 1: Simple HI × He
model1 <- lm(predicted_weight_loss ~ HI * He, data = immune_data)

# Model 2: Main effects only
model2 <- lm(predicted_weight_loss ~ Sex + HI + He + infection, data = immune_data)

# Model 3: Two-way interactions
model3 <- lm(predicted_weight_loss ~ Sex * HI + Sex * He + HI * He +
               Sex * infection + HI * infection + He * infection,
             data = immune_data)

# Model 4: Three-way interactions (full model)
model4 <- lm(predicted_weight_loss ~ Sex * HI * He * infection, data = immune_data)

# Model 5: Focused model (biologically motivated)
model5 <- lm(predicted_weight_loss ~ Sex * (HI + He) * infection + HI:He,
             data = immune_data)

# Compare models
cat("Model comparison (AIC):\n")
model_comparison <- data.frame(
  Model = c("HI × He only", "Main effects", "Two-way", "Full three-way", "Focused"),
  AIC = c(AIC(model1), AIC(model2), AIC(model3), AIC(model4), AIC(model5)),
  df = c(model1$df.residual, model2$df.residual, model3$df.residual,
         model4$df.residual, model5$df.residual)
)
print(model_comparison %>% arrange(AIC))

# Select best model
best_model <- model5  # or choose based on AIC
cat("\nUsing focused model for further analysis\n")

# ==============================================================================
# 4. VISUALIZE MODEL PREDICTIONS
# ==============================================================================

cat("\n4. CREATING MODEL VISUALIZATIONS\n")
cat("================================\n")

# A. Main effects plot using ggpredict
pred_main <- ggpredict(best_model, terms = c("HI [all]", "infection", "Sex"))

plot_main_effects <- ggplot(pred_main, aes(x = x, y = predicted, color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = c("Uninfected" = "#4daf4a", "Infected" = "#ff7f00"),
                     name = "Infection Status") +
  scale_fill_manual(values = c("Uninfected" = "#4daf4a", "Infected" = "#ff7f00"),
                    name = "Infection Status") +
  labs(title = "Predicted Weight Loss by Hybrid Index, Sex, and Infection",
       x = "Hybrid Index",
       y = "Predicted Weight Loss (%)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

save_plot_all_formats(plot_main_effects, "weight_loss_main_effects")

# B. Heterozygosity effect visualization
pred_he <- ggpredict(best_model, terms = c("He [all]", "infection", "Sex"))

plot_he_effects <- ggplot(pred_he, aes(x = x, y = predicted, color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094")) +
  scale_fill_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094")) +
  labs(title = "Heterozygosity Effects on Predicted Weight Loss",
       x = "Heterozygosity (He)",
       y = "Predicted Weight Loss (%)") +
  theme_minimal()

save_plot_all_formats(plot_he_effects, "weight_loss_heterozygosity_effects")

# C. Coefficient plot
library(broom)
library(ggplot2)

coef_data <- tidy(best_model, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = str_replace_all(term, ":", " × "),
    significant = p.value < 0.05
  )

plot_coefficients <- ggplot(coef_data, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     name = "Significant\n(p < 0.05)") +
  labs(title = "Model Coefficients for Predicted Weight Loss",
       x = "Coefficient Estimate",
       y = "") +
  theme_minimal()

plot_coefficients
save_plot_all_formats(plot_coefficients, "weight_loss_coefficients")

# ==============================================================================
# 5. PATHWAY ANALYSIS WITH INTERACTIONS
# ==============================================================================

cat("\n5. TESTING PATHWAYS WITH HI × He × SEX × INFECTION\n")
cat("==================================================\n")

# Function to test each pathway
test_pathway_interactions <- function(pathway_name, data) {
  formula <- as.formula(paste(pathway_name, "~ Sex * (HI + He) * infection + HI:He"))
  model <- lm(formula, data = data)

  # Extract key p-values
  anova_result <- anova(model)

  # Get specific interaction p-values
  coef_table <- summary(model)$coefficients

  results <- data.frame(
    Pathway = pathway_name,
    Sex_HI_p = ifelse("SexM:HI" %in% rownames(coef_table),
                      coef_table["SexM:HI", "Pr(>|t|)"], NA),
    Sex_He_p = ifelse("SexM:He" %in% rownames(coef_table),
                      coef_table["SexM:He", "Pr(>|t|)"], NA),
    HI_He_p = ifelse("HI:He" %in% rownames(coef_table),
                     coef_table["HI:He", "Pr(>|t|)"], NA),
    Sex_infection_p = ifelse("SexM:infectionInfected" %in% rownames(coef_table),
                             coef_table["SexM:infectionInfected", "Pr(>|t|)"], NA),
    R_squared = summary(model)$r.squared
  )

  return(results)
}

# Test all pathways
pathway_results <- map_df(pathway_names, ~test_pathway_interactions(.x, immune_data))

# Apply FDR correction
pathway_results <- pathway_results %>%
  mutate(
    Sex_HI_padj = p.adjust(Sex_HI_p, method = "BH"),
    Sex_He_padj = p.adjust(Sex_He_p, method = "BH"),
    HI_He_padj = p.adjust(HI_He_p, method = "BH"),
    Sex_infection_padj = p.adjust(Sex_infection_p, method = "BH")
  ) %>%
  arrange(Sex_He_p)

cat("Pathway interaction results (FDR corrected):\n")
print(pathway_results)

# Visualize top pathway
if (nrow(pathway_results) > 0) {
  top_pathway <- pathway_results$Pathway[1]
  formula_top <- as.formula(paste(top_pathway, "~ Sex * (HI + He) * infection + HI:He"))
  model_top <- lm(formula_top, data = immune_data)

  pred_pathway <- ggpredict(model_top, terms = c("HI [all]", "infection", "Sex"))

  plot_top_pathway <- ggplot(pred_pathway, aes(x = x, y = predicted, color = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.2, color = NA) +
    geom_line(size = 1.2) +
    facet_wrap(~ facet, nrow = 1) +
    scale_color_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094")) +
    scale_fill_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094")) +
    labs(title = paste(top_pathway, "Response by HI, Sex, and Infection"),
         x = "Hybrid Index",
         y = paste(top_pathway, "Score")) +
    theme_minimal()

  save_plot_all_formats(plot_top_pathway, paste0("pathway_", top_pathway))
}

plot_top_pathway
# ==============================================================================
# 6. PERMANOVA WITH INFECTION STATUS
# ==============================================================================

cat("\n6. PERMANOVA - MULTIVARIATE IMMUNE DIFFERENCES\n")
cat("==============================================\n")

# Prepare immune matrix
immune_matrix <- as.matrix(immune_data[, available_genes])

# Full model with infection
permanova_full <- adonis2(
  immune_matrix ~ Sex * HI * He * infection,
  data = immune_data,
  permutations = 999,
  method = "euclidean"
)

cat("PERMANOVA Results (with infection):\n")
print(permanova_full)

# Simplified model for comparison
permanova_simple <- adonis2(
  immune_matrix ~ Sex + HI + He + infection +
    Sex:HI + Sex:He + HI:He + Sex:infection + HI:infection,
  data = immune_data,
  permutations = 999,
  method = "euclidean"
)

cat("\nSimplified PERMANOVA:\n")
print(permanova_simple)

# ==============================================================================
# 7. NETWORK COHERENCE WITH INFECTION AS FACTOR
# ==============================================================================

cat("\n7. NETWORK COHERENCE ANALYSIS\n")
cat("=============================\n")

# Function to calculate network coherence
calculate_coherence <- function(data, genes) {
  cor_matrix <- cor(data[, genes], use = "pairwise.complete.obs")
  mean(abs(cor_matrix[upper.tri(cor_matrix)]), na.rm = TRUE)
}

# Calculate coherence by groups
coherence_data <- immune_data %>%
  mutate(HI_bin = cut(HI, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) %>%
  group_by(Sex, infection, HI_bin) %>%
  summarise(
    n = n(),
    mean_HI = mean(HI),
    mean_He = mean(He),
    coherence = ifelse(n >= 5,
                       calculate_coherence(pick(everything()), available_genes),
                       NA),
    .groups = "drop"
  ) %>%
  filter(!is.na(coherence))

# Model coherence with infection as factor
coherence_model <- lm(coherence ~ Sex * poly(mean_HI, 2) * infection,
                      data = coherence_data)

cat("Network Coherence Model:\n")
print(summary(coherence_model))

# Visualize coherence
plot_coherence <- ggplot(coherence_data,
                         aes(x = mean_HI, y = coherence, color = Sex)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  facet_wrap(~ infection, nrow = 1) +
  scale_color_manual(values = c("F" = "#4daf4a", "M" = "#ff7f00")) +
  labs(title = "Network Coherence Along Hybrid Gradient",
       x = "Hybrid Index",
       y = "Network Coherence") +
  theme_minimal()

save_plot_all_formats(plot_coherence, "network_coherence_by_infection")

# ==============================================================================
# 8. SUPPLEMENTARY ANALYSES
# ==============================================================================

cat("\n8. SUPPLEMENTARY ANALYSES\n")
cat("========================\n")

# A. UMAP based on immune genes (not PC scores)
cat("A. UMAP visualization...\n")

set.seed(42)
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1

umap_result <- umap(immune_matrix, config = umap_config)
immune_data$UMAP1 <- umap_result$layout[,1]
immune_data$UMAP2 <- umap_result$layout[,2]

# UMAP colored by predicted weight loss
plot_umap_weight_loss <- ggplot(immune_data, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(color = predicted_weight_loss, shape = Sex),
             size = 3, alpha = 0.8) +
  scale_color_viridis_c(name = "Predicted\nWeight Loss (%)") +
  facet_wrap(~ infection, nrow = 1) +
  labs(title = "UMAP of Immune Profiles Colored by Predicted Weight Loss") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"))

save_plot_all_formats(plot_umap_weight_loss, "UMAP_by_weight_loss")

# B. Gene-wise variance analysis (Levene's tests)
cat("\nB. Gene variance analysis...\n")

variance_results <- data.frame()

for(gene in available_genes) {
  levene_result <- leveneTest(
    as.formula(paste(gene, "~ Sex * infection * hybrid_category")),
    data = immune_data
  )

  variance_results <- rbind(variance_results, data.frame(
    Gene = gene,
    F_statistic = levene_result$`F value`[1],
    p_value = levene_result$`Pr(>F)`[1]
  ))
}

variance_results <- variance_results %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significant = p_adjusted < 0.05
  ) %>%
  arrange(p_value)

cat("Genes with differential variance:", sum(variance_results$significant), "\n")

# C. Multivariate dispersion (betadisper)
cat("\nC. Multivariate dispersion analysis...\n")

dist_matrix <- dist(scale(immune_matrix))
groups <- interaction(immune_data$Sex, immune_data$infection, immune_data$hybrid_category)

beta_disp <- betadisper(dist_matrix, groups)
disp_test <- permutest(beta_disp, permutations = 999)

cat("Dispersion test results:\n")
print(disp_test)

immune_data$dist_to_centroid <- beta_disp$distances

# Model dispersion
disp_model <- lm(dist_to_centroid ~ Sex * (HI + He) * infection, data = immune_data)
cat("\nDispersion model:\n")
print(summary(disp_model))

# D. MANOVA for trajectory differences
cat("\nD. MANOVA for UMAP trajectories...\n")

manova_result <- manova(cbind(UMAP1, UMAP2) ~ Sex * HI * He * infection,
                        data = immune_data)
cat("MANOVA Results:\n")
print(summary(manova_result, test = "Pillai"))

##################
cat("\n\nENHANCED UMAP ANALYSIS: HYBRID TRAJECTORIES\n")
cat("============================================\n")

# 1. Create trajectory data for UMAP space
trajectory_data <- immune_data %>%
  mutate(HI_bin = cut(HI, breaks = seq(0, 1, by = 0.2),
                      include.lowest = TRUE,
                      labels = c("0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"))) %>%
  group_by(Sex, infection, HI_bin) %>%
  summarise(
    n = n(),
    mean_HI = mean(HI),
    mean_He = mean(He),
    mean_UMAP1 = mean(UMAP1),
    mean_UMAP2 = mean(UMAP2),
    se_UMAP1 = sd(UMAP1)/sqrt(n()),
    se_UMAP2 = sd(UMAP2)/sqrt(n()),
    mean_weight_loss = mean(predicted_weight_loss),
    .groups = "drop"
  ) %>%
  filter(n >= 3)  # Only bins with sufficient samples

# 2. UMAP with hybrid trajectories by sex and infection
plot_umap_trajectories <- ggplot() +
  # Individual points in background
  geom_point(data = immune_data,
             aes(x = UMAP1, y = UMAP2, color = HI),
             alpha = 0.3, size = 1.5) +
  # Trajectory lines
  geom_path(data = trajectory_data,
            aes(x = mean_UMAP1, y = mean_UMAP2,
                group = interaction(Sex, infection)),
            size = 1.5, alpha = 0.8, color = "black") +
  # Trajectory points with error bars
  geom_errorbar(data = trajectory_data,
                aes(x = mean_UMAP1,
                    ymin = mean_UMAP2 - se_UMAP2,
                    ymax = mean_UMAP2 + se_UMAP2),
                width = 0.1, alpha = 0.5) +
  geom_errorbarh(data = trajectory_data,
                 aes(y = mean_UMAP2,
                     xmin = mean_UMAP1 - se_UMAP1,
                     xmax = mean_UMAP1 + se_UMAP1),
                 height = 0.1, alpha = 0.5) +
  # Trajectory points
  geom_point(data = trajectory_data,
             aes(x = mean_UMAP1, y = mean_UMAP2,
                 fill = mean_HI, shape = Sex, size = n),
             color = "black", stroke = 1) +
  # Facet by infection status
  facet_wrap(~ infection, scales = "free") +
  # Color scales
  scale_color_gradient2(low = "blue", mid = "purple", high = "red",
                        midpoint = 0.5, name = "Hybrid Index") +
  scale_fill_gradient2(low = "blue", mid = "purple", high = "red",
                       midpoint = 0.5, name = "Mean HI") +
  scale_shape_manual(values = c("F" = 21, "M" = 24), name = "Sex") +
  scale_size_continuous(range = c(3, 8), name = "Sample Size") +
  labs(title = "UMAP Trajectories: Sex-Specific Paths Along Hybrid Gradient",
       subtitle = "Arrows show progression from low to high hybrid index",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

# Add arrows to show direction
for(s in c("F", "M")) {
  for(inf in c("Uninfected", "Infected")) {
    traj_subset <- trajectory_data %>%
      filter(Sex == s, infection == inf) %>%
      arrange(mean_HI)

    if(nrow(traj_subset) > 1) {
      # Add arrow from second-to-last to last point
      n_points <- nrow(traj_subset)
      plot_umap_trajectories <- plot_umap_trajectories +
        annotate("segment",
                 x = traj_subset$mean_UMAP1[n_points-1],
                 y = traj_subset$mean_UMAP2[n_points-1],
                 xend = traj_subset$mean_UMAP1[n_points],
                 yend = traj_subset$mean_UMAP2[n_points],
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 color = ifelse(s == "F", "#E69F00", "#56B4E9"),
                 size = 1.5, alpha = 0.8)
    }
  }
}

print(plot_umap_trajectories)
save_plot_all_formats(plot_umap_trajectories, "UMAP_hybrid_trajectories")

# 3. UMAP colored by weight loss with trajectory overlay
plot_umap_weight_trajectory <- ggplot() +
  # Points colored by weight loss
  geom_point(data = immune_data,
             aes(x = UMAP1, y = UMAP2, color = predicted_weight_loss),
             size = 2, alpha = 0.6) +
  # Trajectory overlay
  geom_path(data = trajectory_data,
            aes(x = mean_UMAP1, y = mean_UMAP2,
                group = interaction(Sex, infection),
                linetype = Sex),
            size = 1.2, color = "black", alpha = 0.8) +
  facet_wrap(~ infection) +
  scale_color_viridis_c(name = "Predicted\nWeight Loss (%)",
                        option = "plasma") +
  scale_linetype_manual(values = c("F" = "solid", "M" = "dashed")) +
  labs(title = "Weight Loss Landscape in Immune Space",
       subtitle = "Black lines show sex-specific trajectories along hybrid gradient",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"))

print(plot_umap_weight_trajectory)
save_plot_all_formats(plot_umap_weight_trajectory, "UMAP_weight_loss_trajectories")

# 4. Statistical test of trajectory differences
cat("\nTesting trajectory differences between sexes:\n")
cat("=============================================\n")

# Calculate trajectory length for each sex/infection combination
trajectory_lengths <- trajectory_data %>%
  arrange(Sex, infection, mean_HI) %>%
  group_by(Sex, infection) %>%
  mutate(
    dist_to_next = sqrt((lead(mean_UMAP1) - mean_UMAP1)^2 +
                          (lead(mean_UMAP2) - mean_UMAP2)^2)
  ) %>%
  summarise(
    total_length = sum(dist_to_next, na.rm = TRUE),
    n_segments = sum(!is.na(dist_to_next)),
    .groups = "drop"
  )

cat("\nTrajectory lengths:\n")
print(trajectory_lengths)

# 5. Test if trajectories differ in direction
# Using angle analysis
trajectory_angles <- trajectory_data %>%
  arrange(Sex, infection, mean_HI) %>%
  group_by(Sex, infection) %>%
  mutate(
    # Calculate angle of movement
    dx = lead(mean_UMAP1) - mean_UMAP1,
    dy = lead(mean_UMAP2) - mean_UMAP2,
    angle = atan2(dy, dx) * 180/pi  # Convert to degrees
  ) %>%
  filter(!is.na(angle))

# Circular statistics for angle differences
library(circular)
angle_summary <- trajectory_angles %>%
  group_by(Sex, infection) %>%
  summarise(
    mean_angle = mean.circular(circular(angle, units = "degrees")),
    var_angle = var.circular(circular(angle, units = "degrees")),
    n_angles = n(),
    .groups = "drop"
  )

cat("\nTrajectory directions (mean angles):\n")
print(angle_summary)

# 6. Quantify divergence between male and female trajectories
divergence_data <- trajectory_data %>%
  dplyr::select(infection, HI_bin, Sex, mean_UMAP1, mean_UMAP2) %>%
  pivot_wider(names_from = Sex,
              values_from = c(mean_UMAP1, mean_UMAP2)) %>%
  mutate(
    divergence = sqrt((mean_UMAP1_F - mean_UMAP1_M)^2 +
                        (mean_UMAP2_F - mean_UMAP2_M)^2)
  ) %>%
  filter(!is.na(divergence))

plot_divergence <- ggplot(divergence_data,
                          aes(x = HI_bin, y = divergence,
                              fill = infection, group = infection)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Uninfected" = "#4daf4a",
                               "Infected" = "#ff7f00")) +
  labs(title = "Sex Divergence in Immune Space Along Hybrid Gradient",
       x = "Hybrid Index Bin",
       y = "Distance Between Male and Female Trajectories") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_divergence)
save_plot_all_formats(plot_divergence, "UMAP_sex_divergence")

# 7. Summary statistics
cat("\n\nSUMMARY OF TRAJECTORY ANALYSIS:\n")
cat("================================\n")

# Test if divergence increases with HI
divergence_model <- lm(divergence ~ as.numeric(HI_bin) * infection,
                       data = divergence_data)
cat("\nDivergence model:\n")
print(summary(divergence_model))

# Report key findings
max_divergence <- divergence_data %>%
  group_by(infection) %>%
  slice_max(divergence) %>%
  ungroup()

cat("\nMaximum sex divergence:\n")
print(max_divergence)

cat("\n✓ Enhanced UMAP analysis complete!\n")
cat("Key findings:\n")
cat("- Sex-specific trajectories visualized in immune space\n")
cat("- Divergence patterns quantified along hybrid gradient\n")
cat("- Weight loss landscape mapped to immune trajectories\n")

# ==============================================================================
# 9. CREATE SUMMARY TABLES
# ==============================================================================

cat("\n9. CREATING SUMMARY TABLES\n")
cat("=========================\n")

# Table 1: Model comparison
model_table <- model_comparison %>%
  gt() %>%
  tab_header(
    title = "Model Comparison for Predicted Weight Loss",
    subtitle = "Testing different interaction structures"
  ) %>%
  fmt_number(columns = AIC, decimals = 2) %>%
  tab_style(
    style = list(cell_fill(color = "lightblue")),
    locations = cells_body(rows = AIC == min(AIC))
  )

save_table_all_formats(model_table, "model_comparison")

pathway_table <- pathway_results %>%
  dplyr::select(Pathway, Sex_HI_padj, Sex_He_padj, HI_He_padj,
                Sex_infection_padj, R_squared) %>%
  # Ensure all numeric columns are actually numeric
  mutate(across(c(Sex_HI_padj, Sex_He_padj, HI_He_padj,
                  Sex_infection_padj, R_squared), as.numeric)) %>%
  gt() %>%
  tab_header(
    title = "Pathway Interaction Analysis",
    subtitle = "FDR-adjusted p-values"
  ) %>%
  fmt_number(columns = c(Sex_HI_padj, Sex_He_padj, HI_He_padj,
                         Sex_infection_padj, R_squared),
             decimals = 3) %>%
  tab_style(
    style = list(cell_fill(color = "lightgreen")),
    locations = cells_body(
      columns = c(Sex_HI_padj, Sex_He_padj, HI_He_padj, Sex_infection_padj),
      rows = (Sex_HI_padj < 0.05 | Sex_He_padj < 0.05 |
                HI_He_padj < 0.05 | Sex_infection_padj < 0.05)
    )
  )

save_table_all_formats(pathway_table, "pathway_interactions")

# ==============================================================================
# 10. SAVE RESULTS
# ==============================================================================

cat("\n10. SAVING ANALYSIS RESULTS\n")
cat("==========================\n")

# Save all results
saveRDS(list(
  immune_data = immune_data,
  models = list(
    best_model = best_model,
    all_models = list(model1, model2, model3, model4, model5)
  ),
  pathway_results = pathway_results,
  permanova = list(full = permanova_full, simple = permanova_simple),
  coherence = list(data = coherence_data, model = coherence_model),
  variance_results = variance_results,
  dispersion = list(betadisper = beta_disp, model = disp_model),
  umap = umap_result,
  manova = manova_result
), file = "results/immune_analysis_complete.rds")

cat("\n✓ Analysis complete!\n")
cat("All figures saved in results/figures/\n")
cat("All tables saved in results/tables/\n")
cat("Results saved to results/immune_analysis_complete.rds\n")
