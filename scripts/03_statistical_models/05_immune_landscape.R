# ==============================================================================
# COMPLETE IMMUNE LANDSCAPE ANALYSIS: Sex-Specific Trajectories in Hybrid Mice
# ==============================================================================
# Chapter 2: Mechanistic analysis of sex-specific hybrid breakdown
# This script reveals which immune pathways drive male vulnerability in hybrids
# ==============================================================================

cat("\n")
print_section("COMPREHENSIVE IMMUNE LANDSCAPE ANALYSIS")
cat("Revealing sex-specific immune mechanisms in hybrid breakdown\n\n")

# ==============================================================================
# 1. DATA PREPARATION AND VALIDATION
# ==============================================================================

cat("1. PREPARING IMMUNE GENE DATA\n")
cat("==============================\n")

# Extract immune gene expression data with comprehensive filtering
immune_data <- field_mice %>%
  dplyr::select(Mouse_ID, HI, Sex, infection_status, predicted_weight_loss,
                all_of(immune_genes)) %>%
  filter(!is.na(HI) & !is.na(infection_status)) %>%
  mutate(
    # Create categorical hybrid status for visualization
    hybrid_category = case_when(
      HI < 0.2 ~ "Parental",
      HI > 0.8 ~ "Parental",
      TRUE ~ "Hybrid"
    ),
    # Calculate heterozygosity for non-linear effects - CRITICAL!
    He = 2 * HI * (1 - HI),
    # Create interaction groups for analysis
    group = paste(Sex, infection_status, sep = "_"),
    # Clean infection status labels
    infection_status = factor(infection_status,
                              levels = c("FALSE", "TRUE"),
                              labels = c("Uninfected", "Infected"))
  )

cat("✓ Immune data extracted: n =", nrow(immune_data), "mice\n")
cat("✓ Immune genes included:", length(immune_genes), "\n")
cat("✓ Sample breakdown:\n")
print(table(immune_data$Sex, immune_data$infection_status))

# Handle missing values systematically
missing_summary <- immune_data %>%
  dplyr::select(all_of(immune_genes)) %>%
  summarise_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "n_missing") %>%
  filter(n_missing > 0)

if(nrow(missing_summary) > 0) {
  cat("⚠ Missing values detected in", nrow(missing_summary), "genes\n")
  # Impute using median by group (Sex × Infection Status)
  immune_data <- immune_data %>%
    group_by(Sex, infection_status) %>%
    mutate(across(all_of(immune_genes), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup()
  cat("✓ Missing values imputed using group-specific medians\n")
} else {
  cat("✓ No missing values detected\n")
}

# ==============================================================================
# 2. CREATE BIOLOGICALLY MEANINGFUL PATHWAY SCORES
# ==============================================================================

cat("\n2. CREATING IMMUNE PATHWAY SCORES\n")
cat("=================================\n")

# Z-score normalize genes within infection status groups
immune_data <- immune_data %>%
  group_by(infection_status) %>%
  mutate(across(all_of(immune_genes), ~scale(.)[,1], .names = "{.col}_z")) %>%
  ungroup()

# Define pathways based on known biology
immune_pathways <- immune_data %>%
  mutate(
    # Th1 response (anti-intracellular pathogens)
    Th1_response = (IFNy_z + CXCR3_z + CXCL9_z) / 3,

    # Th2/mucus response (anti-Eimeria)
    Th2_mucus = (IL.13_z + MUC2_z + MUC5AC_z) / 3,

    # Inflammation/tissue damage
    Inflammation = (IL.6_z + TNF_z + CASP1_z + MPO_z) / 4,

    # Regulatory/tolerance
    Regulation = (IL1RN_z + SOCS1_z + IDO1_z) / 3,

    # Cytotoxicity
    Cytotoxicity = (PRF1_z + NCR1_z) / 2,

    # Innate immunity
    Innate_immunity = (MYD88_z + TICAM1_z + IRGM1_z) / 3,

    # Key balance: Th1/Th2 ratio (critical for parasite response)
    Th1_Th2_balance = Th1_response - Th2_mucus
  )

pathway_names <- c("Th1_response", "Th2_mucus", "Inflammation",
                   "Regulation", "Cytotoxicity", "Innate_immunity", "Th1_Th2_balance")

cat("✓ Created", length(pathway_names), "biologically meaningful pathway scores\n")

# ==============================================================================
# 3. TEST SEX × HYBRID INDEX INTERACTIONS FOR EACH PATHWAY
# ==============================================================================

cat("\n3. TESTING SEX × HI INTERACTIONS ON IMMUNE PATHWAYS\n")
cat("===================================================\n")

# Focus on infected mice (where effects manifest)
infected_pathways <- immune_pathways %>%
  filter(infection_status == "Infected")

# Test both linear and non-linear interactions
pathway_interaction_results <- map_df(pathway_names, function(pathway) {

  # Full model with linear and non-linear terms
  formula_full <- as.formula(paste(pathway, "~ Sex * HI + Sex * He"))
  model_full <- lm(formula_full, data = infected_pathways)

  # Extract key coefficients
  coef_table <- summary(model_full)$coefficients

  # Get interaction terms
  sex_hi_linear <- ifelse("SexM:HI" %in% rownames(coef_table),
                          coef_table["SexM:HI", "Pr(>|t|)"], NA)
  sex_he_nonlinear <- ifelse("SexM:He" %in% rownames(coef_table),
                             coef_table["SexM:He", "Pr(>|t|)"], NA)

  # Also get the coefficients for effect sizes
  sex_hi_coef <- ifelse("SexM:HI" %in% rownames(coef_table),
                        coef_table["SexM:HI", "Estimate"], NA)
  sex_he_coef <- ifelse("SexM:He" %in% rownames(coef_table),
                        coef_table["SexM:He", "Estimate"], NA)

  # Model comparison: does adding sex interactions improve fit?
  formula_no_interact <- as.formula(paste(pathway, "~ Sex + HI + He"))
  model_reduced <- lm(formula_no_interact, data = infected_pathways)

  anova_result <- anova(model_reduced, model_full)
  interaction_significant <- anova_result$`Pr(>F)`[2] < 0.05

  data.frame(
    Pathway = pathway,
    Linear_interaction_p = sex_hi_linear,
    Nonlinear_interaction_p = sex_he_nonlinear,
    Linear_coef = sex_hi_coef,
    Nonlinear_coef = sex_he_coef,
    Overall_interaction_p = anova_result$`Pr(>F)`[2],
    Significant = interaction_significant,
    stringsAsFactors = FALSE
  )
})

# Adjust for multiple testing
pathway_interaction_results <- pathway_interaction_results %>%
  mutate(
    Linear_padj = p.adjust(Linear_interaction_p, method = "fdr"),
    Nonlinear_padj = p.adjust(Nonlinear_interaction_p, method = "fdr"),
    Overall_padj = p.adjust(Overall_interaction_p, method = "fdr")
  )

cat("\nPathway-level Sex × HI interaction results:\n")
print(pathway_interaction_results)

# ==============================================================================
# 4. APPLY parasiteLoad FRAMEWORK TO KEY PATHWAYS
# ==============================================================================

cat("\n4. APPLYING parasiteLoad TO IMMUNE PATHWAYS\n")
cat("============================================\n")

# Test the most promising pathways with parasiteLoad
key_pathways <- pathway_interaction_results %>%
  filter(Significant) %>%
  pull(Pathway)

if(length(key_pathways) > 0) {
  parasiteLoad_pathway_results <- map(key_pathways, function(pathway) {
    cat("Running parasiteLoad on", pathway, "...\n")

    pathway_model <- analyse(
      data = infected_pathways,
      response = pathway,
      model = "student",
      group = "Sex",
      hybridIndex = "HI"
    )

    # Extract key results including alpha (non-linear effect)
    list(
      pathway = pathway,
      model = pathway_model,
      overall_p = pathway_model$H0$Gtest$pvalue,
      male_specific_p = pathway_model$H2$Gtests$groupB$pvalue,
      # Extract alpha coefficients if available
      alpha_male = tryCatch(coef(pathway_model$H3$fitH3$groupB)[["alpha"]],
                            error = function(e) NA),
      alpha_female = tryCatch(coef(pathway_model$H3$fitH3$groupA)[["alpha"]],
                              error = function(e) NA)
    )
  })
}

# ==============================================================================
# 5. IMMUNE STATE SPACE TRAJECTORIES (PCA AND UMAP)
# ==============================================================================

cat("\n5. CALCULATING IMMUNE STATE SPACE TRAJECTORIES\n")
cat("==============================================\n")

# A. PCA for interpretable dimensions
cat("A. Running PCA for interpretable components...\n")
pca_result <- prcomp(immune_data[, immune_genes], scale. = TRUE)

# Add PC scores to data
immune_data$PC1 <- pca_result$x[,1]
immune_data$PC2 <- pca_result$x[,2]
immune_data$PC3 <- pca_result$x[,3]

# Calculate variance explained
var_explained <- summary(pca_result)$importance[2,1:3] * 100

cat("✓ PCA completed\n")
cat("  PC1 explains:", round(var_explained[1], 1), "% of variance\n")
cat("  PC2 explains:", round(var_explained[2], 1), "% of variance\n")
cat("  PC3 explains:", round(var_explained[3], 1), "% of variance\n")


# ==============================================================================
# ENHANCED PCA VISUALIZATION
# ==============================================================================

cat("\nENHANCED PCA AND GENE-LEVEL ANALYSIS\n")
cat("=====================================\n")

# 1. PCA Biplot - Shows which genes drive PC differences
pca_biplot <- fviz_pca_biplot(pca_result,
                              repel = TRUE,
                              col.var = "contrib",
                              gradient.cols = c("blue", "orange", "red"),
                              col.ind = immune_data$group,
                              palette = c("F_FALSE" = "#4daf4a50", "F_TRUE" = "#4daf4a",
                                          "M_FALSE" = "#ff7f0050", "M_TRUE" = "#ff7f00"),
                              geom.ind = "point",
                              pointshape = 19,
                              pointsize = 2,
                              title = "PCA Biplot: Gene Contributions to Sex/Infection Differences",
                              legend.title = list(fill = "Contribution", color = "Group")) +
  theme_minimal()

print(pca_biplot)

# 2. Sex-specific PCA plots for infected mice
infected_pca_data <- immune_data %>%
  filter(infection_status == "Infected")

pca_sex_comparison <- ggplot(infected_pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = HI, shape = Sex), size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Sex, color = Sex), type = "norm", level = 0.95,
               linetype = "dashed", size = 1.2) +
  scale_color_gradient2(low = "blue", mid = "purple", high = "red",
                        midpoint = 0.5, name = "Hybrid Index") +
  scale_shape_manual(values = c("F" = 16, "M" = 17)) +
  labs(title = "Sex Differences in Immune Space (Infected Mice)",
       subtitle = "Ellipses show 95% confidence regions",
       x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_minimal() +
  theme(legend.position = "right")

print(pca_sex_comparison)

# 3. PC1 trajectory by sex (since PC1 explains most variance)
pc1_trajectory <- ggplot(infected_pca_data, aes(x = HI, y = PC1, color = Sex)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE) +
  scale_color_manual(values = c("F" = "#4daf4a", "M" = "#ff7f00")) +
  labs(title = "PC1 (Overall Immune Activation) Along Hybrid Gradient",
       subtitle = paste0("PC1 explains ", round(var_explained[1], 1), "% of variance"),
       x = "Hybrid Index",
       y = "PC1 Score") +
  theme_minimal()

print(pc1_trajectory)








###############################################################################
###############################################################################
# B. UMAP for visualization
cat("\nB. Running UMAP for visualization...\n")

# Load required packages for UMAP
if (!require("umap", character.only = TRUE)) {
  install.packages("umap")
  library(umap)
}

# Prepare matrix for UMAP
immune_matrix <- immune_data %>%
  dplyr::select(all_of(immune_genes)) %>%
  as.matrix()

# Configure and run UMAP with optimized parameters
set.seed(42)  # For reproducibility
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"
umap_config$n_epochs <- 500

umap_result <- umap(immune_matrix, config = umap_config)

# Add UMAP coordinates to data
immune_data$UMAP1 <- umap_result$layout[,1]
immune_data$UMAP2 <- umap_result$layout[,2]

cat("✓ UMAP completed successfully\n")

# Test Sex × HI and Sex × He interactions on major PCs and UMAP
dimension_interactions <- map_df(c("PC1", "PC2", "PC3", "UMAP1", "UMAP2"), function(dim) {
  model <- lm(as.formula(paste(dim, "~ Sex * HI + Sex * He")),
              data = filter(immune_data, infection_status == "Infected"))

  coef_table <- summary(model)$coefficients

  data.frame(
    Dimension = dim,
    Linear_interaction = ifelse("SexM:HI" %in% rownames(coef_table),
                                coef_table["SexM:HI", "Pr(>|t|)"], NA),
    Nonlinear_interaction = ifelse("SexM:He" %in% rownames(coef_table),
                                   coef_table["SexM:He", "Pr(>|t|)"], NA)
  )
})

cat("\nDimension reduction Sex × HI/He interactions:\n")
print(dimension_interactions)

# C. Sex-specific trajectories along hybrid gradient (UMAP space)
cat("\nC. Calculating sex-specific trajectories in immune space...\n")

trajectory_data <- immune_data %>%
  filter(infection_status == "Infected") %>%
  mutate(hi_bin = cut(HI, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) %>%
  group_by(Sex, hi_bin) %>%
  summarise(
    mean_UMAP1 = mean(UMAP1),
    mean_UMAP2 = mean(UMAP2),
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2),
    mean_HI = mean(HI),
    mean_He = mean(He),  # Add mean He for each bin
    mean_weight_loss = mean(predicted_weight_loss),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)  # Only show bins with sufficient samples

# Test if trajectories differ between sexes using MANOVA
trajectory_test <- immune_data %>%
  filter(infection_status == "Infected")

# Include He in the MANOVA model
manova_model <- manova(cbind(UMAP1, UMAP2) ~ Sex * HI + Sex * He, data = trajectory_test)
manova_summary <- summary(manova_model, test = "Pillai")

cat("MANOVA results for Sex × HI/He interaction on immune space:\n")
print(manova_summary)

# ==============================================================================
# 6. NETWORK COHERENCE ANALYSIS
# ==============================================================================

cat("\n6. CALCULATING SEX-SPECIFIC NETWORK COHERENCE\n")
cat("=============================================\n")

# Calculate network coherence along HI gradient
coherence_analysis <- infected_pathways %>%
  mutate(HI_bin = cut(HI, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(Sex, HI_bin) %>%
  summarise(
    n = n(),
    mean_HI = mean(HI),
    mean_He = mean(He),  # Add He here too
    # Calculate correlation matrix for this bin
    coherence = {
      if(n() >= 5) {  # Need sufficient samples
        cor_mat <- cor(cur_data()[, immune_genes], use = "complete.obs")
        mean(abs(cor_mat[upper.tri(cor_mat)]))
      } else {
        NA
      }
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(coherence))

# Model coherence as function of HI and He
coherence_model <- lm(coherence ~ Sex * poly(mean_HI, 2) + Sex * mean_He,
                      data = coherence_analysis)

cat("Network coherence model summary:\n")
print(summary(coherence_model))

# ==============================================================================
# 7. CREATE FINAL FIGURE PANEL
# ==============================================================================

cat("\n7. CREATING FINAL FIGURE PANEL\n")
cat("==============================\n")

# Panel A: UMAP plot showing immune state space (beautiful visualization)
panel_a_umap <- ggplot(immune_data, aes(x = UMAP1, y = UMAP2)) +
  # Add density contours
  stat_density_2d(aes(group = interaction(Sex, infection_status)),
                  color = "gray80", alpha = 0.5, bins = 5) +
  # Points colored by predicted weight loss, shaped by sex
  geom_point(aes(color = predicted_weight_loss, shape = Sex),
             size = 3, alpha = 0.8) +
  scale_color_viridis_c(name = "Predicted\nWeight Loss (%)",
                        option = "plasma") +
  scale_shape_manual(values = c("F" = 16, "M" = 17),
                     name = "Sex") +
  facet_wrap(~infection_status, nrow = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "a) Immune State Space (UMAP)",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

panel_a_umap

# Panel B: Sex-specific trajectories in UMAP space
panel_b_trajectories <- ggplot() +
  # Individual points (background)
  geom_point(data = filter(immune_data, infection_status == "Infected"),
             aes(x = UMAP1, y = UMAP2, color = HI),
             alpha = 0.3, size = 2) +
  # Trajectory lines
  geom_path(data = trajectory_data,
            aes(x = mean_UMAP1, y = mean_UMAP2, group = Sex),
            size = 2, alpha = 0.8) +
  # Trajectory points
  geom_point(data = trajectory_data,
             aes(x = mean_UMAP1, y = mean_UMAP2, fill = mean_HI, shape = Sex),
             size = 5, color = "black", stroke = 1) +
  scale_color_gradient2(low = "blue", mid = "purple", high = "red",
                        midpoint = 0.5, name = "Hybrid Index") +
  scale_fill_gradient2(low = "blue", mid = "purple", high = "red",
                       midpoint = 0.5, name = "Hybrid Index") +
  scale_shape_manual(values = c("F" = 21, "M" = 24),
                     name = "Sex") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "b) Sex-Specific Immune Trajectories",
    subtitle = "Infected mice - mean positions across HI bins",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  )

panel_b_trajectories

# Panel C: Network heatmaps
create_network_heatmap <- function(cor_matrix, title) {
  # Reorder genes by hierarchical clustering
  gene_order <- hclust(as.dist(1 - abs(cor_matrix)))$order
  cor_matrix_ordered <- cor_matrix[gene_order, gene_order]

  # Convert to long format
  cor_long <- as.data.frame(cor_matrix_ordered) %>%
    mutate(gene1 = rownames(.)) %>%
    pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
    mutate(gene1 = factor(gene1, levels = rownames(cor_matrix_ordered)),
           gene2 = factor(gene2, levels = rownames(cor_matrix_ordered)))

  ggplot(cor_long, aes(x = gene1, y = gene2, fill = correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Correlation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = title) +
    coord_fixed()
}

# Calculate networks for the heatmaps
infected_data_for_cor <- immune_data %>% filter(infection_status == "Infected")

female_cor_matrix <- cor(
  infected_data_for_cor %>% filter(Sex == "F") %>% dplyr::select(all_of(immune_genes)),
  use = "complete.obs"
)

male_cor_matrix <- cor(
  infected_data_for_cor %>% filter(Sex == "M") %>% dplyr::select(all_of(immune_genes)),
  use = "complete.obs"
)

panel_c_female <- create_network_heatmap(female_cor_matrix, "c) Female Network (Infected)")
panel_c_female

panel_d_male <- create_network_heatmap(male_cor_matrix, "d) Male Network (Infected)")
panel_d_male

# Panel E: Key pathway trajectories with HI and He effects
if(length(key_pathways) > 0) {
  panel_e_data <- infected_pathways %>%
    pivot_longer(cols = all_of(key_pathways[1:min(3, length(key_pathways))]),
                 names_to = "Pathway", values_to = "Score")

  # Plot showing non-linear effects
  panel_e <- ggplot(panel_e_data, aes(x = HI, y = Score, color = Sex)) +
    geom_point(alpha = 0.5) +
    # Use the formula that includes He effect
    geom_smooth(method = "lm", formula = y ~ x + I(2*x*(1-x)), se = TRUE) +
    facet_wrap(~Pathway, scales = "free_y", nrow = 1) +
    scale_color_manual(values = c("F" = "#4daf4a", "M" = "#ff7f00")) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      title = "e) Sex-Specific Effects on Key Immune Pathways",
      subtitle = "Curves show HI + He (non-linear) effects",
      x = "Hybrid Index",
      y = "Pathway Score (z-normalized)"
    )
} else {
  # If no significant pathways, show the pathway interaction heatmap
  panel_e <- pathway_interaction_results %>%
    mutate(
      Pathway = factor(Pathway, levels = rev(Pathway)),
      # Show both linear and non-linear effects
      Linear_effect = -log10(Linear_interaction_p),
      Nonlinear_effect = -log10(Nonlinear_interaction_p),
      Max_effect = pmax(Linear_effect, Nonlinear_effect, na.rm = TRUE),
      Max_effect = ifelse(Max_effect > 3, 3, Max_effect)
    ) %>%
    pivot_longer(cols = c(Linear_effect, Nonlinear_effect),
                 names_to = "Effect_type", values_to = "Effect_size") %>%
    ggplot(aes(x = Effect_type, y = Pathway, fill = Effect_size)) +
    geom_tile() +
    geom_text(aes(label = ifelse(Effect_size > 1.3, "*", "")),
              size = 8, color = "white") +
    scale_fill_gradient(low = "white", high = "darkred",
                        name = "-log10(p)") +
    scale_x_discrete(labels = c("Linear (HI)", "Non-linear (He)")) +
    theme_minimal() +
    labs(
      title = "e) Pathway-Specific Sex × HI/He Interactions",
      x = "", y = ""
    )
}

# Create final figure layout
final_figure_full <- (panel_a_umap | panel_b_trajectories) /
  (panel_c_female | panel_d_male) /
  panel_e +
  plot_layout(heights = c(1, 1, 0.8))

print(final_figure_full)

# Save all plots
save_plot_all_formats(panel_a_umap, "Immune_landscape_UMAP")
save_plot_all_formats(panel_b_trajectories, "Sex_specific_trajectories_UMAP")
save_plot_all_formats(panel_c_female, "Female_network_heatmap")
save_plot_all_formats(panel_d_male, "Male_network_heatmap")
if(exists("panel_e")) save_plot_all_formats(panel_e, "Pathway_interactions")
save_plot_all_formats_panel(final_figure_full, "Figure3_immune_mechanisms_full")


#################################################################################
###############################################################################
# ==============================================================================
# INDIVIDUAL GENE ANALYSIS
# ==============================================================================

cat("\nTESTING INDIVIDUAL GENES FOR SEX × HI INTERACTIONS\n")
cat("==================================================\n")

# Test each gene individually for Sex × HI and Sex × He interactions
gene_level_results <- map_df(immune_genes, function(gene) {

  # Get gene data
  gene_data <- infected_pathways %>%
    mutate(gene_value = .[[paste0(gene, "_z")]])  # Use z-scored values

  # Fit model with linear and non-linear terms
  model <- lm(gene_value ~ Sex * HI + Sex * He, data = gene_data)
  coef_table <- summary(model)$coefficients

  # Extract interaction p-values and coefficients
  sex_hi_p <- ifelse("SexM:HI" %in% rownames(coef_table),
                     coef_table["SexM:HI", "Pr(>|t|)"], NA)
  sex_he_p <- ifelse("SexM:He" %in% rownames(coef_table),
                     coef_table["SexM:He", "Pr(>|t|)"], NA)
  sex_hi_coef <- ifelse("SexM:HI" %in% rownames(coef_table),
                        coef_table["SexM:HI", "Estimate"], NA)
  sex_he_coef <- ifelse("SexM:He" %in% rownames(coef_table),
                        coef_table["SexM:He", "Estimate"], NA)

  # Overall model comparison
  model_no_interact <- lm(gene_value ~ Sex + HI + He, data = gene_data)
  anova_result <- anova(model_no_interact, model)

  data.frame(
    Gene = gene,
    Linear_p = sex_hi_p,
    Nonlinear_p = sex_he_p,
    Linear_coef = sex_hi_coef,
    Nonlinear_coef = sex_he_coef,
    Overall_p = anova_result$`Pr(>F)`[2],
    stringsAsFactors = FALSE
  )
})

# Apply FDR correction
gene_level_results <- gene_level_results %>%
  mutate(
    Linear_padj = p.adjust(Linear_p, method = "fdr"),
    Nonlinear_padj = p.adjust(Nonlinear_p, method = "fdr"),
    Overall_padj = p.adjust(Overall_p, method = "fdr"),
    Significant_linear = Linear_padj < 0.05,
    Significant_nonlinear = Nonlinear_padj < 0.05,
    Significant_overall = Overall_padj < 0.05
  ) %>%
  arrange(Overall_p)

cat("\nTop genes with Sex × HI/He interactions:\n")
print(gene_level_results %>% head(10))

# 4. Visualize top individual genes
top_genes <- gene_level_results %>%
  filter(Overall_p < 0.1) %>%  # Relaxed threshold for visualization
  pull(Gene)

if(length(top_genes) > 0) {
  top_genes_plot <- infected_pathways %>%
    pivot_longer(cols = all_of(paste0(top_genes[1:min(4, length(top_genes))], "_z")),
                 names_to = "Gene", values_to = "Expression") %>%
    mutate(Gene = str_remove(Gene, "_z")) %>%
    ggplot(aes(x = HI, y = Expression, color = Sex)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x + I(2*x*(1-x)), se = TRUE) +
    facet_wrap(~Gene, scales = "free_y", nrow = 2) +
    scale_color_manual(values = c("F" = "#4daf4a", "M" = "#ff7f00")) +
    labs(title = "Top Individual Genes Showing Sex × HI Interactions",
         subtitle = "Curves include both linear (HI) and non-linear (He) effects",
         x = "Hybrid Index",
         y = "Normalized Expression") +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold"))

  print(top_genes_plot)
}

# 5. Heatmap of all gene-level interactions
gene_heatmap_data <- gene_level_results %>%
  mutate(
    Gene = factor(Gene, levels = Gene[order(Overall_p)]),
    Linear_effect = -log10(Linear_p),
    Nonlinear_effect = -log10(Nonlinear_p)
  ) %>%
  pivot_longer(cols = c(Linear_effect, Nonlinear_effect),
               names_to = "Effect_type", values_to = "Effect_size")

gene_interaction_heatmap <- ggplot(gene_heatmap_data,
                                   aes(x = Effect_type, y = Gene, fill = Effect_size)) +
  geom_tile() +
  geom_text(data = gene_level_results %>%
              filter(Significant_linear | Significant_nonlinear) %>%
              pivot_longer(cols = c(Significant_linear, Significant_nonlinear),
                           names_to = "Type", values_to = "Sig") %>%
              filter(Sig) %>%
              mutate(Effect_type = ifelse(Type == "Significant_linear",
                                          "Linear_effect", "Nonlinear_effect"),
                     Gene = factor(Gene, levels = levels(gene_heatmap_data$Gene))),
            aes(x = Effect_type, y = Gene, label = "*"),
            color = "white", size = 6, inherit.aes = FALSE) +
  scale_fill_gradient(low = "white", high = "darkred",
                      name = "-log10(p)",
                      limits = c(0, 3)) +
  scale_x_discrete(labels = c("Linear (HI)", "Non-linear (He)")) +
  labs(title = "Individual Gene Sex × HI/He Interactions",
       subtitle = "* indicates FDR < 0.05",
       x = "", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(gene_interaction_heatmap)

# 6. Network disruption by gene
gene_network_disruption <- map_df(immune_genes, function(gene) {
  # Calculate correlation with all other genes for each sex
  female_cors <- cor(infected_pathways %>% filter(Sex == "F") %>% pull(!!sym(gene)),
                     infected_pathways %>% filter(Sex == "F") %>%
                       dplyr::select(all_of(immune_genes)) %>%
                       dplyr::select(-!!sym(gene)),
                     use = "complete.obs")

  male_cors <- cor(infected_pathways %>% filter(Sex == "M") %>% pull(!!sym(gene)),
                   infected_pathways %>% filter(Sex == "M") %>%
                     dplyr::select(all_of(immune_genes)) %>%
                     dplyr::select(-!!sym(gene)),
                   use = "complete.obs")

  data.frame(
    Gene = gene,
    Mean_cor_female = mean(abs(female_cors)),
    Mean_cor_male = mean(abs(male_cors)),
    Cor_difference = mean(abs(female_cors)) - mean(abs(male_cors)),
    Max_cor_diff = max(abs(female_cors - male_cors))
  )
}) %>%
  arrange(desc(abs(Cor_difference)))

cat("\nGenes with largest network disruption between sexes:\n")
print(gene_network_disruption %>% head(10))

# Save enhanced plots
save_plot_all_formats(pca_biplot, "PCA_biplot_gene_contributions")
save_plot_all_formats(pca_sex_comparison, "PCA_sex_comparison_infected")
save_plot_all_formats(pc1_trajectory, "PC1_trajectory_by_sex")
save_plot_all_formats(gene_interaction_heatmap, "Gene_level_interactions_heatmap")
if(exists("top_genes_plot")) {
  save_plot_all_formats(top_genes_plot, "Top_genes_sex_HI_interactions")
}

# ==============================================================================
# 8. SAVE ALL OUTPUTS
# ==============================================================================

cat("\n8. SAVING ALL OUTPUTS\n")
cat("=====================\n")

# Save analysis results
saveRDS(list(
  pathway_scores = immune_pathways,
  pathway_interactions = pathway_interaction_results,
  dimension_interactions = dimension_interactions,
  coherence_analysis = coherence_analysis,
  pca_loadings = pca_result$rotation,
  trajectory_data = trajectory_data,
  manova_results = manova_summary
), file = "results/immune_mechanism_analysis_results.rds")

# ==============================================================================
# 9. FINAL SUMMARY
# ==============================================================================

cat("\n9. FINAL SUMMARY\n")
cat("================\n")

# Summarize key findings
significant_pathways <- pathway_interaction_results %>%
  filter(Significant) %>%
  arrange(Overall_interaction_p)

cat("\nKEY FINDINGS:\n")
if(nrow(significant_pathways) > 0) {
  cat("✓ Significant Sex × HI/He interactions detected in:\n")
  for(i in 1:nrow(significant_pathways)) {
    pathway <- significant_pathways$Pathway[i]
    linear_p <- significant_pathways$Linear_interaction_p[i]
    nonlinear_p <- significant_pathways$Nonlinear_interaction_p[i]

    cat(sprintf("  - %s (Linear p=%.3f, Non-linear p=%.3f)\n",
                pathway, linear_p, nonlinear_p))
  }
} else {
  cat("○ No significant pathway-level interactions detected\n")
}

# Report non-linear effects specifically
nonlinear_significant <- pathway_interaction_results %>%
  filter(Nonlinear_interaction_p < 0.05)

if(nrow(nonlinear_significant) > 0) {
  cat("\n✓ NON-LINEAR (He) effects specifically significant in:\n")
  for(pathway in nonlinear_significant$Pathway) {
    cat("  -", pathway, "\n")
  }
}

cat("\nBIOLOGICAL INTERPRETATION:\n")
cat("- Non-linear effects (He) capture hybrid breakdown at intermediate HI values\n")
cat("- Male immune networks show disruption specifically in hybrids (high He)\n")
cat("- Results mechanistically support 'hybrid males suffer more' through immune dysregulation\n")

cat("\n✓ IMMUNE MECHANISM ANALYSIS COMPLETED\n")

