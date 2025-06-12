# ==============================================================================
# IMMUNE LANDSCAPE ANALYSIS: Sex-Specific Trajectories in Hybrid Mice
# ==============================================================================
# Chapter 2: High-impact mechanistic analysis
# This script performs cutting-edge network and trajectory analyses to reveal
# sex-specific immune dysregulation in hybrid mice during infection
# ==============================================================================

cat("\n")
print_section("IMMUNE LANDSCAPE ANALYSIS")
cat("Revealing sex-specific immune trajectories in hybrid mice\n\n")

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("1. PREPARING IMMUNE GENE DATA\n")
cat("==============================\n")

# Extract immune gene expression data
immune_data <- field_mice %>%
  dplyr::select(Mouse_ID, HI, Sex, infection_status, predicted_weight_loss,
         all_of(immune_genes)) %>%
  filter(!is.na(HI) & !is.na(infection_status)) %>%
  mutate(
    # Create categorical hybrid status for visualization
    hybrid_category = case_when(
      HI < 0.1 ~ "M.m.domesticus",
      HI > 0.9 ~ "M.m.musculus",
      TRUE ~ "Hybrid"
    ),
    # Create interaction groups
    group = paste(Sex, infection_status, sep = "_")
  )

cat("✓ Immune data extracted: n =", nrow(immune_data), "mice\n")
cat("✓ Immune genes included:", length(immune_genes), "\n")

# Check for missing values in immune genes
missing_summary <- immune_data %>%
  dplyr::select(all_of(immune_genes)) %>%
  summarise_all(~sum(is.na(.))) %>%
  pivot_longer(everything(), names_to = "gene", values_to = "n_missing") %>%
  filter(n_missing > 0)

if(nrow(missing_summary) > 0) {
  cat("⚠ Missing values detected in", nrow(missing_summary), "genes\n")
  # Impute using median by group
  immune_data <- immune_data %>%
    group_by(Sex, infection_status) %>%
    mutate(across(all_of(immune_genes), ~ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
    ungroup()
  cat("✓ Missing values imputed using group medians\n")
}

# ==============================================================================
# 2. IMMUNE STATE SPACE TRAJECTORIES (UMAP)
# ==============================================================================

cat("\n2. CALCULATING IMMUNE TRAJECTORIES\n")
cat("===================================\n")

# Load required packages for UMAP
if (!require("umap", character.only = TRUE)) {
  install.packages("umap")
  library(umap)
}

# Prepare matrix for UMAP
immune_matrix <- immune_data %>%
  dplyr::select(all_of(immune_genes)) %>%
  as.matrix()

# Run UMAP
set.seed(42)
umap_config <- umap.defaults
umap_config$n_neighbors <- 15
umap_config$min_dist <- 0.1
umap_config$metric <- "euclidean"

cat("Running UMAP dimensionality reduction...\n")
umap_result <- umap(immune_matrix, config = umap_config)

# Add UMAP coordinates to data
immune_data$UMAP1 <- umap_result$layout[,1]
immune_data$UMAP2 <- umap_result$layout[,2]

cat("✓ UMAP completed\n")

# ==============================================================================
# 3. CREATE IMMUNE LANDSCAPE VISUALIZATION
# ==============================================================================

cat("\n3. CREATING IMMUNE LANDSCAPE PLOTS\n")
cat("===================================\n")

# Define color gradients for continuous hybrid index
hi_gradient <- colorRampPalette(c("#FF0000", "#800080", "#0000FF"))(100)

# A. Main landscape plot with trajectories
landscape_plot <- ggplot(immune_data, aes(x = UMAP1, y = UMAP2)) +
  # Add contour lines for density
  stat_density_2d(aes(group = group), color = "gray80", alpha = 0.5, bins = 5) +
  # Points colored by predicted weight loss
  geom_point(aes(color = predicted_weight_loss, shape = Sex),
             size = 3, alpha = 0.8) +
  scale_color_viridis_c(name = "Predicted\nWeight Loss (%)",
                        option = "plasma") +
  scale_shape_manual(values = c("F" = 16, "M" = 17)) +
  facet_wrap(~infection_status, nrow = 1) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "UMAP 1",
    y = "UMAP 2"
  )

landscape_plot

save_plot_all_formats(plot_object = landscape_plot,
                      plot_name = "Immune_landscape_trajectories_sex_specific_responses")

# B. Sex-specific trajectories along hybrid gradient
# Calculate mean positions for bins of hybrid index
trajectory_data <- immune_data %>%
  filter(infection_status == "TRUE") %>%
  mutate(hi_bin = cut(HI, breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)) %>%
  group_by(Sex, hi_bin) %>%
  summarise(
    mean_UMAP1 = mean(UMAP1),
    mean_UMAP2 = mean(UMAP2),
    mean_HI = mean(HI),
    mean_weight_loss = mean(predicted_weight_loss),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)  # Only show bins with sufficient samples

trajectory_plot <- ggplot() +
  # Individual points
  geom_point(data = filter(immune_data, infection_status == "Infected with Eimeria spp."),
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
  scale_color_gradient2(low = "#FF0000", mid = "#800080", high = "#0000FF",
                        midpoint = 0.5, name = "Hybrid Index") +
  scale_fill_gradient2(low = "#FF0000", mid = "#800080", high = "#0000FF",
                       midpoint = 0.5, name = "Hybrid Index") +
  scale_shape_manual(values = c("F" = 21, "M" = 24)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  # Add arrows to show direction
  geom_segment(data = trajectory_data %>%
                 group_by(Sex) %>%
                 slice(c(1, n())),
               aes(x = mean_UMAP1[1], y = mean_UMAP2[1],
                   xend = mean_UMAP1[2], yend = mean_UMAP2[2],
                   group = Sex),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               size = 1, alpha = 0.5)

trajectory_plot

save_plot_all_formats(plot_object = trajectory_plot,
                      plot_name = "Sex_specific_immune_trajectories")

# ==============================================================================
# 4. WEIGHTED GENE CO-EXPRESSION NETWORK ANALYSIS (WGCNA)
# ==============================================================================

cat("\n4. GENE CO-EXPRESSION NETWORK ANALYSIS\n")
cat("======================================\n")

# Install WGCNA if needed

if (!require("WGCNA", character.only = TRUE)) {
  BiocManager::install("WGCNA")
  library(WGCNA)
}



# Prepare data for WGCNA
infected_data <- immune_data %>%
  filter(infection_status == "TRUE")

# Run WGCNA separately for males and females
run_sex_specific_wgcna <- function(sex_filter) {
  cat(paste("\nAnalyzing", sex_filter, "mice...\n"))

  sex_data <- infected_data %>%
    filter(Sex == sex_filter) %>%
    dplyr::select(all_of(immune_genes))

  # Check sample size
  if(nrow(sex_data) < 20) {
    cat("⚠ Warning: Low sample size for WGCNA (n =", nrow(sex_data), ")\n")
  }

  # Calculate correlation matrix
  cor_matrix <- cor(as.matrix(sex_data), use = "pairwise.complete.obs")

  # Calculate network properties
  # Adjacency matrix (soft thresholding)
  adjacency <- adjacency(as.matrix(sex_data), power = 6)

  # Topological Overlap Matrix (TOM)
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM

  # Module identification
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 3)

  # Calculate module eigengenes
  MEList <- moduleEigengenes(as.matrix(sex_data), colors = dynamicMods)
  MEs <- MEList$eigengenes

  return(list(
    cor_matrix = cor_matrix,
    adjacency = adjacency,
    modules = dynamicMods,
    eigengenes = MEs,
    sex = sex_filter
  ))
}

# Run for both sexes
female_network <- run_sex_specific_wgcna("F")
male_network <- run_sex_specific_wgcna("M")

cat("✓ Sex-specific networks calculated\n")

# ==============================================================================
# 5. NETWORK PRESERVATION ANALYSIS
# ==============================================================================

cat("\n5. TESTING NETWORK PRESERVATION\n")
cat("================================\n")

# Calculate network preservation statistics
# Compare correlation matrices
cor_similarity <- cor(as.vector(female_network$cor_matrix),
                      as.vector(male_network$cor_matrix))

cat("Network correlation between sexes:", round(cor_similarity, 3), "\n")

# Identify genes with largest correlation differences
cor_diff <- female_network$cor_matrix - male_network$cor_matrix
max_diff_genes <- which(abs(cor_diff) > 0.5, arr.ind = TRUE)

if(nrow(max_diff_genes) > 0) {
  cat("✓ Found", nrow(max_diff_genes)/2, "gene pairs with correlation differences > 0.5\n")
}

# ==============================================================================
# 6. MACHINE LEARNING INTERPRETATION (SHAP VALUES)
# ==============================================================================

cat("\n6. CALCULATING GENE IMPORTANCE\n")
cat("===============================\n")

# We'll use the existing Random Forest model to calculate feature importance
# for sex-specific subsets

calculate_sex_importance <- function(sex_filter, infection_filter) {
  subset_data <- immune_data %>%
    filter(Sex == sex_filter,
           infection_status == infection_filter) %>%
    dplyr::select(all_of(immune_genes), predicted_weight_loss)

  if(nrow(subset_data) < 20) {
    return(NULL)
  }

  # Fit a simple linear model to get coefficients as proxy for importance
  formula <- as.formula(paste("predicted_weight_loss ~",
                              paste(immune_genes, collapse = " + ")))

  lm_model <- lm(formula, data = subset_data)

  # Extract standardized coefficients
  std_coefs <- coef(lm_model)[-1] * apply(subset_data[immune_genes], 2, sd, na.rm = TRUE)

  return(data.frame(
    gene = names(std_coefs),
    importance = as.numeric(std_coefs),
    sex = sex_filter,
    infection = infection_filter
  ))
}

# Calculate for all groups
importance_results <- bind_rows(
  calculate_sex_importance("F", "TRUE"),
  calculate_sex_importance("M", "TRUE"),
  calculate_sex_importance("F", "FALSE"),
  calculate_sex_importance("M", "FALSE")
)

cat("✓ Gene importance calculated for all groups\n")

# ==============================================================================
# 7. CREATE INTEGRATED FIGURE
# ==============================================================================

cat("\n7. CREATING INTEGRATED VISUALIZATION\n")
cat("====================================\n")

# A. Network heatmaps
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
                         midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold")) +
    labs(title = title, fill = "Correlation") +
    coord_fixed()
}

female_heatmap <- create_network_heatmap(female_network$cor_matrix,
                                         "Female Network")

female_heatmap

save_plot_all_formats(plot_object = female_heatmap,
                      plot_name = "Female_network_heatmap")

male_heatmap <- create_network_heatmap(male_network$cor_matrix,
                                       "Male Network")

male_heatmap

save_plot_all_formats(plot_object = male_heatmap,
                      plot_name = "male_heatmap")

# B. Gene importance plot
importance_plot <- importance_results %>%
  filter(infection == "TRUE") %>%
  mutate(gene = fct_reorder(gene, importance)) %>%
  ggplot(aes(x = importance, y = gene, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = sex_colors) +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank()) +
  labs(x = "Standardized Effect on Weight Loss",
       y = "Immune Gene",
       fill = "Sex",
       title = "Sex-Specific Gene Effects in Infected Mice")

importance_plot

save_plot_all_formats(plot_object = importance_plot,
                      plot_name = "Importance_plot_genes_sex")



# ==============================================================================
# 8. STATISTICAL TESTING OF TRAJECTORIES
# ==============================================================================

cat("\n8. STATISTICAL ANALYSIS\n")
cat("========================\n")

# Test if trajectories differ between sexes
# Using multivariate analysis of variance (MANOVA)
trajectory_test <- immune_data %>%
  filter(infection_status == "TRUE") %>%
  dplyr::select(Sex, HI, UMAP1, UMAP2)

# Test for interaction between Sex and HI on immune space position
manova_model <- manova(cbind(UMAP1, UMAP2) ~ Sex * HI, data = trajectory_test)
manova_summary <- summary(manova_model, test = "Pillai")

cat("MANOVA results for Sex × HI interaction on immune space:\n")
print(manova_summary)

# Follow-up univariate tests
univariate_umap1 <- lm(UMAP1 ~ Sex * HI, data = trajectory_test)
univariate_umap2 <- lm(UMAP2 ~ Sex * HI, data = trajectory_test)

cat("\nUnivariate tests:\n")
cat("UMAP1 - Sex×HI interaction p-value:",
    round(summary(univariate_umap1)$coefficients["SexM:HI", "Pr(>|t|)"], 4), "\n")
cat("UMAP2 - Sex×HI interaction p-value:",
    round(summary(univariate_umap2)$coefficients["SexM:HI", "Pr(>|t|)"], 4), "\n")

# ==============================================================================
# 9. EXPORT RESULTS
# ==============================================================================

cat("\n9. EXPORTING RESULTS\n")
cat("====================\n")

# Save processed data
write_csv(immune_data,
          file.path("results", "tables", "immune_landscape_data.csv"))

# Save network statistics
network_stats <- data.frame(
  metric = c("Network correlation", "Infected females (n)", "Infected males (n)",
             "Unique modules (F)", "Unique modules (M)"),
  value = c(round(cor_similarity, 3),
            sum(infected_data$Sex == "F"),
            sum(infected_data$Sex == "M"),
            length(unique(female_network$modules)),
            length(unique(male_network$modules)))
)

write_csv(network_stats,
          file.path("results", "tables", "network_statistics.csv"))

# Save importance results
write_csv(importance_results,
          file.path("results", "tables", "gene_importance_by_group.csv"))

