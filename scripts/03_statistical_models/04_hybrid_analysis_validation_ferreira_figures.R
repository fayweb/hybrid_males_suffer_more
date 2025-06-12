# ==============================================================================
# FERREIRA VALIDATION - FIGURE GENERATION
# ==============================================================================
# Purpose: Create publication-ready figures for Ferreira methodology validation
# to complement parasiteLoad banana plots in Figure 2
#
# Creates: Effect plots, validation comparison panels, and combined figures
# Style: Matching parasiteLoad aesthetic for seamless integration
# ==============================================================================

cat("\n=== FERREIRA VALIDATION FIGURE GENERATION ===\n")
cat("Creating publication-ready figures for Figure 2 panel\n")
cat("Matching parasiteLoad aesthetic for seamless integration\n\n")

# ==============================================================================
# 1. FIGURE SETUP & THEME CONFIGURATION
# ==============================================================================

cat("1. SETTING UP FIGURE THEMES\n")
cat("===========================\n")

# Load required packages for themes
if (!require("ggplot2", quietly = TRUE)) {
  library(ggplot2)
}

# Publication theme matching your parasiteLoad plots
ferreira_theme <- theme_classic() +
  theme(
    # Text elements
    text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, color = "gray30"),

    # Legend
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.direction = "horizontal",

    # Panel styling
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.3),

    # Strip text for facets
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),

    # Margins (using unit function)
    plot.margin = unit(c(10, 10, 10, 10), "pt")
  )

# Color schemes matching your parasiteLoad plots
ferreira_colors <- list(
  # Analysis types
  analysis = c(
    "Complete Dataset" = "#EEEE00",     # Blue
    "Uninfected Subset" = "#EE7600",    # Purple
    "Infected Subset" = "#54FF9F"       # Orange
  ),

  # Effect significance
  significance = c(
    "Significant" = "#D62828",          # Red
    "Not Significant" = "#6C757D"       # Gray
  ),

  # Sex colors (matching your scheme)
  sex = c(
    "Female" = "#4daf4a",               # Green
    "Male" = "#ff7f00"                  # Orange
  ),

  # Infection status (matching your scheme)
  infection = c(
    "Uninfected" = "#A6CEE3",           # Light blue
    "Infected" = "#FF7094"              # Pink
  )
)

cat("✓ Publication themes and colors configured\n\n")

# ==============================================================================
# 2. EFFECT SIZE PLOTS - MAIN FERREIRA RESULTS
# ==============================================================================

cat("2. CREATING EFFECT SIZE PLOTS\n")
cat("=============================\n")

# Function to create effect size forest plots
create_ferreira_effect_plot <- function(results_list, subtitle = NULL) {

  # Combine all results with analysis labels
  combined_results <- bind_rows(
    results_list$complete %>% mutate(Analysis = "Complete Dataset"),
    results_list$uninfected %>% mutate(Analysis = "Uninfected Subset"),
    results_list$infected %>% mutate(Analysis = "Infected Subset")
  )

  # Clean up effect names for plotting
  combined_results <- combined_results %>%
    filter(!is.na(Estimate)) %>%
    mutate(
      # Clean effect names - concise for axes
      Effect_Clean = case_when(
        str_detect(Interpretation, "Subspecies") ~ "Subspecies distance",
        str_detect(Interpretation, "Hybridization.*hHe-dist") ~ "Hybridization distance",
        str_detect(Interpretation, "Mean hybridization") ~ "Mean hybridization",
        str_detect(Interpretation, "Sex") ~ "Sex difference",
        str_detect(Interpretation, "Interaction") ~ "Subspecies × Hybridization",
        TRUE ~ Interpretation
      ),

      # Analysis factor with proper order
      Analysis = factor(Analysis, levels = c("Complete Dataset", "Uninfected Subset", "Infected Subset")),

      # Significance styling
      Sig_Status = ifelse(Significant, "Significant", "Not Significant"),

      # Effect direction for coloring
      Effect_Direction = case_when(
        Significant & Estimate > 0 ~ "Positive",
        Significant & Estimate < 0 ~ "Negative",
        !Significant ~ "Not Significant",
        TRUE ~ "Not Significant"
      )
    )

  # Create forest plot
  p <- ggplot(combined_results, aes(x = Estimate, y = Effect_Clean)) +

    # Add zero reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +

    # Error bars (confidence intervals)
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Effect_Direction),
                   height = 0.3, linewidth = 1, alpha = 0.8) +

    # Point estimates
    geom_point(aes(color = Effect_Direction, size = Sig_Status), alpha = 0.9) +

    # Facet by analysis type
    facet_wrap(~ Analysis, ncol = 3, scales = "free_x") +

    # Scales
    scale_color_manual(
      values = c(
        "Positive" = "#CD1076",
        "Negative" = "#7FFF00",
        "Not Significant" = "#6C757D"
      ),
      name = "Effect"
    ) +

    scale_size_manual(
      values = c("Significant" = 3, "Not Significant" = 2),
      name = "Significance"
    ) +

    # Clean labels - let the axes speak
    labs(
      x = "Effect size (coefficient ± 95% CI)",
      y = "Predictor variable"
    ) +

    # Apply theme
    ferreira_theme +
    theme(
      strip.background = element_rect(fill = "gray90", color = "black"),
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "gray95"),
      panel.grid.major.y = element_blank(),
      plot.title = element_blank(),  # Remove title
      plot.subtitle = element_blank()  # Remove subtitle
    )

  return(p)
}

# Create main effect size plot
ferreira_results_combined <- list(
  complete = complete_results,
  uninfected = uninfected_results,
  infected = infected_results
)

main_effects_plot <- create_ferreira_effect_plot(ferreira_results_combined)

print(main_effects_plot)

cat("✓ Main effects plot created\n")

# ==============================================================================
# 3. VALIDATION COMPARISON PLOTS
# ==============================================================================

cat("\n3. CREATING VALIDATION COMPARISON PLOTS\n")
cat("=======================================\n")

# Create comparison summary data
create_validation_comparison_data <- function() {

  # parasiteLoad results (from your script output)
  parasiteload_results <- data.frame(
    Analysis = c("Complete Dataset", "Uninfected Subset", "Male-Specific", "Female-Specific"),
    Method = "parasiteLoad",
    P_Value = c(0.01742, 0.5446, 0.03808, 0.1893),
    Significant = c(TRUE, FALSE, TRUE, FALSE),
    Effect_Type = c("Overall Hybrid", "Constitutive Costs", "Sex-Specific", "Sex-Specific"),
    stringsAsFactors = FALSE
  )

  # Ferreira results summary
  ferreira_summary <- data.frame(
    Analysis = c("Complete Dataset", "Uninfected Subset", "Infected Subset"),
    Method = "Ferreira",
    P_Value = c(
      min(complete_results$P_Value[complete_results$Significant], na.rm = TRUE),
      min(uninfected_results$P_Value[uninfected_results$Significant], na.rm = TRUE),
      min(infected_results$P_Value[infected_results$Significant], na.rm = TRUE)
    ),
    Significant = c(
      any(complete_results$Significant, na.rm = TRUE),
      any(uninfected_results$Significant, na.rm = TRUE),
      any(infected_results$Significant, na.rm = TRUE)
    ),
    Effect_Type = c("Overall Hybrid", "Constitutive Costs", "Infection-Specific"),
    stringsAsFactors = FALSE
  )

  # Handle infinite p-values
  ferreira_summary$P_Value[is.infinite(ferreira_summary$P_Value)] <- 1.0
  ferreira_summary$P_Value[is.na(ferreira_summary$P_Value)] <- 1.0

  return(list(parasiteload = parasiteload_results, ferreira = ferreira_summary))
}

# Create validation heatmap
create_validation_heatmap <- function() {

  comparison_data <- create_validation_comparison_data()

  # Combine data for heatmap
  heatmap_data <- bind_rows(
    comparison_data$parasiteload %>% dplyr::select(Analysis, Method, Significant, Effect_Type),
    comparison_data$ferreira %>% dplyr::select(Analysis, Method, Significant, Effect_Type)
  ) %>%
    mutate(
      # Create outcome categories
      Outcome = case_when(
        Analysis == "Complete Dataset" & Significant ~ "Overall Effects Detected",
        Analysis == "Uninfected Subset" & !Significant ~ "No Constitutive Costs",
        Analysis == "Infected Subset" & Significant ~ "Infection-Dependent Effects",
        Analysis == "Male-Specific" & Significant ~ "Male-Specific Effects",
        Analysis == "Female-Specific" & !Significant ~ "Female Protection",
        TRUE ~ "Other"
      ),

      # Validation status
      Validation_Status = case_when(
        Outcome %in% c("Overall Effects Detected", "No Constitutive Costs",
                       "Infection-Dependent Effects", "Male-Specific Effects") ~ "Supports Hypothesis",
        Outcome == "Female Protection" ~ "Supports Hypothesis",
        TRUE ~ "Neutral"
      ),

      # Method factor
      Method = factor(Method, levels = c("parasiteLoad", "Ferreira"))
    )

  # Create heatmap
  p <- ggplot(heatmap_data, aes(x = Method, y = Analysis)) +
    geom_tile(aes(fill = Validation_Status), color = "white", linewidth = 1) +
    geom_text(aes(label = ifelse(Significant, "✓", "○")),
              size = 8, color = "white", fontface = "bold") +

    scale_fill_manual(
      values = c(
        "Supports Hypothesis" = "#2E8B57",    # Green
        "Neutral" = "#6C757D"                 # Gray
      ),
      name = "Validation"
    ) +

    # Clean labels - no redundant titles
    labs(
      x = "Statistical method",
      y = "Analysis type"
    ) +

    ferreira_theme +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  return(p)
}

validation_heatmap <- create_validation_heatmap()
print(validation_heatmap)

cat("✓ Validation comparison heatmap created\n")

# ==============================================================================
# 4. MODEL PERFORMANCE & DIAGNOSTICS
# ==============================================================================

cat("\n4. CREATING MODEL PERFORMANCE PLOTS\n")
cat("===================================\n")

# Model R-squared comparison
create_model_performance_plot <- function() {

  # Extract R-squared values
  performance_data <- data.frame(
    Analysis = c("Complete Dataset", "Uninfected Subset", "Infected Subset"),
    R_squared = c(
      summary(ferreira_complete_model)$r.squared,
      summary(ferreira_uninfected_model)$r.squared,
      summary(ferreira_infected_model)$r.squared
    ),
    N_Pairs = c(
      nrow(complete_pairwise),
      nrow(uninfected_pairwise),
      nrow(infected_pairwise)
    ),
    Significant_Effects = c(
      sum(complete_results$Significant, na.rm = TRUE),
      sum(uninfected_results$Significant, na.rm = TRUE),
      sum(infected_results$Significant, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Analysis = factor(Analysis, levels = c("Complete Dataset", "Uninfected Subset", "Infected Subset")),
      Performance_Category = case_when(
        R_squared > 0.01 ~ "Good Fit",
        R_squared > 0.005 ~ "Moderate Fit",
        TRUE ~ "Weak Fit"
      )
    )

  # Create performance plot
  p <- ggplot(performance_data, aes(x = Analysis, y = R_squared)) +
    geom_col(aes(fill = Performance_Category), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = paste0("R² = ", round(R_squared, 4), "\n",
                                 Significant_Effects, " sig. effects")),
              vjust = -0.5, size = 3.5, fontface = "bold") +

    scale_fill_manual(
      values = c(
        "Good Fit" = "#FF4500",
        "Moderate Fit" = "#FFA500",
        "Weak Fit" = "yellow"
      ),
      name = "Model Fit"
    ) +

    scale_y_continuous(
      labels = scales::percent_format(accuracy = 0.1),
      expand = expansion(mult = c(0, 0.15))
    ) +

    labs(
      title = "",
      subtitle = "",
      x = "Analysis Type",
      y = "Variance Explained (R²)"
    ) +

    ferreira_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

model_performance_plot <- create_model_performance_plot()
print(model_performance_plot)

cat("✓ Model performance plot created\n")

# ==============================================================================
# 5. HYBRID INDEX EFFECT VISUALIZATION
# ==============================================================================

cat("\n5. CREATING HYBRID INDEX EFFECT PLOTS\n")
cat("=====================================\n")

# Create effect plots showing hybrid index relationships
create_hybrid_effect_plots <- function() {

  # Extract significant effects from infected subset (strongest signals)
  significant_effects <- infected_results %>%
    filter(Significant) %>%
    arrange(P_Value)

  if (nrow(significant_effects) == 0) {
    cat("No significant effects found for visualization\n")
    return(NULL)
  }

  # Create prediction data for visualization
  pred_data <- expand.grid(
    subspecies_genetic_distance_scaled = seq(0, 1, length.out = 50),
    hHe_distance_scaled = c(0.2, 0.5, 0.8),  # Different hybridization levels
    hHe_mean_scaled = 0.5,  # Average
    sex_distance = 0  # Same sex pairs
  )

  # Generate predictions
  pred_data$predicted_similarity <- predict(ferreira_infected_model, pred_data)



  # Create effect plot
  p <- ggplot(pred_data, aes(x = subspecies_genetic_distance_scaled,
                             y = predicted_similarity,
                             color = factor(hHe_distance_scaled))) +
    geom_line(size = 1.2, alpha = 0.8) +
    scale_color_manual(
      values = c("0.2" = "#FF4500", "0.5" = "#FFA500", "0.8" = "yellow"),
      labels = c("Low (0.2)", "Medium (0.5)", "High (0.8)"),
      name = "Hybridization\nDistance"
    ) +

    scale_x_continuous(
      breaks = seq(0, 1, 0.25),
      labels = c("0", "0.25", "0.5", "0.75", "1")
    ) +

    labs(
      title = "",
      subtitle = "",
      x = "Subspecies Genetic Distance (HI difference)",
      y = "Predicted Health Similarity"
    ) +

    ferreira_theme

  return(p)
}

hybrid_effect_plot <- create_hybrid_effect_plots()
if (!is.null(hybrid_effect_plot)) {
  print(hybrid_effect_plot)
  cat("✓ Hybrid effect visualization created\n")
}

# ==============================================================================
# 6. COMBINED FIGURE PANEL FOR MANUSCRIPT
# ==============================================================================

cat("\n6. CREATING COMBINED FIGURE PANEL\n")
cat("=================================\n")

# Create combined figure for Figure 2
create_figure2_panel <- function() {

  # Panel A: Main effect sizes
  panel_a <- main_effects_plot +
    labs(tag = "A") +
    theme(plot.tag = element_text(size = 18, face = "bold"))

  # Panel B: Validation comparison
  panel_b <- validation_heatmap +
    labs(tag = "B") +
    theme(
      plot.tag = element_text(size = 18, face = "bold"),
      legend.position = "right"
    )

  # Panel C: Model performance
  panel_c <- model_performance_plot +
    labs(tag = "C") +
    theme(plot.tag = element_text(size = 18, face = "bold"))

  # Combine panels using patchwork
  if (!is.null(hybrid_effect_plot)) {
    # Include hybrid effect plot as Panel D
    panel_d <- hybrid_effect_plot +
      labs(tag = "D") +
      theme(plot.tag = element_text(size = 18, face = "bold"))

    combined_figure <- (panel_a | panel_b) / (panel_c | panel_d)
  } else {
    # Three-panel layout
    combined_figure <- panel_a / (panel_b | panel_c)
  }

  # Add overall title
  final_figure <- combined_figure +
    plot_annotation(
      title = "",
      subtitle = "",
      theme = theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 16, hjust = 0.5, color = "gray30")
      )
    )

  return(final_figure)
}

figure2_panel <- create_figure2_panel()
print(figure2_panel)

cat("✓ Combined Figure 2 panel created\n")

# ==============================================================================
# 7. SAVE ALL FIGURES
# ==============================================================================

cat("\n7. SAVING PUBLICATION-READY FIGURES\n")
cat("===================================\n")

# Individual plots
save_plot_all_formats_wide(main_effects_plot, "Ferreira_main_effects")
save_plot_all_formats(validation_heatmap, "Ferreira_validation_comparison")
save_plot_all_formats(model_performance_plot, "Ferreira_model_performance")

if (!is.null(hybrid_effect_plot)) {
  save_plot_all_formats(hybrid_effect_plot, "Ferreira_hybrid_effects")
}

# Combined figure panel
save_plot_all_formats(figure2_panel, "Figure2_Ferreira_validation_panel")

# High-resolution version for publication
ggsave(
  filename = file.path("results", "figures", "Figure2_Ferreira_validation_panel", "Figure2_publication_ready.pdf"),
  plot = figure2_panel,
  width = 16, height = 12, dpi = 300, device = "pdf"
)

ggsave(
  filename = file.path("results", "figures", "Figure2_Ferreira_validation_panel", "Figure2_publication_ready.png"),
  plot = figure2_panel,
  width = 16, height = 12, dpi = 600, bg = "white"
)

cat("✓ All figures saved in multiple formats\n")
cat("✓ Publication-ready versions created\n\n")

cat("=== FERREIRA FIGURE GENERATION COMPLETE ===\n")
cat("Beautiful publication-ready figures created!\n")
cat("Ready for Figure 2 manuscript integration! 🎨\n")

