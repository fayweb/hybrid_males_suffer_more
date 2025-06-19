

# DISTANCE-BASED KEY FINDINGS FIGURE
# Purpose: Highlight the most important sex-specific predictability patterns
# FERREIRA-STYLE CLEAR DISTANCE PLOTS
# Purpose: Create clean, interpretable forest plots like Ferreira et al.

# ==============================================================================
# FERREIRA-STYLE FOREST PLOT
# ==============================================================================

# Prepare your data in Ferreira format
ferreira_style_data <- bind_rows(
  # Complete dataset
  complete_results %>% mutate(Analysis = "Complete Dataset"),
  uninfected_results %>% mutate(Analysis = "Uninfected Subset"),
  infected_results %>% mutate(Analysis = "Infected Subset"),
  male_infected_results %>% mutate(Analysis = "Males Infected"),
  female_infected_results %>% mutate(Analysis = "Females Infected")
) %>%
  filter(!is.na(Estimate)) %>%
  filter(Effect %in% c("subspecies_distance", "hHe_distance", "hHe_mean")) %>%
  mutate(
    # Clean effect names
    Effect_Clean = case_when(
      Effect == "subspecies_distance" ~ "Subspecies genetic distance",
      Effect == "hHe_distance" ~ "Hybridization distance",
      Effect == "hHe_mean" ~ "Mean hybridization level",
      TRUE ~ Effect
    ),

    # Analysis colors (like Ferreira's microbiome components)
    Analysis_Color = case_when(
      Analysis == "Complete Dataset" ~ "#1f77b4",      # Blue
      Analysis == "Uninfected Subset" ~ "#ff7f0e",     # Orange
      Analysis == "Infected Subset" ~ "#2ca02c",       # Green
      Analysis == "Males Infected" ~ "#d62728",        # Red
      Analysis == "Females Infected" ~ "#9467bd",      # Purple
      TRUE ~ "#7f7f7f"
    ),

    # Calculate confidence intervals
    CI_Lower = Estimate - 1.96 * Std_Error,
    CI_Upper = Estimate + 1.96 * Std_Error,

    # Factor levels for proper ordering
    Effect_Clean = factor(Effect_Clean,
                          levels = c("Subspecies genetic distance",
                                     "Hybridization distance",
                                     "Mean hybridization level")),
    Analysis = factor(Analysis,
                      levels = c("Complete Dataset", "Uninfected Subset",
                                 "Infected Subset", "Males Infected", "Females Infected"))
  )

# Create Ferreira-style forest plot
ferreira_plot <- ggplot(ferreira_style_data, aes(x = Estimate, y = Effect_Clean)) +
  # Reference line at zero
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.5) +

  # Error bars
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper, color = Analysis),
                 height = 0.15, size = 0.8, alpha = 0.7,
                 position = position_dodge(width = 0.6)) +

  # Point estimates
  geom_point(aes(color = Analysis, size = Significant),
             alpha = 0.8, position = position_dodge(width = 0.6)) +

  # Clean color scheme
  scale_color_manual(values = c("Complete Dataset" = "#1f77b4",
                                "Uninfected Subset" = "#00FFFF",
                                "Infected Subset" = "#FF7094",
                                "Males Infected" = "#ff7f00",
                                "Females Infected" = "#4daf4a")) +

  scale_size_manual(values = c("TRUE" = 2.5, "FALSE" = 1.8), guide = "none") +

  # Clean labels
  labs(
    x = "Effect size estimate (± 95% CI)",
    y = "Genetic predictor variable",
    color = "Analysis type",
    title = ""
  ) +

  # Ferreira-style theme
  theme_classic() +
  theme(
    # Text sizes
    text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),

    # Legend
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),

    # Grid
    panel.grid.major.x = element_line(color = "gray95", size = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    # Clean appearance
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(ferreira_plot)


cat("\n=== DISTANCE-BASED KEY FINDINGS FIGURE ===\n")

save_plot_all_formats_extra_wide(plot_object = ferreira_plot, plot_name = "Distance_based_validation_sex_specific_hybrid")

# Load your results
load("ferreira_sex_stratified_results.RData")

# ==============================================================================
# PANEL A: SEX-SPECIFIC HYBRIDIZATION EFFECTS (THE MAIN STORY)
# ==============================================================================

# Extract key hybridization distance effects
key_effects_data <- bind_rows(
  # Overall infected for reference
  infected_results %>%
    filter(Effect == "hHe_distance") %>%
    mutate(Group = "Overall Infected", Sex = "Combined"),

  # Sex-specific effects
  male_infected_results %>%
    filter(Effect == "hHe_distance") %>%
    mutate(Group = "Males", Sex = "Male"),

  female_infected_results %>%
    filter(Effect == "hHe_distance") %>%
    mutate(Group = "Females", Sex = "Female")
) %>%
  mutate(
    Group = factor(Group, levels = c("Overall Infected", "Males", "Females")),
    CI_Lower = Estimate - 1.96 * Std_Error,
    CI_Upper = Estimate + 1.96 * Std_Error,
    Significance_Label = ifelse(Significant,
                                paste0("β = ", sprintf("%.3f", Estimate),
                                       "\np < 0.001"),
                                paste0("β = ", sprintf("%.3f", Estimate),
                                       "\np = ", sprintf("%.3f", P_Value)))
  )

# Panel A: Main effects
panel_a <- ggplot(key_effects_data, aes(x = Group, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper, color = Group),
                width = 0.2, size = 1.2, alpha = 0.8) +
  geom_point(aes(color = Group, size = Significant), alpha = 0.9) +
  geom_text(aes(label = Significance_Label, color = Group),
            vjust = -0.5, hjust = 0.5, size = 3.5, fontface = "bold") +
  scale_color_manual(values = c("Overall Infected" = "#FF7094",
                                "Males" = "#ff7f0e",
                                "Females" = "#4daf4a")) +
  scale_size_manual(values = c("TRUE" = 4, "FALSE" = 3), guide = "none") +
  labs(
    title = "",
    subtitle = "Hybridization distance effects on predicted weight loss similarity",
    x = "Analysis group",
    y = "Effect size (β ± 95% CI)",
    color = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray95", size = 0.3)
  )

panel_a

save_plot_all_formats(plot_object = panel_a, plot_name = "main_effects")
# ==============================================================================
# PANEL B: INTERPRETATION DIAGRAM
# ==============================================================================

# Create interpretation data
interpretation_data <- data.frame(
  Sex = c("Females", "Males"),
  Effect_Size = c(0.146, -0.001),  # Your actual values
  P_Value = c(0.000, 0.962),       # Your actual p-values
  Interpretation = c("Structured\nVulnerability", "Chaotic\nVulnerability"),
  Description = c("Genetics predicts\nhealth outcomes", "Genetics doesn't predict\nhealth outcomes"),
  Significance = c(TRUE, FALSE)
)

# Panel B: Interpretation
panel_b <- ggplot(interpretation_data, aes(x = Sex, y = Effect_Size)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_col(aes(fill = Sex, alpha = Significance), width = 0.6) +
  geom_text(aes(label = Interpretation),
            y = 0.08, vjust = 0, size = 4, fontface = "bold") +
  geom_text(aes(label = Description),
            y = 0.05, vjust = 0, size = 3.5, color = "gray30") +
  scale_fill_manual(values = c("Females" = "#1f77b4", "Males" = "#ff7f0e")) +
  scale_alpha_manual(values = c("TRUE" = 0.8, "FALSE" = 0.4)) +
  labs(
    title = "Biological interpretation",
    subtitle = "Sex-specific patterns of health predictability",
    x = "Sex",
    y = "Hybridization distance effect"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "none"
  )

panel_b
# ==============================================================================
# PANEL C: MODEL PERFORMANCE COMPARISON
# ==============================================================================

# Extract R-squared values for comparison
performance_data <- data.frame(
  Analysis = c("Complete Dataset", "Uninfected Subset", "Infected Subset",
               "Males Infected", "Females Infected"),
  R_squared = c(
    unique(complete_results$R_squared)[1],
    unique(uninfected_results$R_squared)[1],
    unique(infected_results$R_squared)[1],
    unique(male_infected_results$R_squared)[1],
    unique(female_infected_results$R_squared)[1]
  ),
  Significant_Effects = c(
    sum(complete_results$Significant, na.rm = TRUE),
    sum(uninfected_results$Significant, na.rm = TRUE),
    sum(infected_results$Significant, na.rm = TRUE),
    sum(male_infected_results$Significant, na.rm = TRUE),
    sum(female_infected_results$Significant, na.rm = TRUE)
  ),
  Group_Type = c("Overall", "Overall", "Overall", "Sex-specific", "Sex-specific")
) %>%
  mutate(
    Analysis = factor(Analysis, levels = c("Complete Dataset", "Uninfected Subset",
                                           "Infected Subset", "Males Infected", "Females Infected")),
    Highlight = Analysis %in% c("Males Infected", "Females Infected")
  )

# Panel C: Model performance
panel_c <- ggplot(performance_data, aes(x = Analysis, y = R_squared)) +
  geom_col(aes(fill = Highlight, alpha = Group_Type), width = 0.7) +
  geom_text(aes(label = paste0("R² = ", sprintf("%.4f", R_squared),
                               "\n", Significant_Effects, " effects")),
            vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "#d62728")) +
  scale_alpha_manual(values = c("Overall" = 0.6, "Sex-specific" = 0.9)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "Model performance across analyses",
    subtitle = "Sex-specific models show distinct patterns",
    x = "Analysis type",
    y = "Variance explained (R²)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

panel_c
# ==============================================================================

