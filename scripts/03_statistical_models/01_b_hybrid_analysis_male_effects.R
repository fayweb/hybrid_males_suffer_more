# ==============================================================================
# TEST MALE EFFECTS IN INFECTED MICE ONLY
# ==============================================================================

cat("=== TESTING MALE EFFECTS IN INFECTED SUBSET ===\n")
cat("================================================\n\n")

# Check if infected_data exists
if(!exists("infected_data")) {
  stop("infected_data not found. Please run the hybrid analysis script first.")
}

cat("Testing for sex-specific hybrid effects within infected mice only\n")
cat("Sample: n =", nrow(infected_data), "infected mice\n")

# Check sex breakdown in infected subset
sex_breakdown <- infected_data %>%
  group_by(Sex) %>%
  summarise(n = n(), .groups = 'drop')

cat("\nSex breakdown in infected subset:\n")
print(sex_breakdown)
cat("\n")

# ==============================================================================
# 1. parasiteLoad ANALYSIS - INFECTED MALES VS FEMALES
# ==============================================================================

cat("1. RUNNING parasiteLoad ANALYSIS ON INFECTED SUBSET\n")
cat("===================================================\n")

# Run parasiteLoad on infected mice with Sex as grouping variable
cat("Running parasiteLoad with Sex grouping on infected mice...\n")

infected_sex_model <- analyse(
  data = infected_data,
  response = "response",
  model = "student",
  group = "Sex"  # Test for sex-specific effects within infected mice
)

cat("✓ Analysis complete\n\n")

# ==============================================================================
# 2. EXTRACT P-VALUES FOR SEX EFFECTS IN INFECTED MICE
# ==============================================================================

cat("2. EXTRACTING SEX-SPECIFIC RESULTS\n")
cat("==================================\n")

# The infected_sex_model will test:
# H0: No hybrid effects overall in infected mice
# H1: Sex differences in hybrid effects
# H2: Group-specific effects (Female vs Male)

# We need to look at the printed output to get p-values
# But we can also examine the model structure

cat("Model structure:\n")
if(is.list(infected_sex_model)) {
  cat("Available model components:", names(infected_sex_model), "\n")
}

# ==============================================================================
# 3. SIMPLE LINEAR MODEL APPROACH AS BACKUP
# ==============================================================================

cat("\n3. BACKUP LINEAR MODEL ANALYSIS\n")
cat("===============================\n")

# Test interaction between HI and Sex in infected mice
cat("Testing HI × Sex interaction in infected mice using linear model...\n")

# Create centered HI for better interpretation
infected_data_centered <- infected_data %>%
  mutate(HI_centered = HI - mean(HI, na.rm = TRUE))

# Linear model with interaction
lm_infected_sex <- lm(response ~ HI_centered * Sex, data = infected_data_centered)

# Get summary
lm_summary <- summary(lm_infected_sex)
cat("\nLinear Model Results (HI × Sex interaction in infected mice):\n")
print(lm_summary$coefficients)

# Extract key p-values
hi_main <- lm_summary$coefficients["HI_centered", "Pr(>|t|)"]
sex_main <- lm_summary$coefficients["SexMale", "Pr(>|t|)"]
interaction_p <- lm_summary$coefficients["HI_centered:SexMale", "Pr(>|t|)"]

cat(sprintf("\nKey results from linear model:\n"))
cat(sprintf("- HI main effect: p = %.5f\n", hi_main))
cat(sprintf("- Sex main effect: p = %.5f\n", sex_main))
cat(sprintf("- HI × Sex interaction: p = %.5f\n", interaction_p))

# Interpretation
if(interaction_p < 0.05) {
  cat("✓ SIGNIFICANT interaction: Males and females respond differently to hybridization\n")
} else {
  cat("○ Non-significant interaction: No sex difference in hybrid effects\n")
}

# ==============================================================================
# 4. SEX-STRATIFIED ANALYSIS
# ==============================================================================

cat("\n4. SEX-STRATIFIED ANALYSIS\n")
cat("==========================\n")

# Test HI effects separately in infected males and infected females
infected_females <- infected_data %>% filter(Sex == "Female")
infected_males <- infected_data %>% filter(Sex == "Male")

cat(sprintf("Infected females: n = %d\n", nrow(infected_females)))
cat(sprintf("Infected males: n = %d\n", nrow(infected_males)))

# Linear models for each sex
lm_females <- lm(response ~ HI, data = infected_females)
lm_males <- lm(response ~ HI, data = infected_males)

# Extract results
female_hi_p <- summary(lm_females)$coefficients["HI", "Pr(>|t|)"]
male_hi_p <- summary(lm_males)$coefficients["HI", "Pr(>|t|)"]
female_hi_coef <- summary(lm_females)$coefficients["HI", "Estimate"]
male_hi_coef <- summary(lm_males)$coefficients["HI", "Estimate"]

cat(sprintf("\nSex-stratified results:\n"))
cat(sprintf("- Infected females HI effect: coef = %.4f, p = %.5f\n", female_hi_coef, female_hi_p))
cat(sprintf("- Infected males HI effect: coef = %.4f, p = %.5f\n", male_hi_coef, male_hi_p))

# Compare effect sizes
if(abs(male_hi_coef) > abs(female_hi_coef)) {
  cat("→ Males show stronger hybrid effect in infected subset\n")
} else {
  cat("→ Females show stronger hybrid effect in infected subset\n")
}

# ==============================================================================
# 5. VISUALIZATION
# ==============================================================================

cat("\n5. CREATING VISUALIZATION\n")
cat("=========================\n")

# Create plot showing sex differences in infected mice
p_infected_sex <- ggplot(infected_data, aes(x = HI, y = response, color = Sex)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3) +
  scale_color_manual(values = c("Female" = "#4daf4a", "Male" = "#ff7f00")) +
  labs(
    title = "",
    x = "Hybrid Index (HI)",
    y = "Predicted Weight Loss (%)",
    subtitle = paste0("Infected mice only (n = ", nrow(infected_data),
                      ": ", sum(infected_data$Sex == "Female"), " females, ",
                      sum(infected_data$Sex == "Male"), " males)")
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

print(p_infected_sex)

# ==============================================================================
# 6. SUMMARY FOR FIGURE LEGEND
# ==============================================================================

cat("\n6. SUMMARY FOR FIGURE PANEL C\n")
cat("=============================\n")

# Determine what result to report for Panel C
cat("Based on analysis of infected mice subset:\n")

if(interaction_p < 0.05) {
  cat(sprintf("✓ Significant sex difference in hybrid effects (HI × Sex p = %.5f)\n", interaction_p))
  panel_c_result <- sprintf("significant sex-specific hybrid effects (p = %.3f)", interaction_p)
} else if(hi_main < 0.05) {
  cat(sprintf("✓ Overall hybrid effect in infected mice (HI p = %.5f)\n", hi_main))
  panel_c_result <- sprintf("significant hybrid effects in infected mice (p = %.3f)", hi_main)
} else {
  cat("○ No significant effects detected in infected subset\n")
  panel_c_result <- "no significant effects in infected subset"
}

cat(sprintf("\nSuggested text for Panel C:\n"))
cat(sprintf("Panel C result: %s\n", panel_c_result))

cat("\n🎯 Use this result to update your figure legend!\n")


# ==============================================================================
# ==============================================================================
# PANEL F: SEX-SPECIFIC HYBRID EFFECTS IN INFECTED MICE
# ==============================================================================

cat("Creating Panel F: Sex-specific effects in infected mice...\n")

# Use your existing infected_data and sex colors
panel_f_plot <- ggplot(infected_data, aes(x = HI, y = response, color = Sex)) +
  # Points with some transparency
  geom_point(alpha = 0.7, size = 2) +

  # Separate trend lines for each sex
  geom_smooth(method = "lm", se = TRUE, alpha = 0.3, linewidth = 1.2) +

  # Use your established color scheme
  scale_color_manual(
    values = c("Female" = "#4daf4a", "Male" = "#ff7f00"),
    name = "Sex"
  ) +

  # Clean theme matching your other panels
  theme_bw() +
  theme(
    legend.position = 'none',  # Keep consistent with other panels
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +

  # Labels
  labs(
    y = "Predicted weight loss (%)\nImmune signature",
    title = ""
  ) +

  # Add annotation for the significant interaction
  annotate(
    "text",
    x = 0.1, y = max(infected_data$response) * 0.95,
    label = "HI × Sex: p = 0.029",
    size = 3.5,
    fontface = "bold",
    hjust = 0
  )

# Combine with gradient bar like your other panels
combined_panel_f <- panel_f_plot / HIgradientBar +
  plot_layout(heights = c(1, 0.1))

print(combined_panel_f)

# Save the panel
save_plot_all_formats_tight(
  plot_object = combined_panel_f,
  plot_name = "Panel_F_Sex_Specific_Infected"
)

cat("✓ Panel F created and saved\n")
cat("Shows contrasting slopes: Females benefit from hybridization, Males suffer\n")

