# ==============================================================================
# DISTRIBUTION ANALYSIS: Understanding our data before modeling
# ==============================================================================
#
# Purpose: Examine distributions of predicted weight loss and other key variables
# to make informed decisions about statistical models for parasiteLoad analysis
#
# Author: Fay Webster
# Date: June 2025
# ==============================================================================

# Load required packages
library(fitdistrplus)  # For distribution fitting
library(moments)       # For skewness and kurtosis
library(VineCopula)    # For distribution testing
library(gridExtra)     # For plot arrangement

if (!exists("field_mice") && !exists("MASTER_SCRIPT_RUNNING")) {
  stop("Please run master script first to load field_mice data")
}

cat("=== DISTRIBUTION ANALYSIS ===\n")
cat("Understanding our data before statistical modeling\n\n")

# ==============================================================================
# 1. DATA OVERVIEW
# ==============================================================================

cat("1. DATA OVERVIEW\n")
cat("================\n")

# Check data structure
cat("Dataset dimensions:", dim(field_mice), "\n")
cat("Column names:\n")
print(names(field_mice))

# Key variables for analysis
key_vars <- c("HI", "Sex", "predicted_weight_loss", "MC.Eimeria", "infection_status")

cat("\nKey variables availability:\n")
for (var in key_vars) {
  if (var %in% names(field_mice)) {
    missing_n <- sum(is.na(field_mice[[var]]))
    missing_pct <- round(missing_n / nrow(field_mice) * 100, 1)
    cat(sprintf("✓ %s: %d missing (%.1f%%)\n", var, missing_n, missing_pct))
  } else {
    cat(sprintf("✗ %s: NOT FOUND\n", var))
  }
}

# ==============================================================================
# 2. PREDICTED WEIGHT LOSS DISTRIBUTION
# ==============================================================================

cat("\n\n2. PREDICTED WEIGHT LOSS DISTRIBUTION\n")
cat("=====================================\n")

# Check if predicted_weight_loss exists and clean it
if ("predicted_weight_loss" %in% names(field_mice)) {

  # Remove missing values for analysis
  wl_data <- field_mice$predicted_weight_loss[!is.na(field_mice$predicted_weight_loss)]

  cat("Sample size with predicted weight loss:", length(wl_data), "\n")
  cat("Range:", round(range(wl_data), 3), "\n")
  cat("Mean ± SD:", round(mean(wl_data), 3), "±", round(sd(wl_data), 3), "\n")
  cat("Median (IQR):", round(median(wl_data), 3),
      "(", round(quantile(wl_data, 0.25), 3), "-", round(quantile(wl_data, 0.75), 3), ")\n")

  # Check for negative values (shouldn't happen with weight loss)
  negative_count <- sum(wl_data < 0)
  cat("Negative values:", negative_count, "\n")

  # Check for extreme values
  q99 <- quantile(wl_data, 0.99)
  extreme_count <- sum(wl_data > q99)
  cat("Values above 99th percentile (", round(q99, 3), "):", extreme_count, "\n")

  # Basic distribution shape
  cat("Skewness:", round(skewness(wl_data), 3), "\n")
  cat("Kurtosis:", round(kurtosis(wl_data), 3), "\n")

} else {
  stop("predicted_weight_loss column not found!")
}

# ==============================================================================
# 3. VISUAL EXAMINATION OF DISTRIBUTIONS
# ==============================================================================

cat("\n\n3. CREATING DISTRIBUTION PLOTS\n")
cat("==============================\n")

# Create clean dataset for plotting
plot_data <- field_mice %>%
  filter(!is.na(predicted_weight_loss), !is.na(HI)) %>%
  mutate(
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "Male")),
    infection_status = ifelse(is.na(infection_status), "Unknown",
                              ifelse(infection_status, "Infected", "Uninfected"))
  )

cat("Clean dataset for plotting: n =", nrow(plot_data), "\n")

# Plot 1: Overall distribution
# Simple histogram that should work
p1 <- ggplot(plot_data, aes(x = predicted_weight_loss)) +
  geom_histogram(bins = 30, fill = "skyblue", alpha = 0.7, color = "white") +
  labs(title = "A)",
       x = "Predicted Weight Loss (%)",
       y = "Count") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))

print(p1)

save_plot_all_formats(plot_object = p1, plot_name = "Overall_distribution of predicted weight loss")

# Plot 2: By sex
p2 <- ggplot(plot_data, aes(x = predicted_weight_loss, fill = Sex)) +
  geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
  facet_wrap(~Sex, ncol = 1) +
  scale_fill_manual(values = c("Female" = "#E69F00", "Male" = "#0072B2")) +
  labs(title = "B)",
       x = "Predicted Weight Loss (%)",
       y = "Count") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

print(p2)

save_plot_all_formats(plot_object = p2, plot_name = "Distribution_weight_loss_sex")

# Plot 3: By infection status
p3 <- ggplot(plot_data, aes(x = predicted_weight_loss, fill = infection_status)) +
  geom_histogram(bins = 25, alpha = 0.7, position = "identity") +
  facet_wrap(~infection_status, ncol = 1) +
  scale_fill_manual(values = c("Uninfected" = "#2E8B57", "Infected" = "#DC143C", "Unknown" = "gray")) +
  labs(title = "C)",
       x = "Predicted Weight Loss (%)",
       y = "Count") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "none")

print(p3)

save_plot_all_formats(plot_object = p3, plot_name = "Distribution_infection_status")

# Plot 4: Q-Q plot for normality
p4 <- ggplot(plot_data, aes(sample = predicted_weight_loss)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "D)",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))

print(p4)

save_plot_all_formats(plot_object = p4, plot_name = "QQ_plot_distribution_weight_loss")

# Combine plots
distribution_grid <- grid.arrange(p1, p2, p3, p4, ncol = 2)

print(distribution_grid)
save_plot_all_formats_panel(plot_object = distribution_grid, plot_name = "Distribution_grid")
# ==============================================================================
# 4. FORMAL DISTRIBUTION TESTING
# ==============================================================================

cat("\n\n4. FORMAL DISTRIBUTION TESTING\n")
cat("==============================\n")

# Test different distributions
wl_clean <- plot_data$predicted_weight_loss

cat("Testing various distributions...\n\n")

# Normal distribution
norm_fit <- fitdist(wl_clean, "norm")
cat("NORMAL DISTRIBUTION:\n")
cat("Parameters: mean =", round(norm_fit$estimate[1], 3),
    ", sd =", round(norm_fit$estimate[2], 3), "\n")
cat("AIC:", round(norm_fit$aic, 2), "\n")
cat("BIC:", round(norm_fit$bic, 2), "\n\n")

# Gamma distribution (for positive continuous data)
if (all(wl_clean > 0)) {
  gamma_fit <- fitdist(wl_clean, "gamma")
  cat("GAMMA DISTRIBUTION:\n")
  cat("Parameters: shape =", round(gamma_fit$estimate[1], 3),
      ", rate =", round(gamma_fit$estimate[2], 3), "\n")
  cat("AIC:", round(gamma_fit$aic, 2), "\n")
  cat("BIC:", round(gamma_fit$bic, 2), "\n\n")
}

# Log-normal distribution (if all positive)
if (all(wl_clean > 0)) {
  lnorm_fit <- fitdist(wl_clean, "lnorm")
  cat("LOG-NORMAL DISTRIBUTION:\n")
  cat("Parameters: meanlog =", round(lnorm_fit$estimate[1], 3),
      ", sdlog =", round(lnorm_fit$estimate[2], 3), "\n")
  cat("AIC:", round(lnorm_fit$aic, 2), "\n")
  cat("BIC:", round(lnorm_fit$bic, 2), "\n\n")
}

# Weibull distribution
if (all(wl_clean > 0)) {
  weibull_fit <- fitdist(wl_clean, "weibull")
  cat("WEIBULL DISTRIBUTION:\n")
  cat("Parameters: shape =", round(weibull_fit$estimate[1], 3),
      ", scale =", round(weibull_fit$estimate[2], 3), "\n")
  cat("AIC:", round(weibull_fit$aic, 2), "\n")
  cat("BIC:", round(weibull_fit$bic, 2), "\n\n")
}

# Shapiro-Wilk test for normality (if n < 5000)
if (length(wl_clean) < 5000) {
  shapiro_test <- shapiro.test(wl_clean)
  cat("SHAPIRO-WILK NORMALITY TEST:\n")
  cat("W =", round(shapiro_test$statistic, 4), "\n")
  cat("p-value =", format(shapiro_test$p.value, scientific = TRUE), "\n")
  if (shapiro_test$p.value < 0.05) {
    cat("Result: Data significantly deviates from normal distribution\n\n")
  } else {
    cat("Result: Data consistent with normal distribution\n\n")
  }
} else {
  cat("Sample size too large for Shapiro-Wilk test\n\n")
}

# ==============================================================================
# 5. HYBRID INDEX DISTRIBUTION
# ==============================================================================

cat("5. HYBRID INDEX DISTRIBUTION\n")
cat("============================\n")

hi_clean <- plot_data$HI
cat("Sample size:", length(hi_clean), "\n")
cat("Range:", round(range(hi_clean), 3), "\n")
cat("Mean ± SD:", round(mean(hi_clean), 3), "±", round(sd(hi_clean), 3), "\n")
cat("Median (IQR):", round(median(hi_clean), 3),
    "(", round(quantile(hi_clean, 0.25), 3), "-", round(quantile(hi_clean, 0.75), 3), ")\n")

# Plot hybrid index distribution
hi_plot <- ggplot(plot_data, aes(x = HI)) +
  geom_histogram(bins = 30, fill = "purple", alpha = 0.7, color = "white") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Hybrid Index Distribution",
       x = "Hybrid Index (0 = M.m.domesticus, 1 = M.m.musculus)",
       y = "Count",
       caption = "Dashed line: hybrid center (HI = 0.5)") +
  theme_classic() +
  theme(plot.title = element_text(size = 12, face = "bold"))

print(hi_plot)

# ==============================================================================
# 6. RECOMMENDATIONS FOR PARASITELOAD MODELS
# ==============================================================================

cat("\n\n6. MODEL RECOMMENDATIONS\n")
cat("========================\n")

cat("Based on distribution analysis:\n\n")

# Decision tree for model choice
if (all(wl_clean > 0)) {
  cat("✓ All predicted weight loss values are positive\n")

  # Check if data looks normal
  if (exists("shapiro_test") && shapiro_test$p.value > 0.05) {
    cat("✓ Data appears approximately normal\n")
    cat("RECOMMENDATION: Use model = 'student' in parasiteLoad\n")
    cat("(Student's t-distribution handles slight deviations from normality)\n\n")
  } else {
    cat("✗ Data deviates significantly from normal distribution\n")

    # Check skewness
    if (skewness(wl_clean) > 1) {
      cat("✓ Data is right-skewed\n")
      cat("RECOMMENDATION: Consider log-transformation + 'student' model\n")
      cat("OR use 'negbin' if treating as count-like data\n\n")
    } else {
      cat("RECOMMENDATION: Try 'student' model first, check residuals\n\n")
    }
  }
} else {
  cat("✗ Some predicted weight loss values are negative or zero\n")
  cat("RECOMMENDATION: Investigate data issues before modeling\n\n")
}

# Additional considerations
cat("ADDITIONAL CONSIDERATIONS:\n")
cat("- Sample size: n =", length(wl_clean), "\n")
if (length(wl_clean) > 50) {
  cat("✓ Sample size adequate for complex models\n")
} else {
  cat("! Small sample size - use simpler models\n")
}

cat("- Range of values:", round(range(wl_clean), 3), "\n")
cat("- Variance:", round(var(wl_clean), 3), "\n")

if (var(wl_clean) > mean(wl_clean)) {
  cat("! Variance > mean suggests overdispersion\n")
  cat("  Consider negative binomial if treating as count data\n")
}

cat("\nREADY TO PROCEED WITH PARASITELOAD ANALYSIS!\n")

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================

# Save the clean dataset for subsequent analysis
analysis_ready_data <- plot_data %>%
  filter(!is.na(predicted_weight_loss), !is.na(HI)) %>%
  dplyr::select(all_of(c("HI", "Sex", "predicted_weight_loss", "infection_status")),
         everything())

cat("\nSaving analysis-ready dataset with", nrow(analysis_ready_data), "complete cases\n")

# Save plots
ggsave("results/figures/Distribution_Analysis.pdf", distribution_grid,
       width = 12, height = 10, dpi = 300)
ggsave("results/figures/Hybrid_Index_Distribution.pdf", hi_plot,
       width = 8, height = 6, dpi = 300)

cat("Distribution analysis complete!\n")
cat("Files saved:\n")
cat("- results/figures/Distribution_Analysis.pdf\n")
cat("- results/figures/Hybrid_Index_Distribution.pdf\n")

# ==============================================================================
# END OF DISTRIBUTION ANALYSIS
# ==============================================================================

