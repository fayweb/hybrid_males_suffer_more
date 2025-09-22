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

# ==============================================================================
# 1. DATA PREPARATION
# ==============================================================================

cat("1. PREPARING COMPLETE DATASET\n")
cat("=============================\n")

# Filter to available genes
available_genes <- immune_genes[immune_genes %in% names(field_mice)]
cat("Available immune genes:", length(available_genes), "\n")

# Create complete dataset with all interactions
field_mice <- field_mice %>%
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
  ) %>%
  drop_na(HI)




# Create complete dataset with all interactions
immune_data <- field_mice %>%
  dplyr::select(Mouse_ID, HI, Sex, infection_status, predicted_weight_loss, species_Eimeria,
                all_of(available_genes), infection) %>%
  filter(!is.na(HI) & !is.na(Sex) & !is.na(infection_status) &
           !is.na(predicted_weight_loss))

immune_data$species_Eimeria <- factor(
  immune_data$species_Eimeria,
  levels = c("Uninfected", "E. ferrisi", "E. falciformis")
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
# 3. TEST MULTIPLE MODEL STRUCTURES FOR PREDICTED WEIGHT LOSS
# ==============================================================================

cat("3. TESTING MODEL STRUCTURES FOR PREDICTED WEIGHT LOSS\n")
cat("====================================================\n")

# Model 1: Simple HI × He
model1 <- lm(predicted_weight_loss ~ HI * He, data = field_mice)
summary(model1)

# Model 2: Main effects only
model2 <- lm(predicted_weight_loss ~ Sex + HI + He + infection, data = field_mice)
summary(model2)

# Model 3: Two-way interactions
model3 <- lm(predicted_weight_loss ~ Sex * HI + Sex * He + HI * He +
               Sex * infection + HI * infection + He * infection,
             data = field_mice)
summary(model3)

# Model 4: Three-way interactions (full model)
model4 <- lm(predicted_weight_loss ~ Sex * HI * He * infection, data = field_mice)
summary(model4)

# Model 5: Focused model (biologically motivated)
model5 <- lm(predicted_weight_loss ~ Sex * (HI + He) * infection + HI:He,
             data = field_mice)

summary(model5)

# Model 3 updated: Two-way interactions using species_Eimeria instead of Sex and infection
modele6 <- lm(predicted_weight_loss ~ species_Eimeria * HI * He,
              data = field_mice)

summary(modele6)


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
best_model <- model3  # or choose based on AIC
cat("\nUsing focused model for further analysis\n")



detailed_models <- data.frame(
  Model = c("Model 1: Genetic baseline", "", "", "", "", "",  # Added 6th blank
            "Model 2: Main effects", "", "", "", "", "", "",  # Added 7th blank
            "Model 3: Two-way interactions", "", "", "", "", "", "", "", "", "", "", "", ""),  # 13 total
  Predictors = c("", "Intercept", "HI", "He", "HI × He", "",
                 "", "Intercept", "Sex (Male)", "HI", "He", "Infection", "",
                 "", "Intercept", "Sex (Male)", "HI", "He", "Infection",
                 "Sex × HI", "Sex × He", "Sex × Infection", "HI × He",
                 "HI × Infection", "He × Infection", ""),  # 26 entries
  Coefficient = c("", "9.570", "-0.154", "-0.392", "4.015", "",
                  "", "9.053", "-0.090", "0.185", "1.556", "1.125", "",
                  "", "9.316", "-1.308", "-0.670", "1.286", "2.320",
                  "2.144", "0.547", "-0.369", "1.692", "-1.299", "-0.908", ""),
  SE = c("", "0.411", "0.666", "2.035", "3.436", "",
         "", "0.372", "0.278", "0.406", "1.014", "0.285", "",
         "", "0.514", "0.736", "0.828", "2.379", "0.795",
         "0.833", "2.076", "0.571", "3.508", "0.841", "2.095", ""),
  t_value = c("", "23.27", "-0.23", "-0.19", "1.17", "",
              "", "24.35", "-0.32", "0.46", "1.53", "3.94", "",
              "", "18.14", "-1.78", "-0.81", "0.54", "2.92",
              "2.57", "0.26", "-0.65", "0.48", "-1.55", "-0.43", ""),
  p_value = c("", "<0.001", "0.818", "0.847", "0.244", "",
              "", "<0.001", "0.746", "0.649", "0.126", "<0.001", "",
              "", "<0.001", "0.077", "0.419", "0.589", "0.004",
              "0.011", "0.793", "0.519", "0.630", "0.123", "0.665", ""),
  Significance = c("", "***", "", "", "", "",
                   "", "***", "", "", "", "***", "",
                   "", "***", ".", "", "", "**",
                   "*", "", "", "", "", "", ""),
  Model_Stats = c("Residual SE: 2.47", "", "", "", "", "df: 300",
                  "Residual SE: 2.42", "", "", "", "", "", "df: 299",
                  "Residual SE: 2.40", "", "", "", "", "",
                  "", "", "", "", "", "", "df: 293")
)

# Check dimensions
cat("Dimensions check:\n")
cat("Rows:", nrow(detailed_models), "\n")
cat("Columns:", ncol(detailed_models), "\n")

# Create detailed supplementary table
table_detailed <- detailed_models %>%
  gt() %>%
  tab_header(
    title = "",
    subtitle = "Complete parameter estimates for models testing genetic and sex-specific effects on immune-predicted weight loss (n = 304)"
  ) %>%
  cols_label(
    Model = "Model",
    Predictors = "Parameters",
    Coefficient = "β",
    SE = "SE",
    t_value = "t",
    p_value = "p",
    Significance = "",
    Model_Stats = "Model Info"
  ) %>%
  # Merge SE into Coefficient column for cleaner look
  cols_merge(
    columns = c(Coefficient, SE),
    pattern = "{1} ({2})"
  ) %>%
  # Model headers with summary stats
  tab_style(
    style = list(
      cell_fill(color = "#e8eaf6"),
      cell_text(weight = "bold", size = 13)
    ),
    locations = cells_body(
      rows = c(1, 7, 14)
    )
  ) %>%
  # Update model names with fit statistics
  text_transform(
    locations = cells_body(columns = Model, rows = 1),
    fn = function(x) paste0(x, "<br><span style='font-size:11px; font-style:italic'>R² = 0.019, F(3,300) = 1.94, p = 0.122, AIC = 1419</span>")
  ) %>%
  text_transform(
    locations = cells_body(columns = Model, rows = 7),
    fn = function(x) paste0(x, "<br><span style='font-size:11px; font-style:italic'>R² = 0.063, F(4,299) = 5.05, p < 0.001, AIC = 1407</span>")
  ) %>%
  text_transform(
    locations = cells_body(columns = Model, rows = 14),
    fn = function(x) paste0(x, "<br><span style='font-size:11px; font-style:italic'>R² = 0.093, F(10,293) = 2.99, p = 0.001, AIC = 1409</span>")
  ) %>%
  # Highlight significant effects
  tab_style(
    style = list(
      cell_text(weight = "bold", color = "#1976d2")
    ),
    locations = cells_body(
      columns = c(Coefficient, t_value, p_value),
      rows = c(12, 19, 20)  # Infection and Sex × HI rows
    )
  ) %>%
  # Add significance stars
  cols_merge(
    columns = c(p_value, Significance),
    pattern = "{1}{2}"
  ) %>%
  # Format numbers for t_value column
  fmt_number(
    columns = t_value,
    decimals = 2,
    rows = c(2:5, 8:12, 15:25)  # Only format numeric rows
  ) %>%
  # Add comprehensive footnotes
  tab_footnote(
    footnote = "HI = Hybrid Index (0 = pure M. m. domesticus, 1 = pure M. m. musculus)",
    locations = cells_body(columns = Predictors, rows = 3)
  ) %>%
  tab_footnote(
    footnote = "He = Expected heterozygosity = 2 × HI × (1 - HI), measuring genetic admixture",
    locations = cells_body(columns = Predictors, rows = 4)
  ) %>%
  tab_footnote(
    footnote = "",
    locations = cells_body(columns = Coefficient, rows = 20)
  ) %>%
  # Add interpretation section
  tab_source_note(
    source_note = md("Significance codes: ***p < 0.001, **p < 0.01, *p < 0.05, . p < 0.1<br>
    Models compared using likelihood ratio tests.")
  ) %>%
  # Column formatting
  cols_align(
    align = "center",
    columns = c(Coefficient, t_value, p_value)
  ) %>%
  cols_align(
    align = "left",
    columns = c(Model, Predictors)
  ) %>%
  cols_width(
    Model ~ px(200),
    Predictors ~ px(180),
    Coefficient ~ px(120),
    t_value ~ px(60),
    p_value ~ px(100),
    Model_Stats ~ px(120)
  ) %>%
  tab_options(
    table.font.size = 11,
    heading.title.font.size = 14,
    heading.subtitle.font.size = 12,
    source_notes.font.size = 10,
    footnotes.font.size = 10
  )

# Save the table
save_table_all_formats(table_detailed, "table_S3_detailed_regression")


# Calculate AIC values and model statistics
# Note: Using the exact models from your output
models <- list(
  model1 = lm(predicted_weight_loss ~ HI * He, data = field_mice),  # n=335
  model2 = lm(predicted_weight_loss ~ Sex + HI + He + infection, data = field_mice),  # n=304 (31 missing)
  model3 = lm(predicted_weight_loss ~ Sex * HI + Sex * He + HI * He +
                Sex * infection + HI * infection + He * infection, data = field_mice)  # n=304 (31 missing)
)

# Extract model information
extract_model_info <- function(model, model_name) {
  n_obs <- nobs(model)
  k <- length(coef(model)) + 1  # +1 for residual variance
  aic_val <- AIC(model)
  loglik <- logLik(model)[1]

  # Get key coefficients with 85% CIs (Proceedings B standard)
  coef_summary <- summary(model)$coefficients

  # Extract Sex×HI interaction if present
  sex_hi_ci <- ""
  if("SexM:HI" %in% rownames(coef_summary)) {
    sex_hi_coef <- coef_summary["SexM:HI", "Estimate"]
    se <- coef_summary["SexM:HI", "Std. Error"]
    ci_lower <- sex_hi_coef - 1.44 * se
    ci_upper <- sex_hi_coef + 1.44 * se
    sex_hi_ci <- sprintf("%.2f (%.2f-%.2f)", sex_hi_coef, ci_lower, ci_upper)
  }

  # Extract infection coefficient if present
  infection_ci <- ""
  if("infectionInfected" %in% rownames(coef_summary)) {
    infection_coef <- coef_summary["infectionInfected", "Estimate"]
    se <- coef_summary["infectionInfected", "Std. Error"]
    ci_lower <- infection_coef - 1.44 * se
    ci_upper <- infection_coef + 1.44 * se
    infection_ci <- sprintf("%.2f (%.2f-%.2f)", infection_coef, ci_lower, ci_upper)
  }

  # Extract HI×He interaction if present
  hi_he_ci <- ""
  if("HI:He" %in% rownames(coef_summary)) {
    hi_he_coef <- coef_summary["HI:He", "Estimate"]
    se <- coef_summary["HI:He", "Std. Error"]
    ci_lower <- hi_he_coef - 1.44 * se
    ci_upper <- hi_he_coef + 1.44 * se
    hi_he_ci <- sprintf("%.2f (%.2f-%.2f)", hi_he_coef, ci_lower, ci_upper)
  }

  return(data.frame(
    model = model_name,
    n = n_obs,
    k = k,
    loglik = loglik,
    aic = aic_val,
    sex_hi_interaction = sex_hi_ci,
    infection_effect = infection_ci,
    hi_he_interaction = hi_he_ci,
    stringsAsFactors = FALSE
  ))
}

# Create model comparison dataframe
model_data <- do.call(rbind, lapply(names(models), function(x) {
  extract_model_info(models[[x]], x)
}))

# Calculate delta AIC and Akaike weights
model_data$delta_aic <- model_data$aic - min(model_data$aic)
model_data$akaike_weight <- exp(-0.5 * model_data$delta_aic) / sum(exp(-0.5 * model_data$delta_aic))

# Create the final table for display
proceedings_table_data <- data.frame(
  Model = c("M1", "M2", "M3"),
  Description = c(
    "HI × He",
    "Sex + HI + He + Infection",
    "Sex × HI + interactions"
  ),
  n = model_data$n,
  K = model_data$k,
  LogLik = sprintf("%.1f", model_data$loglik),
  AIC = sprintf("%.1f", model_data$aic),
  Delta_AIC = sprintf("%.1f", model_data$delta_aic),
  wi = sprintf("%.2f", model_data$akaike_weight),
  Key_Parameters = c(
    ifelse(model_data$hi_he_interaction[1] != "",
           model_data$hi_he_interaction[1],
           "--"),
    ifelse(model_data$infection_effect[2] != "",
           model_data$infection_effect[2],
           "--"),
    ifelse(model_data$sex_hi_interaction[3] != "",
           model_data$sex_hi_interaction[3],
           "--")
  ),
  stringsAsFactors = FALSE
)

# Create the gt table
proceedings_table <- proceedings_table_data %>%
  gt() %>%
  cols_label(
    Model = "Model",
    Description = "Description",
    n = "n",
    K = "K",
    LogLik = "LogLik",
    AIC = "AIC",
    Delta_AIC = "ΔAIC",
    wi = "wi",
    Key_Parameters = "Key Parameters (85% CI)"
  ) %>%
  # Highlight best model
  tab_style(
    style = list(
      cell_fill(color = "#e3f2fd"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(rows = which.min(model_data$delta_aic))
  ) %>%
  # Center numeric columns
  cols_align(
    align = "center",
    columns = c(n, K, LogLik, AIC, Delta_AIC, wi)
  ) %>%
  # Formatting
  cols_width(
    Model ~ px(60),
    Description ~ px(180),
    n ~ px(50),
    K ~ px(50),
    LogLik ~ px(80),
    AIC ~ px(80),
    Delta_AIC ~ px(80),
    wi ~ px(60),
    Key_Parameters ~ px(180)
  ) %>%
  tab_options(
    table.font.size = 11,
    table.border.top.style = "solid",
    table.border.bottom.style = "solid"
  )

# Print the table
proceedings_table

# Save the table using the provided function
save_table_all_formats(proceedings_table, "table_2_proceedings_b_model_selection")

# Also create a simplified version for quick reference
model_summary_supp <- data.frame(
  Comparison = c("Model 1 vs Null",
                 "Model 2 vs Model 1",
                 "Model 3 vs Model 2"),
  Hypothesis = c("Do genetic factors alone explain variation?",
                 "Do sex and infection explain variation?",
                 "Do interactions improve the model?"),
  Result = c("No (p = 0.122)",
             "Yes - Infection matters (p < 0.001)",
             "Yes - Sex × HI interaction (p = 0.011)"),
  Interpretation = c("Genetics alone insufficient",
                     "Infection is primary driver",
                     "Males specifically affected by ancestry")
)

table_summary_supp <- model_summary_supp %>%
  gt() %>%
  tab_header(
    title = "Supplementary Table S4. Model selection summary"
  ) %>%
  tab_style(
    style = cell_fill(color = "#e8f5e9"),
    locations = cells_body(rows = c(2, 3))
  ) %>%
  tab_options(table.font.size = 12)

save_table_all_formats(table_summary_supp, "table_S4_model_summary")

# ==============================================================================
# PUBLICATION FIGURE FOR MODEL RESULTS
# ==============================================================================

model2

# Panel A: Infection effect (main finding)
# Create prediction data
infection_pred <- expand.grid(
  infection = c("Uninfected", "Infected"),
  Sex = c("F", "M"),
  HI = mean(field_mice$HI, na.rm = TRUE),
  He = mean(field_mice$He, na.rm = TRUE)
)

infection_pred$predicted <- predict(model2, newdata = infection_pred)
infection_pred$se <- predict(model2, newdata = infection_pred, se.fit = TRUE)$se.fit

# Calculate means for main effect
infection_means <- infection_pred %>%
  group_by(infection) %>%
  summarise(
    mean_pred = mean(predicted),
    se = mean(se)
  )

panel_a <- ggplot(infection_means, aes(x = infection, y = mean_pred)) +
  geom_bar(stat = "identity", fill = c("#00FFFF", "#FF7094"), alpha = 0.7, width = 0.6) +
  geom_errorbar(aes(ymin = mean_pred - 1.96*se, ymax = mean_pred + 1.96*se),
                width = 0.1, size = 0.8) +
  geom_point(data = field_mice %>% filter(!is.na(infection)),
             aes(x = infection, y = predicted_weight_loss),
             position = position_jitter(width = 0.2), alpha = 0.3, size = 1) +
  labs(x = "Infection Status",
       y = "Predicted Weight Loss (%)",
       title = "") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  # Add statistics annotation
  annotate("text", x = 1.5, y = 11.5,
           label = "β = 1.12, p < 0.001",
           size = 4, fontface = "italic")

panel_a

save_plot_all_formats(plot_object = panel_a, plot_name = "predicted_weight_loss_infection")

# Panel B: Model predictions by sex and infection
panel_b <- ggplot(infection_pred, aes(x = infection, y = predicted, fill = Sex)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, alpha = 0.7) +
  geom_errorbar(aes(ymin = predicted - 1.96*se, ymax = predicted + 1.96*se),
                position = position_dodge(0.8), width = 0.2, size = 0.8) +
  scale_fill_manual(values = c("F" = "#4daf4a", "M" = "#ff7f00"),
                    labels = c("Female", "Male")) +
  labs(x = "Infection Status",
       y = "Predicted Weight Loss (%)",
       title = "",
       fill = "Sex") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = c(0.2, 0.9),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold")
  ) +
  annotate("text", x = 2, y = 11.8,
           label = "No interaction\np = 0.52",
           size = 3.5, fontface = "italic")

panel_b

save_plot_all_formats(plot_object = panel_b, plot_name = "sex_effect_predicted_WL")

# Create predictions across HI range
hi_range <- seq(0, 1, by = 0.01)
hi_pred <- expand.grid(
  HI = hi_range,
  Sex = c("F", "M"),
  He = 2 * hi_range * (1 - hi_range),  # Calculate He for each HI
  infection = "Infected"  # Fix infection status
)

hi_pred$predicted <- predict(model3, newdata = hi_pred)
hi_pred$se <- predict(model3, newdata = hi_pred, se.fit = TRUE)$se.fit

panel_c <- ggplot(hi_pred, aes(x = HI, y = predicted, color = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = predicted - 1.96*se, ymax = predicted + 1.96*se),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  geom_point(data = field_mice %>% filter(!is.na(infection)),
             aes(x = HI, y = predicted_weight_loss, color = Sex),
             alpha = 0.3, size = 1) +
  scale_color_manual(values = c("F" = "#E69F00", "M" = "#0072B2"),
                     labels = c("Female", "Male")) +
  scale_fill_manual(values = c("F" = "#E69F00", "M" = "#0072B2"),
                    labels = c("Female", "Male")) +
  labs(x = "Hybrid Index",
       y = "Predicted Weight Loss (%)",
       title = "C) Sex × HI Interaction") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11),
    legend.position = c(0.2, 0.9),
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(face = "bold")
  ) +
  annotate("text", x = 0.5, y = 13,
           label = "Sex × HI: p = 0.011",
           size = 4, fontface = "italic")

panel_c

# Panel D: Model comparison
model_comp <- data.frame(
  Model = factor(c("M1: HI × He", "M2: Main effects", "M3: Interactions"),
                 levels = c("M1: HI × He", "M2: Main effects", "M3: Interactions")),
  AIC = c(1557, 1407, 1409),
  R2 = c(0.019, 0.063, 0.093)
)

panel_d <- ggplot(model_comp, aes(x = Model, y = AIC)) +
  geom_bar(stat = "identity", fill = c("gray70", "#4CAF50", "#2196F3"), alpha = 0.8) +
  geom_text(aes(label = paste0("R² = ", R2)),
            vjust = -0.5, size = 3.5, fontface = "bold") +
  labs(x = "",
       y = "AIC",
       title = "D) Model Comparison") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 20, hjust = 1),
    axis.text.y = element_text(size = 11)
  ) +
  coord_cartesian(ylim = c(1350, 1600))

# Combine all panels
figure1 <- grid.arrange(panel_a, panel_b, panel_c, panel_d,
                        ncol = 2, nrow = 2,
                        padding = unit(0.5, "line"))

# Model 3: Two-way interactions uninfected

# Create uninfected subset for constitutive cost analysis
uninfected_data <- field_mice %>%
  filter(infection_status == "FALSE")

model3_uninfected <- lm(predicted_weight_loss ~ Sex * HI + Sex * He + HI * He +
               Sex  + HI  + He,
             data = uninfected_data)
summary(model3_uninfected)


library(ggfortify)
p <- autoplot(model3, which = 1:4, ncol = 2)

dir.create("Results/Figures/supplementary", recursive = TRUE, showWarnings = FALSE)

grDevices::cairo_pdf("results/figures/Supp_Fig4_model3_diagnostics.pdf",
                     width = 7.5, height = 7.5)
par(mfrow = c(2,2), mar = c(4,4,2,1))
plot(model3)     # Residuals vs Fitted, QQ, Scale-Location, Residuals vs Leverage
dev.off()

