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
  scale_color_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094"),
                     name = "Infection Status") +
  scale_fill_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094"),
                    name = "Infection Status") +
  labs(title = "Predicted Weight Loss by Hybrid Index, Sex, and Infection",
       x = "Hybrid Index",
       y = "Predicted Weight Loss (%)") +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

plot_main_effects

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

plot_he_effects

save_plot_all_formats(plot_he_effects, "weight_loss_heterozygosity_effects")

# C. Coefficient plot
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


# C. Coefficient plot second best
coef_data_3 <- tidy(model3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = str_replace_all(term, ":", " × "),
    significant = p.value < 0.05
  )

plot_coefficients_3 <- ggplot(coef_data_3, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_point(aes(color = significant), size = 3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"),
                     name = "Significant\n(p < 0.05)") +
  labs(title = "",
       x = "Coefficient Estimate",
       y = "") +
  theme_minimal()

plot_coefficients_3

save_plot_all_formats(plot_coefficients_3, "weight_loss_coefficients_HI")

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

# Select top pathway based on adjusted p-value for Sex × HI interaction
top_pathway <- pathway_results %>%
  arrange(Sex_HI_padj) %>%
  slice(1) %>%
  pull(Pathway)

# Build model
formula_top <- as.formula(paste(top_pathway, "~ Sex * (HI + He) * infection + HI:He"))
model_top <- lm(formula_top, data = immune_data)
summary(model_top)

# Predict
pred_pathway <- ggpredict(model_top, terms = c("HI [all]", "infection", "Sex"))

# 1. Define your hybrid gradient bar as a function per facet
make_HI_gradient_bar <- function(label) {
  ggplot(data.frame(hi = seq(0, 1, 0.001)), aes(x = hi, y = 1, fill = hi)) +
    geom_tile() +
    scale_x_continuous(breaks = seq(0, 1, 0.25),
                       labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_void() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(color = "black", size = 10),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    labs(title = label)
}

# 2. Create both gradients manually for "F" and "M"
HIbar_F <- make_HI_gradient_bar("F")
HIbar_M <- make_HI_gradient_bar("M")

# 3. Patchwork layout: one plot row per panel, and one for each bar
# First, regenerate your plot without a shared x-axis title
main_plot <- plot_top_pathway +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm")
  )

# 4. Assemble all with patchwork
combined_plot <- main_plot /
  (HIbar_F | HIbar_M) +
  plot_layout(heights = c(1, 0.1))

# 5. Print and save
print(combined_plot)

save_plot_all_formats(combined_plot, paste0("pathway_", top_pathway, "_facet_gradient"))
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

plot_coherence

coherence_model <- lm(coherence ~ poly(mean_HI, 2), data = coherence_data)
summary(coherence_model)

save_plot_all_formats(plot_coherence, "network_coherence_by_infection")



# ==============================================================================
# 8. UMAP
# ==============================================================================
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

plot_umap_weight_loss

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

# Create dispersion plot showing HI × infection interaction
# This is highly significant (p = 0.00169)

# Prepare prediction data for dispersion model
pred_dispersion <- ggpredict(disp_model, terms = c("HI [all]", "infection", "Sex"))

# Create the dispersion plot
plot_dispersion <- ggplot(pred_dispersion, aes(x = x, y = predicted, color = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~ facet, nrow = 1) +
  scale_color_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094"),
                     name = "Infection Status") +
  scale_fill_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094"),
                    name = "Infection Status") +
  labs(title = "d) Immune Dispersion Shows HI × Infection Effect",
       subtitle = "HI × infection: p = 0.00169**",
       x = "Hybrid Index",
       y = "Distance to Centroid\n(Immune Dysregulation)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(face = "italic", color = "darkred"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

  plot_dispersion
# Save the plot
save_plot_all_formats(plot_dispersion, "immune_dispersion_HI_infection")

# Alternative: Show actual data points with model predictions
plot_dispersion_points <- ggplot(immune_data, aes(x = HI, y = dist_to_centroid)) +
  geom_point(aes(color = infection, shape = Sex), alpha = 0.6, size = 2) +
  geom_smooth(aes(color = infection), method = "lm", formula = y ~ poly(x, 2),
              se = TRUE, size = 1.2) +
  scale_color_manual(values = c("Uninfected" = "#00FFFF", "Infected" = "#FF7094"),
                     name = "Infection Status") +
  scale_shape_manual(values = c("F" = 16, "M" = 17), name = "Sex") +
  labs(title = "d) Immune Dispersion Shows HI × Infection Effect",
       subtitle = "HI × infection: p = 0.00169**",
       x = "Hybrid Index",
       y = "Distance to Centroid\n(Immune Dysregulation)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(face = "italic", color = "darkred")
  )

print(plot_dispersion_points)

manova_result <- manova(cbind(UMAP1, UMAP2) ~ Sex * HI * He * infection,
                        data = immune_data)
cat("MANOVA Results:\n")
print(summary(manova_result, test = "Pillai"))

##################
cat("\n\nENHANCED UMAP ANALYSIS: HYBRID TRAJECTORIES\n")
cat("============================================\n")

# Define correct sex colors
sex_colors <- c(
  "F" = "#4daf4a",   # Green
  "M" = "#ff7f00"    # Orange
)

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
      n_points <- nrow(traj_subset)
      plot_umap_trajectories <- plot_umap_trajectories +
        annotate("segment",
                 x = traj_subset$mean_UMAP1[n_points-1],
                 y = traj_subset$mean_UMAP2[n_points-1],
                 xend = traj_subset$mean_UMAP1[n_points],
                 yend = traj_subset$mean_UMAP2[n_points],
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
                 color = sex_colors[[s]],  # <- updated here
                 size = 1.5, alpha = 0.8)
    }
  }
}

print(plot_umap_trajectories)
save_plot_all_formats_panel(plot_umap_trajectories, "UMAP_hybrid_trajectories")

# Create a color lookup table
sex_color_map <- c("F" = sex_colors[["F"]], "M" = sex_colors[["M"]])

# Add a column for path color to trajectory_data
trajectory_data <- trajectory_data %>%
  mutate(sex_color = sex_color_map[Sex])

# Plot
plot_umap_weight_trajectory <- ggplot() +
  # Points colored by predicted weight loss (continuous)
  geom_point(data = immune_data,
             aes(x = UMAP1, y = UMAP2, color = predicted_weight_loss),
             size = 2, alpha = 0.6) +
  scale_color_viridis_c(name = "Predicted\nWeight Loss (%)", option = "plasma") +

  # Trajectory overlay with color applied manually (not via aes)
  geom_path(data = trajectory_data,
            aes(x = mean_UMAP1, y = mean_UMAP2,
                group = interaction(Sex, infection)),
            size = 1.2, alpha = 0.9,
            color = NA) +
  geom_path(data = trajectory_data,
            aes(x = mean_UMAP1, y = mean_UMAP2,
                group = interaction(Sex, infection)),
            size = 1.2, alpha = 0.9,
            color = trajectory_data$sex_color,
            inherit.aes = FALSE) +

  facet_wrap(~ infection) +
  labs(title = "Weight Loss Landscape in Immune Space",
       subtitle = "Points: predicted weight loss; Lines: sex-specific hybrid trajectories",
       x = "UMAP Dimension 1",
       y = "UMAP Dimension 2") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA, color = "black"))

print(plot_umap_weight_trajectory)

save_plot_all_formats(plot_umap_weight_trajectory, "UMAP_weight_loss_trajectories")
cat("=============================================\n")

# Best UMAP visualization for Sex × HI × He interaction
plot_umap_trajectories <- ggplot() +
  # Individual points in background
  geom_point(data = immune_data,
             aes(x = UMAP1, y = UMAP2, color = HI),
             alpha = 0.3, size = 2) +
  # Trajectory lines for each sex
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
                       midpoint = 0.5, guide = "none") +
  scale_shape_manual(values = c("F" = 21, "M" = 24), name = "Sex") +
  facet_wrap(~ infection, nrow = 1) +  # Shows infected vs uninfected
  labs(
    title = "c) Sex-Specific Immune Trajectories",
    subtitle = "Sex × HI × He interaction: p = 0.020*",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(face = "italic", color = "darkred"),
    panel.border = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank()
  )

plot_umap_trajectories

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
  scale_fill_manual(values = c("Uninfected" = "#00FFFF",
                               "Infected" = "#FF7094")) +
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

