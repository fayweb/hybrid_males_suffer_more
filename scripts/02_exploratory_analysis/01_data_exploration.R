# ==============================================================================
# FIGURE 1: Experimental Design and Data Overview
# ==============================================================================
#
# Purpose: Create Figure 1 for the manuscript showing:
# - Study design and sample composition
# - Hybrid index distribution across the zone
# - Sex distribution and infection status
# - Predicted weight loss from Chapter 1 model
#
# Author: Fay Webster
# Date: June 2025
# ==============================================================================

# Load required functions and data (assumes master script has been run)
if (!exists("field_mice")) {
  stop("Please run 00_master_script.R first to load data and packages")
}

cat("=== CREATING FIGURE 1: DATA OVERVIEW ===\n")
cat("Dataset: field_mice (n =", nrow(field_mice), ")\n\n")

# ==============================================================================
# 1. EXAMINE KEY VARIABLES
# ==============================================================================

cat("1. EXAMINING KEY VARIABLES\n")
cat("===========================\n")

# Check the actual key variables we need
cat("Hybrid Index (HI) summary:\n")
print(summary(field_mice$HI))

cat("\nSex distribution:\n")
print(table(field_mice$Sex, useNA = "ifany"))

cat("\nPredicted weight loss summary:\n")
print(summary(field_mice$predicted_weight_loss))

cat("\nEimeria infection status (MC.Eimeria):\n")
print(table(field_mice$MC.Eimeria, useNA = "ifany"))

cat("\nEimeria species:\n")
print(table(field_mice$species_Eimeria, useNA = "ifany"))

# ==============================================================================
# 2. DATA COMPLETENESS FOR HYBRID ANALYSIS
# ==============================================================================

cat("\n\n2. DATA COMPLETENESS\n")
cat("====================\n")

# Key variables for hybrid analysis
key_vars <- c("HI", "Sex", "predicted_weight_loss", "MC.Eimeria")

# Check completeness
completeness <- data.frame(
  variable = key_vars,
  missing_n = sapply(key_vars, function(x) sum(is.na(field_mice[[x]]))),
  missing_pct = sapply(key_vars, function(x) round(sum(is.na(field_mice[[x]])) / nrow(field_mice) * 100, 1))
)

print(completeness)

# Complete cases for hybrid analysis
complete_cases <- complete.cases(field_mice[key_vars])
complete_n <- sum(complete_cases)
complete_pct <- round(complete_n / nrow(field_mice) * 100, 1)

cat(sprintf("\nComplete cases for hybrid analysis: %d/%d (%s%%)\n",
            complete_n, nrow(field_mice), complete_pct))

# ==============================================================================
# 3. CREATE FIGURE 1 COMPONENTS
# ==============================================================================

cat("\n\n3. CREATING FIGURE 1 PANELS\n")
cat("=============================\n")

number_hybrids <-
  ggplot(field_mice, aes(x = HI)) +
  geom_histogram(bins = 30, fill = "gray80", color = "white") +
  geom_rug(aes(color = HI), sides = "b", alpha = 0.8, size = 1) +
  geom_vline(xintercept = 0.5, linetype = "dashed", size = 1) +
  scale_color_gradientn(
    colors = c("blue", "purple", "red"),
    values = scales::rescale(c(0, 0.5, 1)),
    name = "Hybrid Index"
  ) +
  annotate("text", x = 0.05, y = Inf, label = "M.m.domesticus",
           color = "blue", size = 4, vjust = 2, hjust = 0, fontface = "italic") +
  annotate("text", x = 0.95, y = Inf, label = "M.m.musculus",
           color = "red", size = 4, vjust = 2, hjust = 1, fontface = "italic") +
  labs(
   # title = "A) Hybrid Index with Individual Gradient",
    x = "Hybrid Index (0 = M.m.domesticus, 1 = M.m.musculus)",
    y = "Number of mice"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "right"
  )

number_hybrids

save_plot_all_formats(number_hybrids, plot_name = "Number_mice_hybridization_barplot")



hybridization_gradient <-
  ggplot(field_mice, aes(x = HI, y = 0, color = HI)) +
  geom_jitter(height = 0.05, width = 0.01, alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black", size = 1) +
  scale_color_gradientn(
    colors = c("blue", "purple", "red"),
    values = scales::rescale(c(0, 0.5, 1)),
    name = "Hybrid Index"
  ) +
  annotate("text", x = 0.05, y = 0.08, label = "M.m.domesticus",
           color = "blue", size = 4, fontface = "italic", hjust = 0) +
  annotate("text", x = 0.95, y = 0.08, label = "M.m.musculus",
           color = "red", size = 4, fontface = "italic", hjust = 1) +
  labs(
   # title = "A) Hybridization Gradient (n = 336)",
    x = "Hybrid Index (0 = M.m.domesticus, 1 = M.m.musculus)",
    y = NULL,
    caption = "Dashed line: hybrid center (HI = 0.5)"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold")
  )

hybridization_gradient

save_plot_all_formats(hybridization_gradient, plot_name = "Hybridization_gradient_barplot")


# Panel B: Infection Detection Success
# Show the infection detection pipeline results
detection_summary <- data.frame(
  Method = c("qPCR Success", "qPCR data missing", "Amplicon Success", "Final Complete"),
  n = c(185, 151, 134, 169),
  Type = c("Primary", "Primary", "Backup", "Combined")
)

panel_b <- ggplot(detection_summary, aes(x = Method, y = n, fill = Type)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = n), vjust = -0.5, size = 4, fontface = "bold") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  ) +
  labs(
    #title = "B) Infection Detection Pipeline",
    x = "Detection Method",
    y = "Number of mice",
    fill = "Method Type"
  ) +
  scale_fill_manual(values = c("Primary" = "#4CAF50", "Backup" = "#FF9800", "Combined" = "#2196F3")) +
  ylim(0, max(detection_summary$n) * 1.1)

print(panel_b)

save_plot_all_formats(panel_b, plot_name = "Eimeria_quantification_method")


# Panel C: Final Dataset Composition (Complete Cases)
# Use infection_status for most complete dataset
final_data_summary <- field_mice %>%
  filter(!is.na(infection_status)) %>%
  mutate(
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "Male")),
    infection_status = factor(infection_status, levels = c("FALSE", "TRUE"),
                              labels = c("Uninfected", "Infected with Eimeria spp."))
  ) %>%
  count(Sex, infection_status)

# Proportional distribution by infection status
# Refined plot of final dataset by sex and infection status
panel_c <- ggplot(final_data_summary, aes(x = Sex, y = n, fill = Sex)) +
  geom_col(width = 0.6, alpha = 0.7, color = "black") +
  geom_text(aes(label = n), vjust = -0.6, size = 4, fontface = "bold") +
  facet_wrap(~infection_status) +
  scale_fill_manual(values = sex_colors, name = NULL) +
  labs(
   # title = paste0("C) Final Analysis Dataset (n = ", sum(final_data_summary$n), ")"),
    x = NULL,
    y = "Number of mice"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

print(panel_c)

infection_status_sex <-
  ggplot(final_data_summary, aes(x = infection_status, y = n, fill = Sex)) +
  geom_col(position = "fill", color = "black", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = sex_colors) +
  labs(
   # title = "Sex Distribution Within Infection Groups",
    x = "Infection Status",
    y = "Proportion of Mice",
    fill = "Sex"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 13),
    strip.text = element_text(size = 12)
  )

save_plot_all_formats(panel_c, plot_name = "Infection_status_sex_distribution")


# Panel D: Predicted Weight Loss by Sex and Final Infection Status
panel_d <- field_mice %>%
  filter(!is.na(infection_status)) %>%
  mutate(
    Sex = factor(Sex, levels = c("F", "M"), labels = c("Female", "Male")),
    infection_status = factor(infection_status, levels = c("FALSE", "TRUE"),
                              labels = c("Uninfected", "Infected with Eimeria spp."))
  ) %>%
  ggplot(aes(x = Sex, y = predicted_weight_loss, fill = Sex)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.6) +
  facet_wrap(~infection_status) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    #title = "D) Predicted Weight Loss by Sex and Infection",
    x = "Sex",
    y = "Predicted weight loss (%)",
    fill = "Sex"
  ) +
  scale_fill_manual(values = sex_colors)

print(panel_d)

save_plot_all_formats(panel_d, plot_name = "Sex_impact_Eimeria")

cat("\n\n4. TESTING SEX AND INFECTION EFFECTS ON PREDICTED WEIGHT LOSS\n")
cat("===============================================================\n")

# Linear model: predicted WL ~ Sex + infection
wl_model <- lm(predicted_weight_loss ~ Sex + infection_status, data = field_mice)
summary(wl_model)

# Optional: include interaction
wl_model_int <- lm(predicted_weight_loss ~ Sex * infection_status, data = field_mice)
anova(wl_model, wl_model_int)  # Test interaction effect

# Compare HI across sex and infection status
hi_model <- lm(HI ~ Sex + infection_status, data = field_mice)
summary(hi_model)

# Nonparametric check (e.g., if residuals look off)
wilcox.test(HI ~ infection_status, data = field_mice)  # Quick test


emmeans(wl_model, pairwise ~ Sex)
emmeans(wl_model, pairwise ~ infection_status)

# Clean model outputs
wl_model_tidy <- tidy(wl_model) %>%
  mutate(model = "Predicted weight loss")

hi_model_tidy <- tidy(hi_model) %>%
  mutate(model = "Hybrid index")

supp_table <- bind_rows(wl_model_tidy, hi_model_tidy) %>%
  dplyr::select(model, term, estimate, std.error, statistic, p.value) %>%
  mutate(
    term = case_when(
      term == "SexM" ~ "Male",
      term == "infection_statusTRUE" ~ "Infected with *Eimeria* spp.",
      TRUE ~ term
    ),
    estimate = round(estimate, 3),
    std.error = round(std.error, 3),
    statistic = round(statistic, 2),
    p.value = signif(p.value, 3)
  )

# Create gt table
gt_supp_table <- supp_table %>%
  gt() %>%
  tab_header(
   title = md("")
  ) %>%
  cols_label(
    model = "Model",
    term = "Term",
    estimate = "Estimate",
    std.error = "SE",
    statistic = "t-value",
    p.value = "p-value"
  ) %>%
  fmt_markdown(columns = term) %>%  # for italic Eimeria
  fmt_number(columns = c(estimate, std.error), decimals = 3) %>%
  fmt_number(columns = statistic, decimals = 2) %>%
  fmt_number(columns = p.value, decimals = 3) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold"
  )

gt_supp_table

save_table_all_formats(gt_supp_table, table_name = "Sex_Infection_model_results")


model_data <- bind_rows(wl_model_tidy, hi_model_tidy) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term = case_when(
      term == "SexM" ~ "Male sex",
      term == "infection_statusTRUE" ~ "Infected (vs. uninfected)",
      TRUE ~ term
    ),
    p_signif = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# Create the coefficient plot
model_data_clean <- model_data %>%
  mutate(
    Response = model,
    Variable = factor(term, levels = c("Infected (vs. uninfected)", "Male sex")),
    Variable = fct_recode(Variable,
                          "Infected with *Eimeria* spp." = "Infected (vs. uninfected)",
                          "Male" = "Male sex"
    ),
    Estimate_label = paste0(round(estimate, 2), " ± ", round(std.error, 2), " ", p_signif),
    term_color = case_when(
      term == "Infected (vs. uninfected)" ~ infection_status_colors["Infected with Eimeria spp."],
      term == "Male sex" ~ sex_colors["Male"],
      TRUE ~ "black"
    )
  )

# Create coefficient plot
coeff_plot <- ggplot(model_data_clean, aes(x = estimate, y = Variable, color = Variable)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = estimate - std.error, xmax = estimate + std.error),
                 height = 0.2, size = 1) +
  geom_point(size = 3) +
  facet_wrap(~Response, scales = "free_x") +
  scale_color_manual(values = c(
    "Infected with *Eimeria* spp." = infection_status_colors["Infected with Eimeria spp."],
    "Male" = sex_colors["Male"]
  )) +
  labs(
    x = "Effect size (estimate ± SE)",
    y = NULL,
    title = "Effects of Sex and Infection on Health Outcomes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 13)
  )

print(coeff_plot)


Sex_Infection_model_plot <-
  ggplot(model_data_clean, aes(x = estimate, y = Variable)) +
  geom_col(aes(fill = term_color), width = 0.6, color = "black") +
  geom_errorbar(
    aes(xmin = estimate - std.error, xmax = estimate + std.error),
    width = 0.2
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30") +
  geom_text(
    aes(label = Estimate_label),
    hjust = -0.1,
    size = 4
  ) +
  facet_wrap(~Response, scales = "free_x") +
  scale_fill_identity(guide = "none") +
  labs(
    x = "Effect size (estimate ± SE)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.y = element_markdown(size = 12),  # Italic rendering for *Eimeria*
    legend.position = "none"
  )


Sex_Infection_model_plot

save_plot_all_formats(Sex_Infection_model_plot, plot_name = "Sex_Infection_model_plot")
# ==============================================================================
# 5. CREATE SUMMARY TABLE
# ==============================================================================

cat("\n5. CREATING SUMMARY TABLE\n")
cat("==========================\n")

# Create Table 1: Sample characteristics for final analysis dataset
final_analysis_data <- field_mice %>%
  filter(!is.na(infection_status))

table1 <- data.frame(
  Characteristic = c(
    "Total sample size (complete data)",
    "Total sample size (all mice)",
    "Females (final dataset)",
    "Males (final dataset)",
    "Hybrid index range",
    "Hybrid index median (IQR)",
    "Eimeria infected (final)",
    "Eimeria uninfected (final)",
    "E. ferrisi infected",
    "E. falciformis infected",
    "Predicted weight loss range (%)",
    "Predicted weight loss median (IQR)",
    "qPCR detection success rate",
    "Amplicon sequencing success rate"
  ),
  Value = c(
    nrow(final_analysis_data),
    nrow(field_mice),
    sum(final_analysis_data$Sex == "F", na.rm = TRUE),
    sum(final_analysis_data$Sex == "M", na.rm = TRUE),
    paste(round(range(final_analysis_data$HI, na.rm = TRUE), 3), collapse = " - "),
    paste0(round(median(final_analysis_data$HI, na.rm = TRUE), 3), " (",
           round(quantile(final_analysis_data$HI, 0.25, na.rm = TRUE), 3), " - ",
           round(quantile(final_analysis_data$HI, 0.75, na.rm = TRUE), 3), ")"),
    sum(final_analysis_data$infection_status == "TRUE", na.rm = TRUE),
    sum(final_analysis_data$infection_status == "FALSE", na.rm = TRUE),
    sum(final_analysis_data$species_Eimeria == "E_ferrisi", na.rm = TRUE),
    sum(final_analysis_data$species_Eimeria == "E_falciformis", na.rm = TRUE),
    paste(round(range(final_analysis_data$predicted_weight_loss, na.rm = TRUE), 2), collapse = " - "),
    paste0(round(median(final_analysis_data$predicted_weight_loss, na.rm = TRUE), 2), " (",
           round(quantile(final_analysis_data$predicted_weight_loss, 0.25, na.rm = TRUE), 2), " - ",
           round(quantile(final_analysis_data$predicted_weight_loss, 0.75, na.rm = TRUE), 2), ")"),
    "55.1% (185/336)",
    "79.3% (134/169)"
  )
)

print(table1)


cat("✓ Saved Table 1: Sample Characteristics\n")

knitr::kable(table1, caption = "Table 1: Summary of the Final Dataset", align = "ll", booktabs = TRUE)

# Replace relevant strings with HTML italic tags
table1$Characteristic <- table1$Characteristic %>%
  str_replace_all("E\\. ferrisi", "<i>E\\. ferrisi</i>") %>%
  str_replace_all("E\\. falciformis", "<i>E\\. falciformis</i>") %>%
  str_replace_all("Eimeria", "<i>Eimeria</i>")

# Create the GT table
gt_table1 <- table1 %>%
  gt() %>%
  tab_header(
    title = html("")
  ) %>%
  cols_label(
    Characteristic = "Characteristic",
    Value = "Value"
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 14,
    heading.title.font.weight = "bold"
  ) %>%
  fmt_markdown(columns = "Characteristic")  # Enable rendering italics

# View the table
gt_table1
# Save all formats
save_table_all_formats(gt_table1, table_name = "Sample_Characteristics_design")


#################### Plot hybrid index across sampling locations
# Sample points to sf
sample_points <- field_mice %>%
  filter(!is.na(HI), !is.na(Longitude), !is.na(Latitude)) %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# Load base map
europe <- ne_countries(continent = "Europe", returnclass = "sf", scale = "medium")

# Custom color scale from blue (musculus) → purple (hybrid) → red (domesticus)
hybrid_gradient <- scale_color_gradientn(
  name = "Hybrid Index",
  colors = c("blue", "purple", "firebrick1"),
  values = scales::rescale(c(1, 0.5, 0)),  # musculus (1) to domesticus (0)
  limits = c(0, 1)
)

ggplot() +
  geom_sf(data = europe, fill = "gray95", color = "gray80") +
  geom_sf(data = sample_points, aes(color = HI), size = 2.5, alpha = 0.9) +
  scale_color_gradientn(
    name = "Hybrid Index",
    colors = c("blue", "purple", "firebrick1"),
    values = scales::rescale(c(1, 0.5, 0)),  # musculus (1) → hybrid (0.5) → domesticus (0)
    limits = c(0, 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(11.5, 15), ylim = c(51, 54), expand = FALSE) +
  theme_minimal(base_size = 13) +
  labs(
    title = "A) Sampling Locations Across the House Mouse Hybrid Zone",
    x = "Longitude", y = "Latitude"
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(size = 15, face = "bold"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

# Convert to EPSG:3857
sample_points_merc <- sample_points %>% st_transform(crs = 3857)

# Plot using OpenStreetMap terrain tiles
hybrid_map <-
  ggplot() +
  annotation_map_tile(type = "osm") +  # Or use "cartolight", "cartodark", "stamenbw", "esri.topo"
  geom_sf(data = sample_points_merc, aes(color = HI), size = 2.5, alpha = 0.9) +
  scale_color_gradientn(
    name = "Hybrid Index",
    colors = c("blue", "purple", "firebrick1"),
    values = scales::rescale(c(1, 0.5, 0)),
    limits = c(0, 1)
  ) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering) +
  labs(
    #title = "A) Sampling Locations Across the House Mouse Hybrid Zone",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal(base_size = 13)

hybrid_map

save_plot_all_formats(hybrid_map, plot_name = "Hybrid_index_locations")
# ==============================================================================
# 6. SUMMARY
# ==============================================================================

cat("\n\n=== FIGURE 1 CREATION COMPLETE ===\n")
cat("Created panels showing:\n")
cat("A) Hybrid index distribution across the zone\n")
cat("B) Sample composition by sex and infection status\n")
cat("C) Distribution of predicted weight loss values\n")
cat("D) Weight loss patterns by sex and infection\n\n")

cat("Files created:\n")
cat("- Figure1_Data_Overview.png/pdf (publication figure)\n")
cat("- Table1_Sample_Characteristics.csv (manuscript table)\n\n")

cat("Ready for hybrid analysis using Alice Balard's framework!\n")

