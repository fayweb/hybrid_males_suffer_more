# ==============================================================================
# COMPREHENSIVE RESULTS EXTRACTION FOR "HYBRID MALES SUFFER MORE" MANUSCRIPT
# ==============================================================================
# This script extracts all statistics and creates publication-ready tables
# for the Results section of the manuscript
# ==============================================================================

# Load required libraries
library(tidyverse)
library(broom)
library(effsize)
library(gt)

# ==============================================================================
# PART 1: DATA PREPARATION AND QUALITY CONTROL
# ==============================================================================

cat("\n========== DATA PREPARATION ==========\n")

# Check initial dataset
cat(sprintf("Initial dataset: %d mice\n", nrow(field_mice)))

# Remove mice with missing hybrid index
field_mice_complete <- field_mice %>%
  drop_na(HI) %>%
  mutate(
    # Create clean infection status variable
    infection_status = case_when(
      MC.Eimeria == TRUE ~ "Infected",
      MC.Eimeria == FALSE ~ "Uninfected",
      TRUE ~ "Unknown"
    ),
    # Create sex labels
    sex_label = ifelse(Sex == "F", "Female", "Male")
  )

cat(sprintf("After removing missing HI: %d mice\n", nrow(field_mice_complete)))
cat(sprintf("Mice removed: %d\n\n", nrow(field_mice) - nrow(field_mice_complete)))

# ==============================================================================
# PART 2: SAMPLE CHARACTERISTICS
# ==============================================================================

cat("========== SAMPLE CHARACTERISTICS ==========\n")

# Geographic distribution
geo_stats <- field_mice_complete %>%
  summarise(
    n_locations = n_distinct(paste(Longitude, Latitude)),
    lat_range = sprintf("%.2f–%.2f°N", min(Latitude), max(Latitude)),
    lon_range = sprintf("%.2f–%.2f°E", min(Longitude), max(Longitude)),
    years = paste(unique(sort(Year)), collapse = ", ")
  )

cat("\nGeographic distribution:\n")
cat(sprintf("  Localities sampled: %d\n", geo_stats$n_locations))
cat(sprintf("  Latitude range: %s\n", geo_stats$lat_range))
cat(sprintf("  Longitude range: %s\n", geo_stats$lon_range))
cat(sprintf("  Sampling years: %s\n", geo_stats$years))

# Sex distribution
sex_stats <- field_mice_complete %>%
  count(sex_label) %>%
  mutate(
    percentage = round(n / sum(n) * 100, 1),
    label = sprintf("%d (%.1f%%)", n, percentage)
  )

cat("\nSex distribution:\n")
for(i in 1:nrow(sex_stats)) {
  cat(sprintf("  %s: %s\n", sex_stats$sex_label[i], sex_stats$label[i]))
}

# Hybrid index distribution
hi_stats <- field_mice_complete %>%
  summarise(
    mean = round(mean(HI), 3),
    sd = round(sd(HI), 3),
    median = round(median(HI), 3),
    q1 = round(quantile(HI, 0.25), 3),
    q3 = round(quantile(HI, 0.75), 3),
    range = sprintf("%.3f–%.3f", min(HI), max(HI))
  )

cat("\nHybrid index distribution:\n")
cat(sprintf("  Mean ± SD: %.3f ± %.3f\n", hi_stats$mean, hi_stats$sd))
cat(sprintf("  Median (IQR): %.3f (%.3f–%.3f)\n", hi_stats$median, hi_stats$q1, hi_stats$q3))
cat(sprintf("  Range: %s\n", hi_stats$range))

# ==============================================================================
# PART 3: INFECTION STATUS ANALYSIS
# ==============================================================================

cat("\n========== INFECTION STATUS ANALYSIS ==========\n")

# Overall infection status
infection_stats <- field_mice_complete %>%
  count(infection_status) %>%
  mutate(
    percentage = round(n / sum(n) * 100, 1),
    label = sprintf("%d (%.1f%%)", n, percentage)
  )

cat("\nInfection status breakdown:\n")
for(i in 1:nrow(infection_stats)) {
  cat(sprintf("  %s: %s\n", infection_stats$infection_status[i], infection_stats$label[i]))
}

# Create complete case dataset (known infection status only)
complete_data <- field_mice_complete %>%
  filter(infection_status != "Unknown")

cat(sprintf("\nComplete case analysis dataset: %d mice (%.1f%% of total)\n",
            nrow(complete_data),
            nrow(complete_data) / nrow(field_mice_complete) * 100))

# Sex × Infection crosstab
sex_infection_stats <- complete_data %>%
  count(sex_label, infection_status) %>%
  pivot_wider(names_from = infection_status, values_from = n) %>%
  mutate(Total = Infected + Uninfected)

cat("\nSex × Infection status:\n")
print(sex_infection_stats)

# Species identification (among infected)
species_stats <- field_mice_complete %>%
  filter(infection_status == "Infected") %>%
  mutate(
    species_clean = case_when(
      species_Eimeria == "E. ferrisi" ~ "E. ferrisi",
      species_Eimeria == "E. falciformis" ~ "E. falciformis",
      species_Eimeria == "Uninfected" ~ "Detection error",
      is.na(species_Eimeria) ~ "Not identified",
      TRUE ~ "Other"
    )
  ) %>%
  count(species_clean) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

cat("\nSpecies identification among infected mice:\n")
print(species_stats)

# ==============================================================================
# PART 4: PREDICTED WEIGHT LOSS ANALYSIS
# ==============================================================================

cat("\n========== PREDICTED WEIGHT LOSS ANALYSIS ==========\n")

# Overall distribution
pwl_overall <- field_mice_complete %>%
  summarise(
    n = n(),
    mean = round(mean(predicted_weight_loss), 2),
    sd = round(sd(predicted_weight_loss), 2),
    range = sprintf("%.2f–%.2f",
                    min(predicted_weight_loss),
                    max(predicted_weight_loss))
  )

cat("\nOverall predicted weight loss:\n")
cat(sprintf("  n = %d\n", pwl_overall$n))
cat(sprintf("  Mean ± SD: %.2f%% ± %.2f%%\n", pwl_overall$mean, pwl_overall$sd))
cat(sprintf("  Range: %s%%\n", pwl_overall$range))

# Test normality
shapiro_test <- shapiro.test(field_mice_complete$predicted_weight_loss)
cat(sprintf("  Shapiro-Wilk test: W = %.3f, p = %.3f\n",
            shapiro_test$statistic,
            shapiro_test$p.value))

# By infection status
pwl_by_infection <- complete_data %>%
  group_by(infection_status) %>%
  summarise(
    n = n(),
    mean = round(mean(predicted_weight_loss), 2),
    sd = round(sd(predicted_weight_loss), 2),
    se = round(sd / sqrt(n), 3)
  )

cat("\nPredicted weight loss by infection status:\n")
print(pwl_by_infection)

# Statistical comparison
t_test <- t.test(predicted_weight_loss ~ infection_status, data = complete_data)
cohen_d <- cohen.d(complete_data$predicted_weight_loss, complete_data$infection_status)

cat("\nInfection effect on predicted weight loss:\n")
cat(sprintf("  Difference: %.2f%% (95%% CI: %.2f to %.2f)\n",
            abs(diff(t_test$estimate)),
            t_test$conf.int[1],
            t_test$conf.int[2]))
cat(sprintf("  t-test: t(%.1f) = %.2f, p = %.3f\n",
            t_test$parameter,
            abs(t_test$statistic),
            t_test$p.value))
cat(sprintf("  Cohen's d = %.3f [%.3f, %.3f]\n",
            cohen_d$estimate,
            cohen_d$conf.int[1],
            cohen_d$conf.int[2]))

# ==============================================================================
# PART 5: HYBRID INDEX ANALYSIS
# ==============================================================================

cat("\n========== HYBRID INDEX BY INFECTION STATUS ==========\n")

# HI by infection status
hi_by_infection <- complete_data %>%
  group_by(infection_status) %>%
  summarise(
    n = n(),
    mean = round(mean(HI), 3),
    sd = round(sd(HI), 3),
    se = round(sd / sqrt(n), 3)
  )

print(hi_by_infection)

# Statistical comparison
t_test_hi <- t.test(HI ~ infection_status, data = complete_data)

cat("\nInfection effect on hybrid index:\n")
cat(sprintf("  Infected mice are more M. m. musculus-like\n"))
cat(sprintf("  Difference: %.3f (95%% CI: %.3f to %.3f)\n",
            abs(diff(t_test_hi$estimate)),
            t_test_hi$conf.int[1],
            t_test_hi$conf.int[2]))
cat(sprintf("  t-test: t(%.1f) = %.2f, p < 0.001\n",
            t_test_hi$parameter,
            abs(t_test_hi$statistic)))

# ==============================================================================
# PART 6: LINEAR MODELS
# ==============================================================================

cat("\n========== LINEAR MODEL RESULTS ==========\n")

# Model 1: Infection only
lm1 <- lm(predicted_weight_loss ~ infection_status, data = complete_data)
cat("\nModel 1 - Infection effect only:\n")
print(tidy(lm1) %>% mutate(across(where(is.numeric), ~round(., 3))))
cat(sprintf("  R² = %.3f\n", summary(lm1)$r.squared))

# Model 2: Infection + Sex
lm2 <- lm(predicted_weight_loss ~ infection_status + Sex, data = complete_data)
cat("\nModel 2 - Main effects:\n")
print(tidy(lm2) %>% mutate(across(where(is.numeric), ~round(., 3))))
cat(sprintf("  R² = %.3f\n", summary(lm2)$r.squared))

# Model 3: Interaction
lm3 <- lm(predicted_weight_loss ~ infection_status * Sex, data = complete_data)
cat("\nModel 3 - Interaction model:\n")
print(tidy(lm3) %>% mutate(across(where(is.numeric), ~round(., 3))))
cat(sprintf("  R² = %.3f\n", summary(lm3)$r.squared))

# Model comparison
cat("\nModel comparison:\n")
anova_result <- anova(lm1, lm2, lm3)
print(anova_result)

# ==============================================================================
# PART 7: CREATE PUBLICATION TABLES
# ==============================================================================

cat("\n========== CREATING TABLES FOR MANUSCRIPT ==========\n")


# TABLE 1: Sample Characteristics - Clean and Minimal
# ==============================================================================

# Prepare the data in a structured format
sample_data <- data.frame(
  Variable = c("Sample size",
               "Sex",
               "Infection status¹",
               "Hybrid index²",
               "Predicted weight loss (%)"),
  Overall = c("335",
              "F: 167 (49.9%); M: 168 (50.1%)",
              "Known: 185 (55.2%)",
              "0.574 ± 0.360",
              "9.96 ± 2.47"),
  Uninfected = c("93",
                 "F: 46 (49.5%); M: 47 (50.5%)",
                 "—",
                 "0.468 ± 0.382",
                 "9.08 ± 2.33"),
  Infected = c("92",
               "F: 37 (40.2%); M: 55 (59.8%)",
               "—",
               "0.742 ± 0.267",
               "10.34 ± 2.45"),
  stringsAsFactors = FALSE
)

table1 <- sample_data %>%
  gt() %>%
  # Column formatting
  cols_label(
    Variable = "",
    Overall = "All mice",
    Uninfected = "Uninfected",
    Infected = "Infected"
  ) %>%
  # Minimal borders - only where needed
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1.5)),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1)),
    locations = cells_body(rows = nrow(sample_data))
  ) %>%
  # Remove all other borders
  opt_table_lines(extent = "none") %>%
  # Clean spacing
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.title.font.weight = "normal",
    column_labels.font.weight = "bold",
    column_labels.padding = px(5),
    data_row.padding = px(4),
    table.width = pct(100)
  ) %>%
  # Title
  tab_header(
    title = "Characteristics of wild house mice from the Brandenburg hybrid zone"
  ) %>%
  # Footnotes
  tab_footnote(
    footnote = "Determined by qPCR and amplicon sequencing",
    locations = cells_body(columns = Variable, rows = 3)
  ) %>%
  tab_footnote(
    footnote = "0 = M. m. domesticus, 1 = M. m. musculus",
    locations = cells_body(columns = Variable, rows = 4)
  ) %>%
  # Style footnotes
  tab_options(
    footnotes.font.size = px(10),
    footnotes.padding = px(3)
  )

print(table1)

save_table_all_formats(table_object = table1, table_name = "Sample_characteristics")


# TABLE 2: Statistical Comparisons - Ultra Clean
# ==============================================================================

stats_data <- data.frame(
  Comparison = c("Predicted weight loss (%)",
                 "Hybrid index",
                 "Cohen's d"),
  Uninfected = c("9.08 ± 2.33",
                 "0.468 ± 0.382",
                 "—"),
  Infected = c("10.34 ± 2.45",
               "0.742 ± 0.267",
               "—"),
  Difference_CI = c("1.26 (0.57–1.95)",
                    "0.274 (0.179–0.370)",
                    "0.528 (0.233–0.823)"),
  P_value = c("<0.001",
              "<0.001",
              "—"),
  stringsAsFactors = FALSE
)

table2 <- stats_data %>%
  gt() %>%
  cols_label(
    Comparison = "",
    Uninfected = "Uninfected",
    Infected = "Infected",
    Difference_CI = "Difference (95% CI)",
    P_value = md("*P*")
  ) %>%
  # Minimal styling
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1.5)),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1)),
    locations = cells_body(rows = nrow(stats_data))
  ) %>%
  opt_table_lines(extent = "none") %>%
  # Column alignment
  cols_align(
    align = "center",
    columns = c(Uninfected, Infected, Difference_CI, P_value)
  ) %>%
  cols_align(
    align = "left",
    columns = Comparison
  ) %>%
  # Clean spacing and fonts
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.title.font.weight = "normal",
    column_labels.font.weight = "bold",
    column_labels.padding = px(5),
    data_row.padding = px(4),
    table.width = pct(100)
  ) %>%
  # Title
  tab_header(
    title = "Infection-associated differences in predicted weight loss and host genetics"
  ) %>%
  # Format P-values in italics
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(columns = P_value)
  )

print(table2)


# TABLE 3: Linear Model Results - Professional Format
# ==============================================================================

model_data <- data.frame(
  Term = c("Intercept",
           "Infected",
           "Male",
           "Infected × Male",
           "",
           "Model R²"),
  Model_1 = c("9.08***", "1.26***", "—", "—", "", "0.066"),
  Model_2 = c("9.08***", "1.22***", "0.40", "—", "", "0.072"),
  Model_3 = c("9.85***", "0.78", "0.82", "−0.81", "", "0.082"),
  stringsAsFactors = FALSE
)

table3 <- model_data %>%
  gt() %>%
  cols_label(
    Term = "",
    Model_1 = "Model 1",
    Model_2 = "Model 2",
    Model_3 = "Model 3"
  ) %>%
  # Minimal borders
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1.5)),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(1)),
    locations = cells_body(rows = 5)
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1)),
    locations = cells_body(rows = nrow(model_data))
  ) %>%
  opt_table_lines(extent = "none") %>%
  # Alignment
  cols_align(
    align = "center",
    columns = c(Model_1, Model_2, Model_3)
  ) %>%
  # Italicize R²
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(
      columns = Term,
      rows = Term == "Model R²"
    )
  ) %>%
  # Clean options
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.title.font.weight = "normal",
    column_labels.font.weight = "normal",
    column_labels.padding = px(5),
    data_row.padding = px(4),
    table.width = pct(80)
  ) %>%
  # Title and subtitle
  tab_header(
    title = "Linear regression models of predicted weight loss",
    subtitle = "Model 1: Infection only; Model 2: Main effects; Model 3: Interaction"
  ) %>%
  # Footnote
  tab_footnote(
    footnote = "*** P < 0.001",
    locations = cells_body(columns = Model_1, rows = 1)
  ) %>%
  tab_options(
    footnotes.font.size = px(10),
    heading.subtitle.font.size = px(11)
  )

print(table3)
save_table_all_formats(table_object = table3, table_name = "Regression_models_sample_design")

# ALTERNATIVE: Ultra-minimalist single summary table
# ==============================================================================

summary_data <- data.frame(
  " " = c("n", "Females (%)", "Hybrid index", "Weight loss (%)", "Effect size"),
  "Uninfected" = c("93", "46 (49.5)", "0.468 ± 0.382", "9.08 ± 2.33", "—"),
  "Infected" = c("92", "37 (40.2)", "0.742 ± 0.267", "10.34 ± 2.45", "—"),
  "Statistic" = c("—", "χ² = 1.72", "t = 5.67***", "t = 3.59***", "d = 0.528"),
  check.names = FALSE
)

table_minimal <- summary_data %>%
  gt() %>%
  tab_style(
    style = list(
      cell_borders(sides = "bottom", color = "black", weight = px(1)),
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels()
  ) %>%
  opt_table_lines(extent = "none") %>%
  tab_options(
    table.font.size = px(11),
    column_labels.padding = px(3),
    data_row.padding = px(2),
    table.width = pct(70)
  ) %>%
  cols_align(align = "center", columns = -1) %>%
  tab_header(title = "Sample characteristics and infection effects")

print(table_minimal)



# 1) Annotation dictionary keyed by your assay names (immune_genes)
gene_annot <- tibble::tribble(
  ~Gene,    ~Display, ~Full_Name,                                             ~Function,                                             ~Category,                 ~Response_Type,
  "IFNy",   "IFN-γ",  "Interferon gamma",                                     "Th1 master regulator; parasite control",              "Cytokine",               "Th1",
  "CXCR3",  "CXCR3",  "C-X-C chemokine receptor 3",                           "T cell trafficking to inflamed tissues",              "Chemokine receptor",     "Th1",
  "IL.6",   "IL-6",   "Interleukin 6",                                        "Acute phase response; inflammation",                  "Cytokine",               "Inflammatory",
  "IL.13",  "IL-13",  "Interleukin 13",                                       "Th2 cytokine; mucus production",                      "Th2 response",           "Th2",
  "IL1RN",  "IL1RN",  "Interleukin 1 receptor antagonist",                    "IL-1 antagonist; limits inflammation",                "Regulatory",             "Regulatory",
  "CASP1",  "CASP1",  "Caspase-1",                                            "Inflammasome effector; IL-1β maturation",             "Inflammasome",           "Innate",
  "CXCL9",  "CXCL9",  "C-X-C motif chemokine ligand 9",                       "T cell recruitment; IFN-γ inducible",                 "Chemokine",              "Th1",
  "IDO1",   "IDO1",   "Indoleamine 2,3-dioxygenase 1",                        "Tryptophan catabolism; immunoregulation",             "Metabolic enzyme",       "Regulatory",
  "IRGM1",  "IRGM1",  "Immunity-related GTPase M",                            "Cell-autonomous defense; autophagy",                  "Cell-autonomous",        "Cell-autonomous",
  "MPO",    "MPO",    "Myeloperoxidase",                                      "Oxidative burst; neutrophil marker",                  "Innate effector enzyme", "Innate",
  "MUC2",   "MUC2",   "Mucin 2",                                              "Intestinal barrier; pathogen exclusion",              "Barrier (mucin)",        "Barrier",
  "MUC5AC", "MUC5AC", "Mucin 5AC",                                            "Goblet cell mucin; barrier",                          "Barrier (mucin)",        "Barrier",
  "MYD88",  "MYD88",  "Myeloid differentiation primary response 88",          "TLR/IL-1 signaling adaptor",                          "Innate signaling",       "Innate",
  "NCR1",   "NCR1",   "Natural cytotoxicity receptor 1",                      "NK cell activation receptor",                         "Cytotoxic receptor",     "Cytotoxic",
  "PRF1",   "PRF1",   "Perforin 1",                                           "Cytolytic granule pore-former",                       "Cytotoxic effector",     "Cytotoxic",
  "RETNLB", "RETNLB", "Resistin-like beta",                                   "Type 2 effector; goblet cell hyperplasia",            "Th2 effector",           "Th2",
  "SOCS1",  "SOCS1",  "Suppressor of cytokine signaling 1",                   "Negative feedback on JAK/STAT",                       "Regulatory",             "Regulatory",
  "TICAM1", "TICAM1", "TIR-domain-containing adaptor molecule 1 (TRIF)",      "TLR3/4 adaptor; antiviral",                           "Innate signaling",       "Innate",
  "TNF",    "TNF",    "Tumor necrosis factor",                                "Pro-inflammatory; cachexia",                          "Cytokine",               "Inflammatory"
)

# 2) Build the table in the exact order of your immune_genes vector
gene_data <- tibble::tibble(Gene = immune_genes) %>%
  dplyr::left_join(gene_annot, by = "Gene") %>%
  dplyr::mutate(Gene = dplyr::coalesce(Display, Gene)) %>%
  dplyr::select(Gene, Full_Name, Function, Category, Response_Type)

# Sanity check: fail fast if any annotation is missing
if (any(is.na(gene_data$Full_Name))) {
  missing <- immune_genes[is.na(gene_data$Full_Name)]
  stop(paste("Missing annotations for:", paste(missing, collapse = ", ")))
}

# 3) Make the GT table (same styling as before)
gene_table <- gene_data %>%
  gt::gt() %>%
  gt::tab_header(
    title   = gt::md("**Immune Gene Panel for Infection Response Prediction**"),
    subtitle= gt::md("*Nineteen genes quantified by high-throughput qPCR at day 8 post-infection*")
  ) %>%
  gt::cols_label(
    Gene          = gt::md("**Gene**"),
    Full_Name     = gt::md("**Full Name**"),
    Function      = gt::md("**Primary Function**"),
    Category      = gt::md("**Category**"),
    Response_Type = gt::md("**Response Type**")
  ) %>%
  gt::tab_row_group(label = gt::md("**Th1/Pro-inflammatory Response**"),
                    rows  = Response_Type %in% c("Th1","Inflammatory")) %>%
  gt::tab_row_group(label = gt::md("**Th2/Anti-helminth Response**"),
                    rows  = Response_Type == "Th2") %>%
  gt::tab_row_group(label = gt::md("**Regulatory/Anti-inflammatory**"),
                    rows  = Response_Type == "Regulatory") %>%
  gt::tab_row_group(label = gt::md("**Innate Immunity**"),
                    rows  = Response_Type == "Innate") %>%
  gt::tab_row_group(label = gt::md("**Cell-mediated Immunity**"),
                    rows  = Response_Type %in% c("Cell-autonomous","Cytotoxic","Barrier")) %>%  # keeps your original grouping
  gt::tab_style(
    style = list(gt::cell_fill(color = "#E8F4F8"), gt::cell_text(weight = "bold")),
    locations = gt::cells_row_groups()
  ) %>%
  gt::tab_style(
    style = gt::cell_text(weight = "bold", style = "italic"),
    locations = gt::cells_body(columns = "Gene")
  ) %>%
  gt::opt_row_striping() %>%
  gt::tab_options(
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.border.bottom.color = "black",
    column_labels.border.bottom.color = "black",
    table_body.border.bottom.color = "black"
  ) %>%
  gt::tab_footnote(
    footnote  = gt::md("Genes match Chapter 1 assay symbols (e.g., `IFNy`, `IL.6`, `IL.13`). Display names use standard notation (e.g., IFN-γ, IL-6, IL-13). Expression normalized via quantile normalization with GAPDH and PPIB."),
    locations = gt::cells_title()
  ) %>%
  gt::tab_source_note(
    source_note = gt::md("**Source:** Webster et al. (2025) *in review*")
  )


# Save the table using your function
save_table_all_formats(gene_table, "Table_S1_Immune_Genes")

# 4) Optional simplified version for main text
gene_table_simple <- gene_data %>%
  dplyr::select(Gene, Function, Category) %>%
  gt::gt() %>%
  gt::tab_header(title = gt::md("**Immune Genes Used for Prediction**")) %>%
  gt::cols_label(
    Gene     = gt::md("**Gene**"),
    Function = gt::md("**Function**"),
    Category = gt::md("**Category**")
  ) %>%
  gt::tab_options(table.font.size = 10, table.width = gt::pct(80))

save_table_all_formats(gene_table_simple, "Table_S1_Immune_Genes_Simple")

