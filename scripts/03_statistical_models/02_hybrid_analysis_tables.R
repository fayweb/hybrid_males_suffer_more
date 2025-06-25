# ==============================================================================
# TABLE WITH COLORS MATCHING YOUR FIGURE
# ==============================================================================

create_formatted_parasiteload_table_with_figure_colors <- function(results_table) {

  formatted_table <- results_table %>%
    dplyr::select(Hypothesis, Complete_result, Constitutive_result, Infection_result) %>%
    gt() %>%

    # Column headers with EXPLICIT group definitions
    cols_label(
      Hypothesis = "Statistical Test",
      Complete_result = html("Complete Dataset<br><em>(n = 304)</em><br><strong>A = Female, B = Male</strong>"),
      Constitutive_result = html("Constitutive Costs<br><em>(n = 171, uninfected)</em><br><strong>A = Female, B = Male</strong>"),
      Infection_result = html("Infection Dominance<br><em>(n = 304, by status)</em><br><strong>A = Uninfected, B = Infected</strong>")
    ) %>%

    # Table header
    tab_header(
      title = "",
      subtitle = html("")
    ) %>%

    # Group rows by test type
    tab_row_group(
      label = "Individual Hypothesis Tests",
      rows = 1:6
    ) %>%
    tab_row_group(
      label = "Model Comparisons",
      rows = 7:10
    ) %>%

    # COLOR SCHEME MATCHING YOUR FIGURE
    # H0 rows - Orange/Peach color (outermost ring)
    tab_style(
      style = cell_fill(color = "#FFD4A3", alpha = 0.7),  # Light orange/peach
      locations = cells_body(rows = c(1))  # H0: Overall hybrid effect
    ) %>%

    # H1 rows - Pink/Purple color (second ring)
    tab_style(
      style = cell_fill(color = "#E6B8E6", alpha = 0.7),  # Light purple/pink
      locations = cells_body(rows = c(2, 7))  # H1: Group difference + H1 vs H0 comparison
    ) %>%

    # H2 rows - Light blue color (third ring)
    tab_style(
      style = cell_fill(color = "#B8D4E6", alpha = 0.7),  # Light blue
      locations = cells_body(rows = c(3, 4, 8))  # H2: Group A, Group B + H2 vs H0 comparison
    ) %>%

    # H3 rows - Green color (innermost ring)
    tab_style(
      style = cell_fill(color = "#C8E6C8", alpha = 0.7),  # Light green
      locations = cells_body(rows = c(5, 6, 9, 10))  # H3: Advanced models + H3 comparisons
    ) %>%

    # CORRECTED BOLDING - Only bold cells that are actually significant
    # Complete Dataset column - bold only significant p-values
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = Complete_result,
        rows = which(results_table$Complete_p < 0.05)
      )
    ) %>%

    # Constitutive column - bold only significant p-values
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = Constitutive_result,
        rows = which(results_table$Constitutive_p < 0.05)
      )
    ) %>%

    # Infection column - bold only significant p-values
    tab_style(
      style = cell_text(weight = "bold"),
      locations = cells_body(
        columns = Infection_result,
        rows = which(results_table$Infection_p < 0.05)
      )
    ) %>%

    # Enhanced footnotes with group clarification
    tab_footnote(
      footnote = "α represents deviation from additivity (positive = increased costs with M. m. musculus ancestry)",
      locations = cells_body(columns = Hypothesis, rows = 1)
    ) %>%
    tab_footnote(
      footnote = "Group A = Female (Complete/Constitutive) OR Uninfected (Infection analysis)",
      locations = cells_body(columns = Hypothesis, rows = 3)
    ) %>%
    tab_footnote(
      footnote = "Group B = Male (Complete/Constitutive) OR Infected (Infection analysis) - Males show significant effect (p = 0.038*)",
      locations = cells_body(columns = Hypothesis, rows = 4)
    ) %>%
    tab_footnote(
      footnote = "Critical comparison for detecting group-specific effects (H3 vs H2)",
      locations = cells_body(columns = Hypothesis, rows = 10)
    ) %>%
    tab_footnote(
      footnote = "Complete & Constitutive: A = Female, B = Male | Infection: A = Uninfected, B = Infected with Eimeria spp."
    ) %>%
    tab_footnote(
      footnote = "Significance: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1 | Bold = p < 0.05"
    ) %>%
    tab_footnote(
      footnote = html("Color scheme matches nested hypothesis structure: <span style='background-color: #FFD4A3; padding: 2px;'>H0</span> <span style='background-color: #E6B8E6; padding: 2px;'>H1</span> <span style='background-color: #B8D4E6; padding: 2px;'>H2</span> <span style='background-color: #C8E6C8; padding: 2px;'>H3</span>")
    ) %>%

    # Table options
    tab_options(
      table.font.size = 11,
      heading.title.font.size = 14,
      heading.subtitle.font.size = 12,
      row_group.font.weight = "bold",
      table.width = pct(100)
    )

  return(formatted_table)
}

# ==============================================================================
# CREATE TABLE WITH FIGURE-MATCHING COLORS
# ==============================================================================

cat("Creating table with colors matching your nested hypothesis figure...\n")

# Create final table with figure-matching colors
final_parasiteload_table_figure_colors <- create_formatted_parasiteload_table_with_figure_colors(parasiteload_table_final)

# Save table with figure colors
save_table_all_formats(
  table_object = final_parasiteload_table_figure_colors,
  table_name = "parasiteLoad_results_FINAL_figure_colors",
  output_dir = "results/tables"
)

cat("✅ Table with figure-matching colors created and saved!\n")

# ==============================================================================
# COLOR MAPPING VERIFICATION
# ==============================================================================

cat("\n=== COLOR MAPPING VERIFICATION ===\n")
cat("===================================\n")

color_mapping <- data.frame(
  Hypothesis_Level = c("H0", "H1", "H2", "H3"),
  Color_Name = c("Orange/Peach", "Purple/Pink", "Light Blue", "Light Green"),
  Hex_Code = c("#FFD4A3", "#E6B8E6", "#B8D4E6", "#C8E6C8"),
  Rows_Affected = c(
    "Row 1 (H0: Overall hybrid effect)",
    "Rows 2, 7 (H1: Group difference, H1 vs H0)",
    "Rows 3, 4, 8 (H2: Group A, Group B, H2 vs H0)",
    "Rows 5, 6, 9, 10 (H3: Advanced models, H3 comparisons)"
  ),
  stringsAsFactors = FALSE
)

print(color_mapping)

cat("\n🎨 Color scheme now matches your nested hypothesis figure perfectly!\n")
cat("Each hypothesis level has its corresponding color from the concentric circles.\n")
