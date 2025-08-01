# ==============================================================================
# MASTER MODEL DESIGN TABLE - CLEAN WHITE VERSION
# ==============================================================================

create_actual_model_table <- function() {
    
    actual_models <- data.frame(
        Model_ID = c("RF-FIELD", "PL-1", "PL-2", "PL-3", "LM-1", "LM-2", "LM-3"),
        
        Model_Type = c("Random Forest", "parasiteLoad", "parasiteLoad", "parasiteLoad", 
                       "Linear Regression", "Linear Regression", "Linear Regression"),
        
        Formula = c("Predicted_WL ~ 19 immune genes (field application)",
                    "response ~ HI + hHe + α×hHe | Sex", 
                    "response ~ HI + hHe + α×hHe | Sex (uninfected)",
                    "response ~ HI + hHe + α×hHe | Infection_status",
                    "Predicted_WL ~ HI × hHe",
                    "Predicted_WL ~ Sex + HI + hHe + Infection", 
                    "Predicted_WL ~ Sex × HI + interactions"),
        
        Purpose = c("Generate immune-based health predictions",
                    "Test sex-specific hybrid dysfunction", 
                    "Test constitutive hybrid costs",
                    "Test infection dominance over genetics",
                    "Test genetic baseline effects",
                    "Test main effects model",
                    "Test interaction effects"),
        
        Sample_Size = c("n = 335", "n = 335", "n = 171", "n = 304", "n = 335", "n = 304", "n = 304"),
        
        Key_Result = c("Range: 3.98-18.5% predicted weight loss",
                       "Sex-specific effects: p = 0.027",
                       "No constitutive costs: p = 0.544", 
                       "Infection dominance: p = 0.002",
                       "Non-significant: p = 0.44",
                       "Infection effect: p < 0.001",
                       "Sex × HI: p = 0.011"),
        
        Analysis_Phase = c("Field Application", 
                           "Hybrid Analysis", "Hybrid Analysis", "Hybrid Analysis",
                           "Linear Regression", "Linear Regression", "Linear Regression")
    )
    
    return(actual_models)
}

# Create the table
master_table <- create_actual_model_table()

# ==============================================================================
# CLEAN WHITE TABLE - NO TITLES, NO COLORS
# ==============================================================================

create_clean_master_model_table <- function(master_table) {
    
    # Create clean formatted table
    pub_table <- master_table %>%
        gt() %>%
        
        # NO COLORS - all white background (default)
        
        # Bold significant results (p < 0.05) only
        tab_style(
            style = cell_text(weight = "bold"),
            locations = cells_body(
                columns = Key_Result,
                rows = c(2, 4, 6, 7)  # PL-1, PL-3, LM-2, LM-3 (significant results)
            )
        ) %>%
        
        # Format columns
        cols_width(
            Analysis_Phase ~ px(130),
            Model_ID ~ px(100),
            Model_Type ~ px(130),
            Formula ~ px(350),
            Purpose ~ px(250),
            Sample_Size ~ px(80),
            Key_Result ~ px(200)
        ) %>%
        
        # Add footnotes
        tab_footnote(
            footnote = "Random Forest model from Webster et al. (in review)",
            locations = cells_body(columns = Model_Type, rows = 1)
        ) %>%
        tab_footnote(
            footnote = "parasiteLoad package implements Baird et al. (2012) framework", 
            locations = cells_body(columns = Model_Type, rows = 2:4)
        ) %>%
        tab_footnote(
            footnote = "α represents deviation from additivity (positive = increased costs with M. m. musculus ancestry)",
            locations = cells_body(columns = Formula, rows = 2:4)
        ) %>%
        tab_footnote(
            footnote = "Bold = p < 0.05"
        ) %>%
        
        # Style headers
        tab_style(
            style = cell_text(weight = "bold"),
            locations = cells_column_labels()
        ) %>%
        
        # Column labels
        cols_label(
            Analysis_Phase = "Analysis Phase",
            Model_ID = "Model ID",
            Model_Type = "Statistical Method",
            Formula = "Model Formula",
            Purpose = "Analytical Purpose",
            Sample_Size = "Sample Size",
            Key_Result = "Key Result"
        ) %>%
        
        # Table options
        tab_options(
            table.font.size = 11,
            table.width = pct(100)
        )
    
    return(pub_table)
}

# Create clean publication table
publication_table <- create_clean_master_model_table(master_table)

# Display the table
publication_table

# Save in all formats
save_table_all_formats(publication_table, "master_model_design_table_clean")

cat("✅ Clean master model table created - all white, no titles!\n")
cat("📝 Ready for journal submission with separate table legend\n")
