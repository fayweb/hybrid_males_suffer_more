#By convention, the variable contribution plot has a circle around the variables that has a radius of 1. Here’s some code to make one.
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}


##Combining correlogram with the significance test
## Computing the p-value of correlations
## To compute the matrix of p-value, a custom R function is used
# ... : further arguments to pass to the native R cor.test function
#tutorial: http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
# mat : is a matrix of data

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}


# create a function that is the opposite of intersect
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}


# testing distributions
# Define function to be used to test, get the log lik and aic
tryDistrib <- function(x, distrib){
  # deals with fitdistr error:
  fit <-
    tryCatch(MASS::fitdistr(x, distrib), error=function(err) "fit failed")
  return(list(fit = fit,
              loglik = tryCatch(fit$loglik, error=function(err) "no loglik computed"),
              AIC = tryCatch(fit$aic, error=function(err) "no aic computed")))
}


findGoodDist <- function(x, distribs, distribs2){
  l =lapply(distribs, function(i) tryDistrib(x, i))
  names(l) <- distribs
  print(l)
  listDistr <- lapply(distribs2, function(i){
    if (i %in% "t"){
      fitdistrplus::fitdist(x, i, start = list(df =2))
    } else {
      fitdistrplus::fitdist(x,i)
    }}
  )
  par(mfrow=c(2,2))
  denscomp(listDistr, legendtext=distribs2)
  cdfcomp(listDistr, legendtext=distribs2)
  qqcomp(listDistr, legendtext=distribs2)
  ppcomp(listDistr, legendtext=distribs2)
  par(mfrow=c(1,1))
}

# Function to italicize y-axis and legend labels
italics_y <- function(ggplot_object, labels) {

  ggplot_object +
    scale_y_discrete(labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
      axis.text.y = element_markdown(),  # Apply markdown to y-axis text
      legend.text = element_markdown()   # Apply markdown to legend text
    )
}


# Function to italicize x-axis and legend labels
italics_x <- function(ggplot_object, labels) {

  ggplot_object +
    scale_x_discrete(labels = labels) +
    scale_fill_manual(values = color_mapping, labels = labels) +
    theme(
      axis.text.x = element_markdown(),  # Apply markdown to y-axis text
      legend.text = element_markdown()   # Apply markdown to legend text
    )
}

save_table_all_formats <- function(table_object, table_name, output_dir = "results/tables") {

  # Create table-specific folder
  table_folder <- file.path(output_dir, table_name)
  dir.create(table_folder, recursive = TRUE, showWarnings = FALSE)

  # Define base file path
  base_path <- file.path(table_folder, table_name)

  # Save in multiple formats
  tryCatch({
    # HTML
    gtsave(table_object, filename = paste0(base_path, ".html"))
    cat("✓ Saved", table_name, "as HTML\n")

    # DOCX
    gtsave(table_object, filename = paste0(base_path, ".docx"))
    cat("✓ Saved", table_name, "as DOCX\n")

    # PNG
    gtsave(table_object, filename = paste0(base_path, ".png"),
           vwidth = 1200, vheight = 800)
    cat("✓ Saved", table_name, "as PNG\n")

    # PDF
    gtsave(table_object, filename = paste0(base_path, ".pdf"))
    cat("✓ Saved", table_name, "as PDF\n")

    # TEX
    latex_code <- as_latex(table_object)
    writeLines(latex_code, paste0(base_path, ".tex"))
    cat("✓ Saved", table_name, "as TEX\n")

    cat("✅ All formats saved in folder:", table_folder, "\n\n")

  }, error = function(e) {
    cat("❌ Error saving", table_name, ":", e$message, "\n")
  })
}



save_plot_all_formats <- function(plot_object, plot_name, output_dir = "results/figures",
                                  width = 7, height = 5, dpi = 300) {

  # Create figure-specific folder
  plot_folder <- file.path(output_dir, plot_name)
  dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

  # Define base file path (no extension yet)
  base_path <- file.path(plot_folder, plot_name)

  # Try saving all formats
  tryCatch({
    # PDF (vector graphic, ideal for publications)
    ggsave(filename = paste0(base_path, ".pdf"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in", device = cairo_pdf)
    cat("✓ Saved", plot_name, "as PDF\n")

    # JPEG (raster graphic, high-res)
    ggsave(filename = paste0(base_path, ".jpeg"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in")
    cat("✓ Saved", plot_name, "as JPEG\n")

    cat("✅ All formats saved in folder:", plot_folder, "\n\n")

  }, error = function(e) {
    cat("❌ Error saving", plot_name, ":", e$message, "\n")
  })
}

save_plot_all_formats_tight <- function(plot_object, plot_name, output_dir = "results/figures",
                                  width = 5, height = 5, dpi = 300) {

  # Create figure-specific folder
  plot_folder <- file.path(output_dir, plot_name)
  dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

  # Define base file path (no extension yet)
  base_path <- file.path(plot_folder, plot_name)

  # Try saving all formats
  tryCatch({
    # PDF (vector graphic, ideal for publications)
    ggsave(filename = paste0(base_path, ".pdf"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in", device = cairo_pdf)
    cat("✓ Saved", plot_name, "as PDF\n")

    # JPEG (raster graphic, high-res)
    ggsave(filename = paste0(base_path, ".jpeg"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in")
    cat("✓ Saved", plot_name, "as JPEG\n")

    cat("✅ All formats saved in folder:", plot_folder, "\n\n")

  }, error = function(e) {
    cat("❌ Error saving", plot_name, ":", e$message, "\n")
  })
}

save_plot_all_formats_wide <- function(plot_object, plot_name, output_dir = "results/figures",
                                        width = 11, height = 5, dpi = 300) {

  # Create figure-specific folder
  plot_folder <- file.path(output_dir, plot_name)
  dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

  # Define base file path (no extension yet)
  base_path <- file.path(plot_folder, plot_name)

  # Try saving all formats
  tryCatch({
    # PDF (vector graphic, ideal for publications)
    ggsave(filename = paste0(base_path, ".pdf"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in", device = cairo_pdf)
    cat("✓ Saved", plot_name, "as PDF\n")

    # JPEG (raster graphic, high-res)
    ggsave(filename = paste0(base_path, ".jpeg"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in")
    cat("✓ Saved", plot_name, "as JPEG\n")

    cat("✅ All formats saved in folder:", plot_folder, "\n\n")

  }, error = function(e) {
    cat("❌ Error saving", plot_name, ":", e$message, "\n")
  })
}

save_plot_all_formats_panel <- function(plot_object, plot_name, output_dir = "results/figures",
                                  width = 7, height = 7, dpi = 300) {

  # Create figure-specific folder
  plot_folder <- file.path(output_dir, plot_name)
  dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

  # Define base file path (no extension yet)
  base_path <- file.path(plot_folder, plot_name)

  # Try saving all formats
  tryCatch({
    # PDF (vector graphic, ideal for publications)
    ggsave(filename = paste0(base_path, ".pdf"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in", device = cairo_pdf)
    cat("✓ Saved", plot_name, "as PDF\n")

    # JPEG (raster graphic, high-res)
    ggsave(filename = paste0(base_path, ".jpeg"),
           plot = plot_object, width = width, height = height, dpi = dpi, units = "in")
    cat("✓ Saved", plot_name, "as JPEG\n")

    cat("✅ All formats saved in folder:", plot_folder, "\n\n")

  }, error = function(e) {
    cat("❌ Error saving", plot_name, ":", e$message, "\n")
  })
}


# Function to print section headers
print_section <- function(title) {
  cat("\n", rep("=", 60), "\n")
  cat(toupper(title), "\n")
  cat(rep("=", 60), "\n\n")
}

# Function to save plots with consistent formatting

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

# Function to load primary dataset
load_field_data <- function() {
  cat("Loading field mice dataset...\n")

  field_data <- read_csv(
    file.path("data", "processed", "field_with_predictions.csv"),
    show_col_types = FALSE
  )

  cat("✓ Loaded", nrow(field_data), "wild mice with complete data\n")
  cat("✓ Variables:", ncol(field_data), "columns\n")

  # Basic data validation
  if ("predicted_weight_loss" %in% names(field_data)) {
    cat("✓ Predicted weight loss values present\n")
  } else {
    warning("⚠ Predicted weight loss values missing!")
  }

  return(field_data)
}

# Function to load Chapter 1 Random Forest model
load_chapter1_model <- function() {
  cat("Loading Chapter 1 Random Forest model...\n")

  model_path <- file.path("data", "processed", "chapter1_rf_model.rds")

  if (file.exists(model_path)) {
    model <- readRDS(model_path)
    cat("✓ Random Forest model loaded successfully\n")
    cat("  - Model type:", class(model)[1], "\n")
    cat("  - Number of trees:", model$ntree, "\n")
    cat("  - Variables used:", length(model$forest$xlevels), "\n")
    return(model)
  } else {
    warning("⚠ Chapter 1 model not found at:", model_path)
    return(NULL)
  }
}
