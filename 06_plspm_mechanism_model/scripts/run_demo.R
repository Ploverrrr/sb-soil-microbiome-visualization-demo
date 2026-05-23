# Purpose:
#   Build a reproducible PLS-PM mechanism-model demo from shared simulated
#   soil chemistry, microbial abundance, and functional annotation toy data.
#   The analysis follows the original plspm workflow: plspm(), innerplot(),
#   outerplot(loadings), outerplot(weights), and exported model tables.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/environmental_variables.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/plspm_model_input_table.csv
#   - results/plspm_scaled_input_table.csv
#   - results/plspm_path_matrix.csv
#   - results/plspm_latent_variable_blocks.csv
#   - results/plspm_path_coefficients.csv
#   - results/plspm_inner_model.csv
#   - results/plspm_inner_summary.csv
#   - results/plspm_outer_model.csv
#   - results/plspm_effects.csv
#   - results/plspm_latent_scores.csv
#   - results/plspm_model_metrics.csv
#   - figures/plspm_inner_path_model.pdf/png
#   - figures/plspm_outer_loadings.pdf/png
#   - figures/plspm_outer_weights.pdf/png
#   - figures/plspm_total_effects.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   environmental_variables.csv: sample_id plus the selected environmental,
#     nutrient, contaminant, and metal indicators
#   abundance_table.csv: feature_id plus one column per sample_id
#   functional_annotation_table.csv: sample_id, function_category, count
#
# User-editable settings:
#   Edit the settings block below to change paths, latent-variable indicators,
#   path structure, plspm modes, colors, and output file names.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
environmental_file <- file.path(shared_data_dir, "environmental_variables.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

# These blocks follow the original PLS-PM scripts' mechanism narrative:
# Environment -> Nutrient -> Sb_As -> Function -> Alpha, with direct paths.
latent_variable_blocks <- list(
  Environment = c("pH", "EC"),
  Nutrient = c("TOC", "TN", "TP"),
  Sb_As = c("Sb_total", "As"),
  Function = c("MetalRes", "Sulfur", "Nitrogen", "Carbon"),
  Alpha = c("Chao1", "Shannon", "Simpson", "Pielou")
)

Environment <- c(0, 0, 0, 0, 0)
Nutrient <- c(1, 0, 0, 0, 0)
Sb_As <- c(1, 1, 0, 0, 0)
Function <- c(1, 1, 1, 0, 0)
Alpha <- c(1, 1, 1, 1, 0)
path_matrix <- rbind(Environment, Nutrient, Sb_As, Function, Alpha)
colnames(path_matrix) <- rownames(path_matrix)

plspm_modes <- rep("A", length(latent_variable_blocks))

positive_path_color <- "red"
negative_path_color <- "blue"
link_color <- "gray"
box_line_width <- 0
arrow_width <- 0.1
show_path_values <- TRUE

model_input_output <- "plspm_model_input_table.csv"
scaled_input_output <- "plspm_scaled_input_table.csv"
path_matrix_output <- "plspm_path_matrix.csv"
block_table_output <- "plspm_latent_variable_blocks.csv"
path_coefficients_output <- "plspm_path_coefficients.csv"
inner_model_output <- "plspm_inner_model.csv"
inner_summary_output <- "plspm_inner_summary.csv"
outer_model_output <- "plspm_outer_model.csv"
effects_output <- "plspm_effects.csv"
latent_scores_output <- "plspm_latent_scores.csv"
model_metrics_output <- "plspm_model_metrics.csv"

inner_path_pdf <- "plspm_inner_path_model.pdf"
inner_path_png <- "plspm_inner_path_model.png"
outer_loadings_pdf <- "plspm_outer_loadings.pdf"
outer_loadings_png <- "plspm_outer_loadings.png"
outer_weights_pdf <- "plspm_outer_weights.pdf"
outer_weights_png <- "plspm_outer_weights.png"
effects_pdf <- "plspm_total_effects.pdf"
effects_png <- "plspm_total_effects.png"

inner_plot_width <- 7
inner_plot_height <- 5.6
outer_plot_width <- 9
outer_plot_height <- 6.2
effects_plot_width <- 8
effects_plot_height <- 5.8
png_dpi <- 300

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("plspm")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(metadata_file, environmental_file, abundance_file, functional_annotation_file))

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
environmental_variables <- read.csv(environmental_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_table <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", "group"), "sample_metadata.csv")
check_required_columns(environmental_variables, c("sample_id", unlist(latent_variable_blocks[c("Environment", "Nutrient", "Sb_As")])), "environmental_variables.csv")
check_required_columns(abundance_table, "feature_id", "abundance_table.csv")
check_required_columns(functional_annotation, c("sample_id", "function_category", "count"), "functional_annotation_table.csv")
check_sample_ids(sample_metadata, environmental_variables, abundance_table, functional_annotation)

if (!identical(names(latent_variable_blocks), rownames(path_matrix))) {
  stop("latent_variable_blocks names must match path_matrix row names.", call. = FALSE)
}
if (!identical(rownames(path_matrix), colnames(path_matrix))) {
  stop("path_matrix row names and column names must match.", call. = FALSE)
}
if (length(plspm_modes) != length(latent_variable_blocks)) {
  stop("plspm_modes length must match the number of latent-variable blocks.", call. = FALSE)
}

model_input <- build_plspm_input(sample_metadata, environmental_variables, abundance_table, functional_annotation)
indicator_columns <- unique(unlist(latent_variable_blocks))
check_required_columns(model_input, c("sample_id", "group", indicator_columns), "derived PLS-PM model input")
scaled_indicators <- scale_indicator_table(model_input, indicator_columns)

write.csv(model_input, file.path(results_dir, model_input_output), row.names = FALSE)
write.csv(cbind(model_input[, c("sample_id", "group")], scaled_indicators), file.path(results_dir, scaled_input_output), row.names = FALSE)
write.csv(path_matrix, file.path(results_dir, path_matrix_output), row.names = TRUE)
write.csv(make_block_table(latent_variable_blocks), file.path(results_dir, block_table_output), row.names = FALSE)

plspm_model <- plspm::plspm(
  Data = scaled_indicators,
  path_matrix = path_matrix,
  blocks = latent_variable_blocks,
  modes = plspm_modes
)

path_coefficients <- as.data.frame(plspm_model$path_coefs)
path_coefficients$target_latent_variable <- rownames(path_coefficients)
path_coefficients <- path_coefficients[, c("target_latent_variable", setdiff(colnames(path_coefficients), "target_latent_variable")), drop = FALSE]
write.csv(path_coefficients, file.path(results_dir, path_coefficients_output), row.names = FALSE)

inner_model_table <- flatten_inner_model(plspm_model$inner_model)
write.csv(inner_model_table, file.path(results_dir, inner_model_output), row.names = FALSE)
write.csv(plspm_model$inner_summary, file.path(results_dir, inner_summary_output), row.names = FALSE)
write.csv(plspm_model$outer_model, file.path(results_dir, outer_model_output), row.names = FALSE)
write.csv(plspm_model$effects, file.path(results_dir, effects_output), row.names = FALSE)
write.csv(cbind(sample_id = model_input$sample_id, group = model_input$group, as.data.frame(plspm_model$scores)), file.path(results_dir, latent_scores_output), row.names = FALSE)

model_metrics <- data.frame(
  metric = c("gof", "sample_count", "manifest_variable_count", "latent_variable_count"),
  value = c(plspm_model$gof, nrow(model_input), length(indicator_columns), length(latent_variable_blocks)),
  stringsAsFactors = FALSE
)
write.csv(model_metrics, file.path(results_dir, model_metrics_output), row.names = FALSE)

save_plspm_base_plot(
  file.path(figures_dir, inner_path_pdf),
  file.path(figures_dir, inner_path_png),
  inner_plot_width,
  inner_plot_height,
  png_dpi,
  function() {
    plspm::innerplot(
      plspm_model,
      colpos = positive_path_color,
      colneg = negative_path_color,
      show.values = show_path_values,
      lcol = link_color,
      box.lwd = box_line_width
    )
  }
)

save_plspm_base_plot(
  file.path(figures_dir, outer_loadings_pdf),
  file.path(figures_dir, outer_loadings_png),
  outer_plot_width,
  outer_plot_height,
  png_dpi,
  function() {
    plspm::outerplot(
      plspm_model,
      what = "loadings",
      arr.width = arrow_width,
      colpos = positive_path_color,
      colneg = negative_path_color,
      show.values = show_path_values,
      lcol = link_color
    )
  }
)

save_plspm_base_plot(
  file.path(figures_dir, outer_weights_pdf),
  file.path(figures_dir, outer_weights_png),
  outer_plot_width,
  outer_plot_height,
  png_dpi,
  function() {
    plspm::outerplot(
      plspm_model,
      what = "weights",
      arr.width = arrow_width,
      colpos = positive_path_color,
      colneg = negative_path_color,
      show.values = show_path_values,
      lcol = link_color
    )
  }
)

save_effects_plot(
  plspm_model$effects,
  file.path(figures_dir, effects_pdf),
  file.path(figures_dir, effects_png),
  effects_plot_width,
  effects_plot_height,
  png_dpi,
  positive_path_color,
  negative_path_color
)

message("PLS-PM mechanism model demo completed.")
message("Latent variables: ", paste(names(latent_variable_blocks), collapse = " -> "))
message("Goodness-of-fit: ", round(plspm_model$gof, 3))
message("Outputs written to results/ and figures/.")
