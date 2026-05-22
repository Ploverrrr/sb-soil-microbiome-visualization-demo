# Purpose:
#   Build a reproducible random forest importance plus environmental
#   correlation heatmap demo from the shared simulated soil microbiome toy data.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/environmental_variables.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/taxonomy_table.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/relative_abundance_by_feature.csv
#   - results/selected_feature_abundance.csv
#   - results/correlation_results.csv
#   - results/rf_importance.csv
#   - results/rf_model_performance.csv
#   - results/heatmap_plotting_table.csv
#   - figures/correlation_heatmap.pdf
#   - figures/correlation_heatmap.png
#   - figures/rf_importance_plot.pdf
#   - figures/rf_importance_plot.png
#   - figures/rf_correlation_combined.pdf
#   - figures/rf_correlation_combined.png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   environmental_variables.csv: sample_id and selected environmental variables
#   abundance_table.csv: feature_id, one column per sample_id
#   taxonomy_table.csv: feature_id and selected target_taxonomic_level
#   functional_annotation_table.csv: sample_id, pathway, count when feature_source = "function"
#
# User-editable settings:
#   Edit the settings block below to change input paths, feature source,
#   taxonomy level, environmental variables, RF target, top_n, output names,
#   figure size, color palette, and bubble size range.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
environment_file <- file.path(shared_data_dir, "environmental_variables.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

feature_source <- "taxonomy" # options: "taxonomy" or "function"
target_taxonomic_level <- "Genus"
functional_feature_column <- "pathway"

environmental_variables <- c(
  "pH", "EC", "TOC", "TN", "TP", "SO4", "NO3",
  "Sb_total", "Sb_III", "Sb_V", "As", "Cu", "Zn", "Cd", "Fe", "Mn"
)

target_environmental_variable_for_rf <- "Sb_total"

top_n_features <- 12
correlation_method <- "spearman"
p_adjust_method <- "BH"

relative_abundance_output <- "relative_abundance_by_feature.csv"
selected_feature_output <- "selected_feature_abundance.csv"
correlation_output <- "correlation_results.csv"
rf_importance_output <- "rf_importance.csv"
rf_performance_output <- "rf_model_performance.csv"
heatmap_plotting_output <- "heatmap_plotting_table.csv"

correlation_heatmap_pdf <- "correlation_heatmap.pdf"
correlation_heatmap_png <- "correlation_heatmap.png"
rf_importance_pdf <- "rf_importance_plot.pdf"
rf_importance_png <- "rf_importance_plot.png"
combined_figure_pdf <- "rf_correlation_combined.pdf"
combined_figure_png <- "rf_correlation_combined.png"

heatmap_width <- 9.5
heatmap_height <- 6.6
rf_plot_width <- 6.6
rf_plot_height <- 4.8
combined_figure_width <- 12.2
combined_figure_height <- 6.6
png_dpi <- 300

heatmap_color_palette <- c("#00c9c8", "white", "#ed5e93")
bubble_size_range <- c(1.8, 8.2)

rf_ntree <- 500
rf_seed <- 20260522

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("ggplot2", "randomForest", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(patchwork))

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

input_files <- c(metadata_file, environment_file, abundance_file, taxonomy_file, functional_annotation_file)
check_files_exist(input_files)

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
environmental_data <- read.csv(environment_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_table <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
taxonomy_table <- read.csv(taxonomy_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", "group"), "sample_metadata.csv")
check_required_columns(environmental_data, c("sample_id", environmental_variables), "environmental_variables.csv")
check_required_columns(environmental_data, target_environmental_variable_for_rf, "environmental_variables.csv")
check_required_columns(abundance_table, "feature_id", "abundance_table.csv")
check_required_columns(taxonomy_table, c("feature_id", target_taxonomic_level), "taxonomy_table.csv")

check_sample_alignment(sample_metadata, environmental_data, abundance_table)
check_feature_alignment(abundance_table, taxonomy_table)

if (!feature_source %in% c("taxonomy", "function")) {
  stop("feature_source must be either 'taxonomy' or 'function'.", call. = FALSE)
}

if (!top_n_features >= 2) {
  stop("top_n_features must be at least 2.", call. = FALSE)
}

sample_columns <- setdiff(colnames(abundance_table), "feature_id")
sample_metadata <- sample_metadata[match(sample_columns, sample_metadata$sample_id), , drop = FALSE]
environmental_data <- environmental_data[match(sample_columns, environmental_data$sample_id), , drop = FALSE]

relative_abundance <- calculate_relative_abundance(abundance_table)
write.csv(relative_abundance, file.path(results_dir, relative_abundance_output), row.names = FALSE)

if (feature_source == "taxonomy") {
  feature_abundance <- aggregate_taxonomy_abundance(
    relative_abundance = relative_abundance,
    taxonomy = taxonomy_table,
    target_taxonomic_level = target_taxonomic_level
  )
  feature_label <- target_taxonomic_level
} else {
  feature_abundance <- aggregate_function_abundance(
    functional_annotation = functional_annotation,
    feature_column = functional_feature_column
  )
  feature_label <- functional_feature_column
}

available_feature_count <- nrow(feature_abundance)
if (top_n_features > available_feature_count) {
  stop(
    "top_n_features is larger than the number of available features: ",
    available_feature_count,
    call. = FALSE
  )
}

selected_features <- select_top_features(feature_abundance, top_n_features)
write.csv(selected_features, file.path(results_dir, selected_feature_output), row.names = FALSE)

correlation_results <- make_correlation_table(
  feature_abundance = selected_features,
  environmental_data = environmental_data,
  environmental_variables = environmental_variables,
  correlation_method = correlation_method,
  p_adjust_method = p_adjust_method
)
write.csv(correlation_results, file.path(results_dir, correlation_output), row.names = FALSE)

selected_sample_columns <- setdiff(colnames(selected_features), c("feature", "overall_mean_abundance"))
rf_predictors <- t(as.matrix(selected_features[, selected_sample_columns, drop = FALSE]))
colnames(rf_predictors) <- make.names(selected_features$feature, unique = TRUE)
rf_response <- environmental_data[[target_environmental_variable_for_rf]]

set.seed(rf_seed)
rf_model <- randomForest(
  x = as.data.frame(rf_predictors),
  y = rf_response,
  ntree = rf_ntree,
  importance = TRUE
)

importance_matrix <- importance(rf_model)
importance_column <- if ("%IncMSE" %in% colnames(importance_matrix)) "%IncMSE" else colnames(importance_matrix)[1]
rf_importance <- data.frame(
  feature = selected_features$feature,
  rf_predictor_name = colnames(rf_predictors),
  importance = as.numeric(importance_matrix[, importance_column]),
  importance_metric = importance_column,
  stringsAsFactors = FALSE
)
rf_importance$importance_for_size <- pmax(rf_importance$importance, 0)
rf_importance$bubble_size <- scale_to_range(rf_importance$importance_for_size, bubble_size_range)
rf_importance <- rf_importance[order(rf_importance$importance, decreasing = TRUE), , drop = FALSE]
write.csv(rf_importance, file.path(results_dir, rf_importance_output), row.names = FALSE)

rf_performance <- data.frame(
  target_environmental_variable = target_environmental_variable_for_rf,
  ntree = rf_ntree,
  final_oob_mse = tail(rf_model$mse, 1),
  final_oob_pseudo_r2 = tail(rf_model$rsq, 1),
  stringsAsFactors = FALSE
)
write.csv(rf_performance, file.path(results_dir, rf_performance_output), row.names = FALSE)

heatmap_plotting_table <- merge(correlation_results, rf_importance[, c("feature", "importance", "bubble_size")], by = "feature", all.x = TRUE)
heatmap_plotting_table$bubble_environmental_variable <- heatmap_plotting_table$environmental_variable == target_environmental_variable_for_rf
write.csv(heatmap_plotting_table, file.path(results_dir, heatmap_plotting_output), row.names = FALSE)

feature_order <- selected_features$feature
env_order <- environmental_variables
heatmap_plotting_table$feature <- factor(heatmap_plotting_table$feature, levels = rev(feature_order))
heatmap_plotting_table$environmental_variable <- factor(heatmap_plotting_table$environmental_variable, levels = env_order)

heatmap_plot <- ggplot(heatmap_plotting_table, aes(x = environmental_variable, y = feature)) +
  geom_tile(aes(fill = correlation), color = "grey88", linewidth = 0.25) +
  geom_point(
    data = heatmap_plotting_table[heatmap_plotting_table$bubble_environmental_variable, ],
    aes(size = bubble_size),
    shape = 21,
    stroke = 0.85,
    color = "black",
    fill = NA
  ) +
  geom_text(aes(label = significance), size = 3.4, color = "black", na.rm = TRUE) +
  scale_fill_gradient2(
    low = heatmap_color_palette[1],
    mid = heatmap_color_palette[2],
    high = heatmap_color_palette[3],
    midpoint = 0,
    limits = c(-1, 1),
    name = paste0(tools::toTitleCase(correlation_method), " r")
  ) +
  scale_size_identity(guide = "none") +
  labs(
    title = "Feature-environment correlation heatmap",
    subtitle = paste0("RF importance bubbles shown for ", target_environmental_variable_for_rf),
    x = NULL,
    y = feature_label
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black", face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

rf_importance$feature <- factor(rf_importance$feature, levels = rev(rf_importance$feature))
rf_plot <- ggplot(rf_importance, aes(x = importance, y = feature)) +
  geom_col(fill = "#009292", width = 0.72) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.4) +
  labs(
    title = "Random forest feature importance",
    subtitle = paste0(
      "Target: ", target_environmental_variable_for_rf,
      " | OOB pseudo R2 = ", round(rf_performance$final_oob_pseudo_r2, 3)
    ),
    x = importance_column,
    y = feature_label
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black", face = "italic"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

rf_importance_combined <- rf_importance
rf_importance_combined$feature <- factor(rf_importance_combined$feature, levels = levels(heatmap_plotting_table$feature))
rf_panel_combined <- ggplot(rf_importance_combined, aes(x = importance, y = feature)) +
  geom_col(fill = "#009292", width = 0.68) +
  geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
  labs(
    title = "RF importance",
    subtitle = paste0("OOB pseudo R2 = ", round(rf_performance$final_oob_pseudo_r2, 3)),
    x = importance_column,
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.margin = ggplot2::margin(t = 5.5, r = 5.5, b = 5.5, l = 0)
  )

heatmap_panel_combined <- heatmap_plot +
  labs(
    title = "Correlation heatmap",
    subtitle = paste0("Top ", nrow(selected_features), " ", feature_label, " features"),
    y = feature_label
  ) +
  theme(
    plot.margin = ggplot2::margin(t = 5.5, r = 0, b = 5.5, l = 5.5),
    legend.position = "right"
  )

combined_plot <- heatmap_panel_combined + rf_panel_combined +
  plot_layout(widths = c(4.7, 1.35), guides = "collect") +
  plot_annotation(
    title = "Random forest importance and environmental correlation",
    subtitle = paste0(
      "Feature source: ", feature_source,
      " | RF target: ", target_environmental_variable_for_rf
    ),
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11)
    )
  ) &
  theme(legend.position = "right")

ggsave(file.path(figures_dir, correlation_heatmap_pdf), heatmap_plot, width = heatmap_width, height = heatmap_height, units = "in")
ggsave(file.path(figures_dir, correlation_heatmap_png), heatmap_plot, width = heatmap_width, height = heatmap_height, units = "in", dpi = png_dpi)
ggsave(file.path(figures_dir, rf_importance_pdf), rf_plot, width = rf_plot_width, height = rf_plot_height, units = "in")
ggsave(file.path(figures_dir, rf_importance_png), rf_plot, width = rf_plot_width, height = rf_plot_height, units = "in", dpi = png_dpi)
ggsave(file.path(figures_dir, combined_figure_pdf), combined_plot, width = combined_figure_width, height = combined_figure_height, units = "in")
ggsave(file.path(figures_dir, combined_figure_png), combined_plot, width = combined_figure_width, height = combined_figure_height, units = "in", dpi = png_dpi)

message("RF correlation heatmap demo completed.")
message("Feature source: ", feature_source)
message("Selected features: ", nrow(selected_features))
message("RF target: ", target_environmental_variable_for_rf)
message("OOB pseudo R2: ", round(rf_performance$final_oob_pseudo_r2, 4))
message("Outputs written to results/ and figures/.")
