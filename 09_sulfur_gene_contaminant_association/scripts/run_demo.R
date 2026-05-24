# Purpose:
#   Build a reproducible sulfur gene / contaminant association demo from shared
#   simulated environmental and functional annotation toy data. The workflow
#   follows the original sulfur script: target sulfur KO extraction, Pearson
#   correlation table, ggcorrplot heatmap, significant scatter regressions, and
#   multiple linear regression summary plots.
#
# Input files:
#   - ../data/toy_shared/environmental_variables.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/sulfur_gene_abundance_by_sample.csv
#   - results/sulfur_gene_contaminant_analysis_table.csv
#   - results/pearson_results.csv
#   - results/mlr_results.csv
#   - figures/fig1_pearson_heatmap.pdf/png
#   - figures/fig2_pearson_scatter.pdf/png
#   - figures/fig3_mlr_results.pdf/png
#
# Required columns:
#   environmental_variables.csv: sample_id plus selected contaminant variables
#   functional_annotation_table.csv: sample_id, ko_id, count, abundance
#
# User-editable settings:
#   Edit the settings block below to change input paths, target KOs, value
#   column, contaminant variables, scatter-pair selection, output names, plot
#   sizes, and palettes.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

environmental_file <- file.path(shared_data_dir, "environmental_variables.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

target_ko_names <- c(
  K11180 = "dsrA",
  K11181 = "dsrB",
  K17218 = "soxB",
  K17222 = "soxC",
  K17223 = "soxD",
  K00958 = "sat"
)

functional_value_column <- "abundance"
contaminant_vars <- c("Sb_III", "Sb_V", "Sb_total", "SO4", "As", "Cu")
primary_response_vars <- c("Sb_III", "Sb_V")
preferred_gene_order <- c("dsrAB", "soxBCD", "sat", "dsrA", "dsrB", "soxB", "soxC", "soxD")
scatter_gene_candidates <- c("dsrAB", "soxBCD", "sat")
max_scatter_pairs <- 4

sulfur_gene_output <- "sulfur_gene_abundance_by_sample.csv"
analysis_table_output <- "sulfur_gene_contaminant_analysis_table.csv"
pearson_output <- "pearson_results.csv"
mlr_output <- "mlr_results.csv"

pearson_heatmap_pdf <- "fig1_pearson_heatmap.pdf"
pearson_heatmap_png <- "fig1_pearson_heatmap.png"
scatter_pdf <- "fig2_pearson_scatter.pdf"
scatter_png <- "fig2_pearson_scatter.png"
mlr_pdf <- "fig3_mlr_results.pdf"
mlr_png <- "fig3_mlr_results.png"

heatmap_width <- 8
heatmap_height <- 7
scatter_width <- 10
scatter_height <- 7
mlr_width <- 12
mlr_height <- 8
png_dpi <- 300

heatmap_colors <- c("#A9D1E8", "white", "#E3A8A8")
scatter_point_low <- "#5C8CAE"
scatter_point_high <- "#B86B6B"
smooth_line_color <- "#9A4848"
smooth_fill_color <- "#E5B5B5"
positive_beta_color <- "#E3A8A8"
negative_beta_color <- "#A9D1E8"

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("dplyr", "tidyr", "tibble", "ggplot2", "ggcorrplot", "ggrepel", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggcorrplot)
  library(ggrepel)
  library(patchwork)
})

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(environmental_file, functional_annotation_file))

environmental_variables <- read.csv(environmental_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(environmental_variables, c("sample_id", contaminant_vars), "environmental_variables.csv")
check_required_columns(functional_annotation, c("sample_id", "ko_id", functional_value_column), "functional_annotation_table.csv")
check_sample_ids(environmental_variables, functional_annotation)

sulfur_gene_matrix <- make_sulfur_gene_matrix(
  functional_annotation = functional_annotation,
  target_ko_names = target_ko_names,
  value_column = functional_value_column
)
write.csv(sulfur_gene_matrix, file.path(results_dir, sulfur_gene_output), row.names = FALSE)

analysis_table <- merge(sulfur_gene_matrix, environmental_variables, by = "sample_id", sort = FALSE)
write.csv(analysis_table, file.path(results_dir, analysis_table_output), row.names = FALSE)

gene_vars <- intersect(preferred_gene_order, colnames(analysis_table))
gene_vars <- gene_vars[vapply(gene_vars, function(x) any(!is.na(analysis_table[[x]])), logical(1))]
if (length(gene_vars) < 2) {
  stop("At least two sulfur gene variables are required for this demo.", call. = FALSE)
}

pearson_table <- make_pearson_table(analysis_table, gene_vars, contaminant_vars)
write.csv(pearson_table, file.path(results_dir, pearson_output), row.names = FALSE)

cor_data <- analysis_table[, c(gene_vars, contaminant_vars), drop = FALSE]
cor_data <- cor_data[stats::complete.cases(cor_data), , drop = FALSE]
cor_matrix <- stats::cor(cor_data, method = "pearson")
p_matrix <- cor_pmat(cor_data, method = "pearson")

pearson_heatmap <- ggcorrplot::ggcorrplot(
  cor_matrix,
  method = "square",
  type = "full",
  lab = TRUE,
  lab_size = 3.2,
  p.mat = p_matrix,
  sig.level = 0.05,
  insig = "blank",
  colors = heatmap_colors,
  title = "Pearson Correlation: Sulfur Genes x Contaminants",
  ggtheme = theme_minimal(base_size = 12)
) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black")
  )

save_ggplot_pair(
  pearson_heatmap,
  file.path(figures_dir, pearson_heatmap_pdf),
  file.path(figures_dir, pearson_heatmap_png),
  heatmap_width,
  heatmap_height,
  png_dpi
)

scatter_source <- pearson_table[pearson_table$Gene %in% intersect(scatter_gene_candidates, gene_vars), , drop = FALSE]
if (nrow(scatter_source) == 0) scatter_source <- pearson_table
scatter_pairs <- select_scatter_pairs(scatter_source, max_scatter_pairs)
scatter_plots <- list()
for (i in seq_len(nrow(scatter_pairs))) {
  gene <- scatter_pairs$Gene[i]
  contaminant <- scatter_pairs$Contaminant[i]
  subset_data <- analysis_table[stats::complete.cases(analysis_table[, c(gene, contaminant), drop = FALSE]), , drop = FALSE]
  r_value <- stats::cor(subset_data[[gene]], subset_data[[contaminant]], method = "pearson")
  p_value <- stats::cor.test(subset_data[[gene]], subset_data[[contaminant]], method = "pearson")$p.value
  label <- paste0("r = ", round(r_value, 3), "  ", p_to_stars(p_value))

  scatter_plots[[paste(gene, contaminant, sep = "_")]] <- ggplot(subset_data, aes(x = .data[[contaminant]], y = .data[[gene]])) +
    geom_point(aes(color = .data[[contaminant]]), size = 4, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = smooth_line_color, fill = smooth_fill_color, linewidth = 1.1) +
    scale_color_gradient(low = scatter_point_low, high = scatter_point_high) +
    annotate("text", x = -Inf, y = Inf, label = label, hjust = -0.1, vjust = 1.6, size = 4.2, fontface = "italic") +
    labs(
      x = contaminant,
      y = paste(gene, "relative abundance"),
      title = paste(gene, "vs", contaminant)
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

scatter_plot <- patchwork::wrap_plots(scatter_plots, ncol = 2)
save_ggplot_pair(
  scatter_plot,
  file.path(figures_dir, scatter_pdf),
  file.path(figures_dir, scatter_png),
  scatter_width,
  scatter_height,
  png_dpi
)

mlr_predictors <- intersect(c("dsrAB", "soxBCD", "sat"), gene_vars)
if (length(mlr_predictors) < 2) {
  mlr_predictors <- head(gene_vars, 2)
}
mlr_results <- fit_mlr_models(analysis_table, mlr_predictors, primary_response_vars)
write.csv(mlr_results$table, file.path(results_dir, mlr_output), row.names = FALSE)

mlr_plots <- list()
for (response in primary_response_vars) {
  fit_data <- mlr_results$plot_data[[paste0(response, "_fit")]]
  r2_label <- sprintf(
    "R2 = %.3f  (adj.R2 = %.3f)\nModel p = %.4f",
    unique(fit_data$R2),
    unique(fit_data$R2_adj),
    unique(fit_data$model_p)
  )
  mlr_plots[[paste0(response, "_fit")]] <- ggplot(fit_data, aes(x = fitted, y = actual)) +
    geom_point(color = "#8FBBA8", size = 3.5, alpha = 0.85) +
    geom_abline(linetype = "dashed", color = "grey60", linewidth = 0.8) +
    geom_smooth(method = "lm", color = "#D98C8C", fill = "#F0CACA", se = TRUE, linewidth = 1) +
    annotate("text", x = -Inf, y = Inf, label = r2_label, hjust = -0.05, vjust = 1.4, size = 3.5) +
    labs(
      x = paste0("Fitted ", response),
      y = paste0("Measured ", response),
      title = paste("MLR fitted:", response)
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  beta_data <- mlr_results$plot_data[[paste0(response, "_beta")]]
  beta_data$sig <- p_to_stars(beta_data$p)
  beta_data$label <- paste0(round(beta_data$beta, 3), " ", beta_data$sig)
  beta_data$label_y <- beta_data$beta + ifelse(beta_data$beta >= 0, 0.035, -0.035)
  beta_axis_limits <- range(c(0, beta_data$lower, beta_data$upper, beta_data$label_y), na.rm = TRUE)
  beta_axis_padding <- diff(beta_axis_limits) * 0.12
  if (beta_axis_padding == 0) beta_axis_padding <- 0.1

  mlr_plots[[paste0(response, "_beta")]] <- ggplot(beta_data, aes(x = Gene, y = beta, fill = beta > 0)) +
    geom_col(width = 0.5) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 0.8) +
    geom_text(aes(label = label, y = label_y, hjust = ifelse(beta >= 0, 0, 1)), size = 3.6) +
    geom_hline(yintercept = 0, linewidth = 0.7, color = "grey40") +
    scale_fill_manual(values = c(`TRUE` = positive_beta_color, `FALSE` = negative_beta_color), guide = "none") +
    scale_y_continuous(limits = c(beta_axis_limits[1] - beta_axis_padding, beta_axis_limits[2] + beta_axis_padding)) +
    coord_flip() +
    labs(x = NULL, y = "Standardized beta (95% CI)", title = paste("MLR coefficients ->", response)) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}

mlr_plot <- patchwork::wrap_plots(mlr_plots, ncol = 2)
save_ggplot_pair(
  mlr_plot,
  file.path(figures_dir, mlr_pdf),
  file.path(figures_dir, mlr_png),
  mlr_width,
  mlr_height,
  png_dpi
)

message("Sulfur gene / contaminant association demo completed.")
message("Sulfur variables: ", paste(gene_vars, collapse = ", "))
message("Contaminants: ", paste(contaminant_vars, collapse = ", "))
message("MLR predictors: ", paste(mlr_predictors, collapse = ", "))
message("Outputs written to results/ and figures/.")
