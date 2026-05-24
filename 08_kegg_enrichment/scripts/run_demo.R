# Purpose:
#   Build a reproducible KEGG-style enrichment demo from shared simulated
#   functional annotation toy data. The plotting workflow follows the original
#   KEGG enrichment scripts: clusterProfiler enrichment results, barplot(),
#   dotplot(), patchwork combination, and an auxiliary ggplot bubble plot.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/ko_count_matrix.csv
#   - results/differential_ko_statistics.csv
#   - results/foreground_ko_list.csv
#   - results/kegg_pathway_enrichment_result.csv
#   - results/kegg_module_enrichment_result.csv
#   - figures/kegg_pathway_enrichment_combined.pdf/png
#   - figures/kegg_module_enrichment_combined.pdf/png
#   - figures/kegg_pathway_bubble_plot.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   functional_annotation_table.csv: sample_id, ko_id, pathway, count
#
# User-editable settings:
#   Edit the settings block below to change input paths, contrast groups,
#   enrichment backend, thresholds, p-value adjustment, output names, plot
#   sizes, and colors.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

group_column <- "group"
control_group <- "Control"
treatment_group <- "Smelting"
ko_id_pattern <- "^K"

# Backend options:
#   "clusterprofiler_kegg" - original method, uses enrichKEGG() and enrichMKEGG()
#   "toy_offline"         - reproducible demo mode using enricher() + toy TERM2GENE
#   "auto"                - try the original method first, then fall back to toy_offline
# For real data analysis, use "clusterprofiler_kegg" when KEGG access is available.
# You can also override without editing this file:
#   KEGG_ENRICHMENT_BACKEND=clusterprofiler_kegg Rscript scripts/run_demo.R
enrichment_backend <- Sys.getenv("KEGG_ENRICHMENT_BACKEND", unset = "toy_offline")

deseq_fit_type <- "parametric"
foreground_log2fc_threshold <- 0.5
foreground_pvalue_threshold <- 0.1
minimum_foreground_kos <- 5

enrichment_pvalue_cutoff <- 1
p_adjust_method <- "BH"
qvalue_cutoff <- 1
min_gs_size <- 1
max_gs_size <- 500
show_category <- 10
enrichment_plot_color_by <- "pvalue"

ko_count_output <- "ko_count_matrix.csv"
differential_output <- "differential_ko_statistics.csv"
foreground_output <- "foreground_ko_list.csv"
kegg_pathway_output <- "kegg_pathway_enrichment_result.csv"
kegg_module_output <- "kegg_module_enrichment_result.csv"

kegg_pathway_combined_pdf <- "kegg_pathway_enrichment_combined.pdf"
kegg_pathway_combined_png <- "kegg_pathway_enrichment_combined.png"
kegg_module_combined_pdf <- "kegg_module_enrichment_combined.pdf"
kegg_module_combined_png <- "kegg_module_enrichment_combined.png"
kegg_bubble_pdf <- "kegg_pathway_bubble_plot.pdf"
kegg_bubble_png <- "kegg_pathway_bubble_plot.png"

combined_width <- 14.5
combined_height <- 6
bubble_width <- 7
bubble_height <- 5.8
png_dpi <- 300

bubble_color_palette <- c("#ed5e93", "#00c9c8")
bubble_size_range <- c(2, 8)

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("DESeq2", "clusterProfiler", "enrichplot", "ggplot2", "patchwork")
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
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
})

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(metadata_file, functional_annotation_file))

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", group_column), "sample_metadata.csv")
check_required_columns(functional_annotation, c("sample_id", "ko_id", "pathway", "count"), "functional_annotation_table.csv")
check_sample_ids(sample_metadata, functional_annotation)

missing_groups <- setdiff(c(control_group, treatment_group), unique(sample_metadata[[group_column]]))
if (length(missing_groups) > 0) {
  stop("Selected group value(s) not found in metadata: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}
if (!enrichment_backend %in% c("toy_offline", "clusterprofiler_kegg", "auto")) {
  stop("enrichment_backend must be 'toy_offline', 'clusterprofiler_kegg', or 'auto'.", call. = FALSE)
}
if (!enrichment_plot_color_by %in% c("pvalue", "p.adjust", "qvalue")) {
  stop("enrichment_plot_color_by must be 'pvalue', 'p.adjust', or 'qvalue'.", call. = FALSE)
}

ko_counts <- aggregate_ko_counts(functional_annotation, ko_id_pattern)
write.csv(ko_counts, file.path(results_dir, ko_count_output), row.names = FALSE)

differential_ko <- run_deseq2_ko_contrast(
  ko_count_table = ko_counts,
  sample_metadata = sample_metadata,
  control_group = control_group,
  treatment_group = treatment_group,
  group_column = group_column,
  fit_type = deseq_fit_type
)
differential_ko <- classify_differential_ko(differential_ko, foreground_log2fc_threshold, foreground_pvalue_threshold)
write.csv(differential_ko, file.path(results_dir, differential_output), row.names = FALSE)

foreground_ko <- select_foreground_kos(
  differential_ko,
  min_ko_count = minimum_foreground_kos,
  log2fc_threshold = foreground_log2fc_threshold,
  pvalue_threshold = foreground_pvalue_threshold
)
foreground_ko_list <- foreground_ko$ko_id
universe_ko_list <- unique(ko_counts$ko_id)
write.csv(foreground_ko, file.path(results_dir, foreground_output), row.names = FALSE)

kegg_results <- run_enrichment_backend(
  enrichment_backend = enrichment_backend,
  ko_list = foreground_ko_list,
  universe = universe_ko_list,
  pvalue_cutoff = enrichment_pvalue_cutoff,
  p_adjust_method = p_adjust_method,
  qvalue_cutoff = qvalue_cutoff,
  min_gs_size = min_gs_size,
  max_gs_size = max_gs_size
)
kegg_pathway_result <- kegg_results$ko
kegg_module_result <- kegg_results$module

kegg_pathway_table <- data.frame(kegg_pathway_result)
kegg_module_table <- data.frame(kegg_module_result)
write.csv(kegg_pathway_table, file.path(results_dir, kegg_pathway_output), row.names = FALSE)
write.csv(kegg_module_table, file.path(results_dir, kegg_module_output), row.names = FALSE)

if (nrow(kegg_pathway_table) == 0 || nrow(kegg_module_table) == 0) {
  stop(
    "Enrichment returned an empty result. For the small toy dataset, keep enrichment_backend = 'toy_offline' ",
    "and enrichment_pvalue_cutoff/qvalue_cutoff at 1, or adjust the foreground thresholds.",
    call. = FALSE
  )
}

ko_plots <- make_enrichment_plots(kegg_pathway_result, show_category, enrichment_plot_color_by)
module_plots <- make_enrichment_plots(kegg_module_result, show_category, enrichment_plot_color_by)

kegg_pathway_combined <- ko_plots$bar + ko_plots$dot
kegg_module_combined <- module_plots$bar + module_plots$dot

save_ggplot_pair(
  kegg_pathway_combined,
  file.path(figures_dir, kegg_pathway_combined_pdf),
  file.path(figures_dir, kegg_pathway_combined_png),
  combined_width,
  combined_height,
  png_dpi
)
save_ggplot_pair(
  kegg_module_combined,
  file.path(figures_dir, kegg_module_combined_pdf),
  file.path(figures_dir, kegg_module_combined_png),
  combined_width,
  combined_height,
  png_dpi
)

bubble_table <- kegg_pathway_table[order(kegg_pathway_table$p.adjust, -kegg_pathway_table$Count), , drop = FALSE]
bubble_table <- head(bubble_table, show_category)
bubble_plot <- make_enrichment_bubble_plot(
  bubble_table,
  title = "KEGG-style enrichment bubble plot",
  color_palette = bubble_color_palette,
  size_range = bubble_size_range
)
save_ggplot_pair(
  bubble_plot,
  file.path(figures_dir, kegg_bubble_pdf),
  file.path(figures_dir, kegg_bubble_png),
  bubble_width,
  bubble_height,
  png_dpi
)

message("KEGG enrichment demo completed.")
message("Backend: ", enrichment_backend)
message("Contrast: ", treatment_group, " vs ", control_group)
message("Foreground KOs: ", length(foreground_ko_list), " of ", length(universe_ko_list), " KO universe")
message("KEGG pathway terms: ", nrow(kegg_pathway_table))
message("KEGG module-like terms: ", nrow(kegg_module_table))
message("Outputs written to results/ and figures/.")
