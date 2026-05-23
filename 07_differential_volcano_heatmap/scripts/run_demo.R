# Purpose:
#   Build a reproducible differential abundance volcano + heatmap demo from
#   shared simulated soil microbiome taxonomy and abundance toy data. The
#   plotting workflow follows the original scripts: DESeq2 differential
#   analysis, ggplot2 volcano plots with ggrepel labels, patchwork combining,
#   and ComplexHeatmap/circlize heatmaps.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/taxonomy_table.csv
#
# Output files:
#   - results/taxon_count_matrix.csv
#   - results/deseq2_all_contrasts.csv
#   - results/deseq2_<group>_vs_<control>.csv
#   - results/heatmap_zscore_matrix.csv
#   - results/heatmap_selected_taxa.csv
#   - figures/volcano_<group>_vs_<control>.pdf/png
#   - figures/deseq2_volcano_combined.pdf/png
#   - figures/differential_heatmap.pdf/png
#   - figures/differential_circular_heatmap.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   abundance_table.csv: feature_id plus one column per sample_id
#   taxonomy_table.csv: feature_id plus target_taxonomic_level
#
# User-editable settings:
#   Edit the settings block below to change input paths, taxonomic level,
#   control/treatment groups, thresholds, labels, heatmap size, colors, and
#   output file names.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")

group_column <- "group"
control_group <- "Control"
treatment_groups <- c("Tailing", "Mining", "Smelting")
target_taxonomic_level <- "Genus"

deseq_fit_type <- "parametric"
log2fc_threshold <- 1
pvalue_threshold <- 0.05
top_label_n <- 3
volcano_y_axis_cap <- 80
heatmap_focus_group <- "Smelting"
heatmap_top_n_each_direction <- 20

taxon_count_output <- "taxon_count_matrix.csv"
all_contrasts_output <- "deseq2_all_contrasts.csv"
heatmap_matrix_output <- "heatmap_zscore_matrix.csv"
heatmap_selected_taxa_output <- "heatmap_selected_taxa.csv"

combined_volcano_pdf <- "deseq2_volcano_combined.pdf"
combined_volcano_png <- "deseq2_volcano_combined.png"
rectangular_heatmap_pdf <- "differential_heatmap.pdf"
rectangular_heatmap_png <- "differential_heatmap.png"
circular_heatmap_pdf <- "differential_circular_heatmap.pdf"
circular_heatmap_png <- "differential_circular_heatmap.png"

volcano_width <- 5
volcano_height <- 5.5
combined_volcano_width <- 15
combined_volcano_height <- 5.5
rectangular_heatmap_width <- 7
rectangular_heatmap_height <- 7
circular_heatmap_width <- 7
circular_heatmap_height <- 7
png_dpi <- 300

# Original scripts used blue/cyan for down-regulated taxa, grey for stable
# taxa, and red/pink for up-regulated taxa.
volcano_colors <- c(Down = "#00c9c8", Stable = "grey", Up = "#ed5e93")
heatmap_colors <- c("#00c9c8", "white", "#ed5e93")
heatmap_breaks <- c(-2, 0, 3.5)
row_name_color <- "grey"

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("DESeq2", "ggplot2", "ggrepel", "patchwork", "ComplexHeatmap", "circlize")
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
  library(ggrepel)
  library(patchwork)
})

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(metadata_file, abundance_file, taxonomy_file))

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_table <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
taxonomy_table <- read.csv(taxonomy_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", group_column), "sample_metadata.csv")
check_required_columns(abundance_table, "feature_id", "abundance_table.csv")
check_required_columns(taxonomy_table, c("feature_id", target_taxonomic_level), "taxonomy_table.csv")
check_sample_ids(sample_metadata, abundance_table)
check_feature_ids(abundance_table, taxonomy_table)

missing_groups <- setdiff(c(control_group, treatment_groups), unique(sample_metadata[[group_column]]))
if (length(missing_groups) > 0) {
  stop("Selected group value(s) not found in metadata: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}
if (!heatmap_focus_group %in% treatment_groups) {
  stop("heatmap_focus_group must be one of treatment_groups.", call. = FALSE)
}
if (top_label_n < 1) stop("top_label_n must be at least 1.", call. = FALSE)
if (heatmap_top_n_each_direction < 1) stop("heatmap_top_n_each_direction must be at least 1.", call. = FALSE)

taxon_counts <- aggregate_counts_to_taxon(abundance_table, taxonomy_table, target_taxonomic_level)
write.csv(taxon_counts, file.path(results_dir, taxon_count_output), row.names = FALSE)

volcano_plots <- vector("list", length(treatment_groups))
names(volcano_plots) <- treatment_groups
contrast_results <- vector("list", length(treatment_groups))

for (treatment_group in treatment_groups) {
  diff_result <- run_deseq2_contrast(
    count_table = taxon_counts,
    sample_metadata = sample_metadata,
    control_group = control_group,
    treatment_group = treatment_group,
    group_column = group_column,
    fit_type = deseq_fit_type
  )
  diff_result <- classify_change(diff_result, log2fc_threshold, pvalue_threshold)
  diff_result$comparison <- paste(treatment_group, "vs", control_group)
  diff_result$gene <- diff_result$taxon
  diff_result$negative_log10_pvalue <- -log10(pmax(diff_result$pvalue, .Machine$double.xmin))
  diff_result$negative_log10_pvalue_for_plot <- pmin(diff_result$negative_log10_pvalue, volcano_y_axis_cap)
  contrast_results[[treatment_group]] <- diff_result

  output_name <- paste0("deseq2_", tolower(treatment_group), "_vs_", tolower(control_group), ".csv")
  write.csv(diff_result, file.path(results_dir, output_name), row.names = FALSE)

  top_labeled_taxa <- select_top_labeled_taxa(diff_result, top_label_n)

  volcano_plot <- ggplot(diff_result, aes(x = log2FoldChange, y = negative_log10_pvalue_for_plot)) +
    geom_point(aes(color = change), alpha = 1, size = 1.5) +
    scale_color_manual(values = volcano_colors, breaks = c("Down", "Stable", "Up")) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), lty = 4, col = "black", lwd = 0.6) +
    geom_hline(yintercept = -log10(pvalue_threshold), lty = 4, col = "black", lwd = 0.6) +
    labs(
      title = paste(treatment_group, "vs", control_group),
      x = "Log2 (Fold Change)",
      y = "-Log10 (Pvalue)",
      color = "Change"
    ) +
    scale_y_continuous(limits = c(0, volcano_y_axis_cap), expand = expansion(mult = c(0.02, 0.08))) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10),
      text = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    geom_text_repel(
      data = top_labeled_taxa,
      aes(label = gene),
      vjust = -0.5,
      size = 4,
      fontface = "bold",
      box.padding = 0.45,
      point.padding = 0.2,
      max.overlaps = Inf
    )

  volcano_plots[[treatment_group]] <- volcano_plot
  safe_name <- paste0("volcano_", tolower(treatment_group), "_vs_", tolower(control_group))
  save_ggplot_pair(
    volcano_plot,
    file.path(figures_dir, paste0(safe_name, ".pdf")),
    file.path(figures_dir, paste0(safe_name, ".png")),
    volcano_width,
    volcano_height,
    png_dpi
  )
}

all_contrasts <- do.call(rbind, contrast_results)
write.csv(all_contrasts, file.path(results_dir, all_contrasts_output), row.names = FALSE)

combined_volcano <- volcano_plots[[1]] + volcano_plots[[2]] + volcano_plots[[3]]
save_ggplot_pair(
  combined_volcano,
  file.path(figures_dir, combined_volcano_pdf),
  file.path(figures_dir, combined_volcano_png),
  combined_volcano_width,
  combined_volcano_height,
  png_dpi
)

focus_result <- contrast_results[[heatmap_focus_group]]
selected_taxa <- select_heatmap_taxa(focus_result, heatmap_top_n_each_direction)
group_order <- c(control_group, treatment_groups)
sample_order <- sample_metadata$sample_id[order(match(sample_metadata[[group_column]], group_order), sample_metadata$sample_id)]
if ("replicate" %in% colnames(sample_metadata)) {
  sample_order <- sample_metadata$sample_id[order(match(sample_metadata[[group_column]], group_order), sample_metadata$replicate)]
}
heatmap_matrix <- make_zscore_matrix(taxon_counts, selected_taxa, sample_order)
write.csv(data.frame(taxon = rownames(heatmap_matrix), heatmap_matrix, check.names = FALSE), file.path(results_dir, heatmap_matrix_output), row.names = FALSE)
write.csv(focus_result[focus_result$taxon %in% selected_taxa, , drop = FALSE], file.path(results_dir, heatmap_selected_taxa_output), row.names = FALSE)

sample_annotation <- data.frame(
  Group = factor(sample_metadata[[group_column]][match(colnames(heatmap_matrix), sample_metadata$sample_id)], levels = group_order)
)
rownames(sample_annotation) <- colnames(heatmap_matrix)
annotation_colors <- list(
  Group = c(
    Control = "#4db6ac",
    Tailing = "#f6c85f",
    Mining = "#e45756",
    Smelting = "#7b6fd6"
  )
)
color_function <- circlize::colorRamp2(heatmap_breaks, heatmap_colors)

column_annotation <- ComplexHeatmap::HeatmapAnnotation(
  df = sample_annotation,
  col = annotation_colors,
  annotation_name_gp = grid::gpar(fontsize = 9)
)

rectangular_heatmap <- ComplexHeatmap::Heatmap(
  heatmap_matrix,
  name = "Expression",
  col = color_function,
  top_annotation = column_annotation,
  row_names_gp = grid::gpar(fontsize = 5),
  column_names_gp = grid::gpar(fontsize = 7),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_split = sample_annotation$Group
)

save_complex_heatmap_pair(
  rectangular_heatmap,
  file.path(figures_dir, rectangular_heatmap_pdf),
  file.path(figures_dir, rectangular_heatmap_png),
  rectangular_heatmap_width,
  rectangular_heatmap_height,
  png_dpi
)

save_circular_heatmap_pair(
  heatmap_matrix,
  file.path(figures_dir, circular_heatmap_pdf),
  file.path(figures_dir, circular_heatmap_png),
  circular_heatmap_width,
  circular_heatmap_height,
  png_dpi,
  color_function,
  row_name_color
)

message("Differential abundance volcano + heatmap demo completed.")
message("Taxonomic level: ", target_taxonomic_level)
message("Contrasts: ", paste(paste(treatment_groups, control_group, sep = " vs "), collapse = "; "))
message("Heatmap focus group: ", heatmap_focus_group)
message("Outputs written to results/ and figures/.")
