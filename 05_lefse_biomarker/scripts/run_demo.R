# Purpose:
#   Build a reproducible LEfSe biomarker demo from shared simulated soil
#   microbiome taxonomy and abundance toy data, reusing the original
#   microeco LEfSe plotting workflow.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/taxonomy_table.csv
#
# Output files:
#   - results/taxon_relative_abundance_by_sample.csv
#   - results/lefse_candidate_statistics.csv
#   - results/lefse_biomarker_table.csv
#   - results/lefse_barplot_plotting_table.csv
#   - results/kw_abundance_plotting_table.csv
#   - results/biomarker_group_heatmap_table.csv
#   - results/cladogram_node_table.csv
#   - results/cladogram_edge_table.csv
#   - results/cladogram_label_table.csv
#   - figures/lefse_biomarker_barplot.pdf/png
#   - figures/lefse_cladogram.pdf/png
#   - figures/lefse_barplot_cladogram_combined.pdf/png
#   - figures/lefse_kw_abundance_plot.pdf/png
#   - figures/lefse_lda_kw_combined.pdf/png
#   - figures/biomarker_group_heatmap.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   abundance_table.csv: feature_id plus one column per sample_id
#   taxonomy_table.csv: feature_id, Kingdom, Phylum, Class, Order, Family, Genus
#
# User-editable settings:
#   Edit the settings block below to change input paths, LEfSe parameters,
#   taxonomic levels, group order, plot feature counts, output names, figure
#   sizes, and palettes.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")

target_taxonomic_level <- "all"
taxonomic_levels_for_lefse <- c("Phylum", "Class", "Order", "Family", "Genus")
group_column <- "group"
group_order <- c("Control", "Tailing", "Mining", "Smelting")

lefse_alpha <- 0.05
p_adjust_method <- "fdr"
lefse_norm <- 1000000
lda_threshold <- 4
top_n_taxa <- 80
plot_feature_count <- 40
cladogram_use_taxa_num <- 300
cladogram_use_feature_num <- 40
clade_label_level <- 4
cladogram_filter_taxa <- 0.0001

taxon_relative_output <- "taxon_relative_abundance_by_sample.csv"
candidate_statistics_output <- "lefse_candidate_statistics.csv"
biomarker_table_output <- "lefse_biomarker_table.csv"
barplot_table_output <- "lefse_barplot_plotting_table.csv"
kw_abundance_table_output <- "kw_abundance_plotting_table.csv"
heatmap_table_output <- "biomarker_group_heatmap_table.csv"
cladogram_node_output <- "cladogram_node_table.csv"
cladogram_edge_output <- "cladogram_edge_table.csv"
cladogram_label_output <- "cladogram_label_table.csv"

barplot_pdf <- "lefse_biomarker_barplot.pdf"
barplot_png <- "lefse_biomarker_barplot.png"
cladogram_pdf <- "lefse_cladogram.pdf"
cladogram_png <- "lefse_cladogram.png"
barplot_cladogram_pdf <- "lefse_barplot_cladogram_combined.pdf"
barplot_cladogram_png <- "lefse_barplot_cladogram_combined.png"
kw_abundance_pdf <- "lefse_kw_abundance_plot.pdf"
kw_abundance_png <- "lefse_kw_abundance_plot.png"
lda_kw_pdf <- "lefse_lda_kw_combined.pdf"
lda_kw_png <- "lefse_lda_kw_combined.png"
heatmap_pdf <- "biomarker_group_heatmap.pdf"
heatmap_png <- "biomarker_group_heatmap.png"

barplot_width <- 5.6
barplot_height <- 6.2
cladogram_width <- 10
cladogram_height <- 6.5
barplot_cladogram_width <- 14
barplot_cladogram_height <- 6.5
kw_abundance_width <- 7.2
kw_abundance_height <- 8
lda_kw_width <- 10.5
lda_kw_height <- 6.5
heatmap_width <- 7.2
heatmap_height <- 6.8
png_dpi <- 300

group_colors <- c(
  Control = "#4db6ac",
  Tailing = "#f6c85f",
  Mining = "#e45756",
  Smelting = "#7b6fd6"
)

heatmap_palette <- c("#f7fbff", "#9bd4e4", "#f05a8a")

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("ggplot2", "microeco", "aplot", "gridExtra")
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
check_required_columns(
  taxonomy_table,
  c("feature_id", "Kingdom", taxonomic_levels_for_lefse),
  "taxonomy_table.csv"
)
check_sample_ids(sample_metadata, abundance_table)
check_feature_ids(abundance_table, taxonomy_table)

if (anyDuplicated(sample_metadata$sample_id) > 0) {
  stop("sample_metadata.csv contains duplicated sample_id values.", call. = FALSE)
}
if (anyDuplicated(abundance_table$feature_id) > 0) {
  stop("abundance_table.csv contains duplicated feature_id values.", call. = FALSE)
}
if (!target_taxonomic_level %in% c("all", taxonomic_levels_for_lefse)) {
  stop("target_taxonomic_level must be 'all' or one of: ", paste(taxonomic_levels_for_lefse, collapse = ", "), call. = FALSE)
}

missing_groups <- setdiff(group_order, unique(sample_metadata[[group_column]]))
if (length(missing_groups) > 0) {
  stop("Selected group_order value(s) not found in metadata: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}

relative_abundance <- calculate_relative_abundance(abundance_table)
if (identical(target_taxonomic_level, "all")) {
  taxon_relative_abundance <- aggregate_to_all_taxonomic_levels(
    relative_abundance,
    taxonomy_table,
    taxonomic_levels_for_lefse
  )
} else {
  taxon_relative_abundance <- aggregate_to_taxon(relative_abundance, taxonomy_table, target_taxonomic_level)
}

sample_columns <- setdiff(colnames(taxon_relative_abundance), c("taxon", "rank", "taxon_name"))
taxon_relative_abundance$overall_mean_abundance <- rowMeans(taxon_relative_abundance[, sample_columns, drop = FALSE])
taxon_relative_abundance <- taxon_relative_abundance[
  order(taxon_relative_abundance$overall_mean_abundance, decreasing = TRUE),
  ,
  drop = FALSE
]
write.csv(taxon_relative_abundance, file.path(results_dir, taxon_relative_output), row.names = FALSE)

microeco_tables <- make_microeco_tables(
  sample_metadata = sample_metadata,
  abundance_table = abundance_table,
  taxonomy_table = taxonomy_table,
  group_column = group_column,
  group_order = group_order
)

microeco_dataset <- microeco::microtable$new(
  otu_table = microeco_tables$otu_table,
  sample_table = microeco_tables$sample_table,
  tax_table = microeco_tables$tax_table
)

lefse_model <- microeco::trans_diff$new(
  dataset = microeco_dataset,
  method = "lefse",
  group = "Group",
  taxa_level = target_taxonomic_level,
  alpha = lefse_alpha,
  P_adjust_method = p_adjust_method,
  lefse_subgroup = NULL,
  lefse_norm = lefse_norm
)

lefse_results <- lefse_model$res_diff
lefse_results$Taxa_clean <- strip_tax_prefixes(lefse_results$Taxa)
write.csv(lefse_results, file.path(results_dir, candidate_statistics_output), row.names = FALSE)

if (nrow(lefse_results) == 0) {
  stop("microeco LEfSe returned no biomarkers. Adjust lefse_alpha or toy data settings.", call. = FALSE)
}

plot_feature_count <- min(plot_feature_count, nrow(lefse_results))
cladogram_use_feature_num <- min(cladogram_use_feature_num, nrow(lefse_results))

lefse_for_plots <- lefse_model$clone(deep = TRUE)
lefse_for_plots$res_diff <- lefse_results[seq_len(plot_feature_count), , drop = FALSE]

lda_plot <- lefse_for_plots$plot_diff_bar(
  threshold = lda_threshold,
  color_values = unname(group_colors[group_order]),
  width = 0.5
)
lda_plot <- lda_plot +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 9, face = "italic", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

cladogram_plot <- lefse_model$plot_diff_cladogram(
  filter_taxa = cladogram_filter_taxa,
  use_taxa_num = cladogram_use_taxa_num,
  use_feature_num = cladogram_use_feature_num,
  clade_label_level = clade_label_level,
  color = unname(group_colors[group_order]),
  node_size_offset = 2,
  annotation_shape = 21,
  annotation_shape_size = 4,
  alpha = 0.2
)
cladogram_plot <- cladogram_plot +
  theme(
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA)
  )

kw_plot <- lefse_model$plot_diff_abund(
  use_number = seq_len(plot_feature_count),
  select_taxa = lefse_for_plots$plot_diff_bar_taxa,
  width = 0.5,
  color_values = unname(group_colors[group_order]),
  add_sig = TRUE,
  add_sig_label = "Significance",
  coord_flip = TRUE
)
kw_plot <- kw_plot +
  theme_classic(base_size = 11) +
  labs(x = NULL, y = "Relative abundance") +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 8.5, face = "italic", color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA)
  )

biomarker_table <- lefse_results[seq_len(plot_feature_count), , drop = FALSE]
write.csv(biomarker_table, file.path(results_dir, biomarker_table_output), row.names = FALSE)
write.csv(as.data.frame(lda_plot$data), file.path(results_dir, barplot_table_output), row.names = FALSE)
write.csv(as.data.frame(kw_plot$data), file.path(results_dir, kw_abundance_table_output), row.names = FALSE)

cladogram_data <- as.data.frame(cladogram_plot$data)
write.csv(cladogram_data, file.path(results_dir, cladogram_node_output), row.names = FALSE)
cladogram_edges <- cladogram_data[, intersect(c("parent", "node", "branch.length", "x", "y", "branch", "angle"), colnames(cladogram_data)), drop = FALSE]
write.csv(cladogram_edges, file.path(results_dir, cladogram_edge_output), row.names = FALSE)
cladogram_labels <- cladogram_data[
  !is.na(cladogram_data$node_label) & cladogram_data$node_label != "",
  intersect(c("node", "label", "node_label", "node_class", "abd", "x", "y", "angle"), colnames(cladogram_data)),
  drop = FALSE
]
write.csv(cladogram_labels, file.path(results_dir, cladogram_label_output), row.names = FALSE)

heatmap_data <- lefse_model$res_abund[lefse_model$res_abund$Taxa %in% biomarker_table$Taxa, , drop = FALSE]
heatmap_data$Taxa_clean <- strip_tax_prefixes(heatmap_data$Taxa)
heatmap_data$Mean_percent <- heatmap_data$Mean * 100
heatmap_data$Group <- factor(heatmap_data$Group, levels = group_order)
write.csv(heatmap_data, file.path(results_dir, heatmap_table_output), row.names = FALSE)

heatmap_plot <- ggplot(heatmap_data, aes(x = Group, y = Taxa_clean, fill = Mean_percent)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradientn(colors = heatmap_palette, name = "Mean relative\nabundance (%)") +
  labs(
    title = "LEfSe biomarker abundance profile",
    subtitle = "Group mean relative abundance from microeco results",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "black", face = "italic", size = 7.5),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  )

barplot_cladogram_plot <- aplot::insert_left(cladogram_plot, lda_plot, width = 0.4)
lda_kw_plot <- gridExtra::arrangeGrob(
  lda_plot + theme(legend.position = "none"),
  kw_plot,
  ncol = 2,
  widths = c(1, 1)
)

save_plot_pair(lda_plot, file.path(figures_dir, barplot_pdf), file.path(figures_dir, barplot_png), barplot_width, barplot_height, png_dpi)
save_plot_pair(cladogram_plot, file.path(figures_dir, cladogram_pdf), file.path(figures_dir, cladogram_png), cladogram_width, cladogram_height, png_dpi)
save_plot_pair(kw_plot, file.path(figures_dir, kw_abundance_pdf), file.path(figures_dir, kw_abundance_png), kw_abundance_width, kw_abundance_height, png_dpi)
save_plot_pair(barplot_cladogram_plot, file.path(figures_dir, barplot_cladogram_pdf), file.path(figures_dir, barplot_cladogram_png), barplot_cladogram_width, barplot_cladogram_height, png_dpi)
save_plot_pair(lda_kw_plot, file.path(figures_dir, lda_kw_pdf), file.path(figures_dir, lda_kw_png), lda_kw_width, lda_kw_height, png_dpi)
save_plot_pair(heatmap_plot, file.path(figures_dir, heatmap_pdf), file.path(figures_dir, heatmap_png), heatmap_width, heatmap_height, png_dpi)

cleanup_default_rplots()

message("LEfSe biomarker demo completed with native microeco plotting.")
message("Target taxonomic level: ", target_taxonomic_level)
message("microeco LEfSe biomarkers: ", nrow(lefse_results))
message("Plotted features: ", plot_feature_count)
message("Outputs written to results/ and figures/.")
