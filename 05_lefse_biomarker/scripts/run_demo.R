# Purpose:
#   Build a reproducible LEfSe-style biomarker visualization from shared
#   simulated soil microbiome taxonomy and abundance toy data.
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
#   - results/biomarker_group_heatmap_table.csv
#   - figures/lefse_biomarker_barplot.pdf
#   - figures/lefse_biomarker_barplot.png
#   - figures/biomarker_group_heatmap.pdf
#   - figures/biomarker_group_heatmap.png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   abundance_table.csv: feature_id plus one column per sample_id
#   taxonomy_table.csv: feature_id and the selected target_taxonomic_level
#
# User-editable settings:
#   Edit the settings block below to change input paths, taxonomy level,
#   group order, biomarker filters, output names, figure sizes, and palettes.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")

target_taxonomic_level <- "Genus"
group_column <- "group"
group_order <- c("Control", "Tailing", "Mining", "Smelting")

top_n_taxa <- 30
max_biomarkers_to_plot <- 18
minimum_biomarkers_to_plot <- 8
min_prevalence <- 0.15
min_mean_relative_abundance <- 0.001
kruskal_fdr_cutoff <- 0.20
effect_score_cutoff <- 0.20
pseudocount <- 1e-6

taxon_relative_output <- "taxon_relative_abundance_by_sample.csv"
candidate_statistics_output <- "lefse_candidate_statistics.csv"
biomarker_table_output <- "lefse_biomarker_table.csv"
barplot_table_output <- "lefse_barplot_plotting_table.csv"
heatmap_table_output <- "biomarker_group_heatmap_table.csv"

barplot_pdf <- "lefse_biomarker_barplot.pdf"
barplot_png <- "lefse_biomarker_barplot.png"
heatmap_pdf <- "biomarker_group_heatmap.pdf"
heatmap_png <- "biomarker_group_heatmap.png"

barplot_width <- 7.8
barplot_height <- 6.2
heatmap_width <- 6.8
heatmap_height <- 6.2
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

required_packages <- c("ggplot2")
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
check_required_columns(taxonomy_table, c("feature_id", target_taxonomic_level), "taxonomy_table.csv")
check_sample_ids(sample_metadata, abundance_table)
check_feature_ids(abundance_table, taxonomy_table)

if (anyDuplicated(sample_metadata$sample_id) > 0) {
  stop("sample_metadata.csv contains duplicated sample_id values.", call. = FALSE)
}
if (anyDuplicated(abundance_table$feature_id) > 0) {
  stop("abundance_table.csv contains duplicated feature_id values.", call. = FALSE)
}
if (top_n_taxa < 2) stop("top_n_taxa must be at least 2.", call. = FALSE)
if (max_biomarkers_to_plot < 2) stop("max_biomarkers_to_plot must be at least 2.", call. = FALSE)
if (minimum_biomarkers_to_plot < 1) stop("minimum_biomarkers_to_plot must be at least 1.", call. = FALSE)

missing_groups <- setdiff(group_order, unique(sample_metadata[[group_column]]))
if (length(missing_groups) > 0) {
  stop("Selected group_order value(s) not found in metadata: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}

sample_metadata$group <- sample_metadata[[group_column]]
sample_metadata$group <- factor(sample_metadata$group, levels = group_order)

relative_abundance <- calculate_relative_abundance(abundance_table)
taxon_relative_abundance <- aggregate_to_taxon(relative_abundance, taxonomy_table, target_taxonomic_level)

sample_columns <- setdiff(colnames(taxon_relative_abundance), "taxon")
taxon_relative_abundance$overall_mean_abundance <- rowMeans(taxon_relative_abundance[, sample_columns, drop = FALSE])
taxon_relative_abundance <- taxon_relative_abundance[
  order(taxon_relative_abundance$overall_mean_abundance, decreasing = TRUE),
  ,
  drop = FALSE
]

if (top_n_taxa > nrow(taxon_relative_abundance)) {
  message(
    "top_n_taxa (", top_n_taxa,
    ") is larger than the number of available taxa (", nrow(taxon_relative_abundance),
    "); using all available taxa."
  )
  top_n_taxa <- nrow(taxon_relative_abundance)
}

taxon_relative_abundance <- head(taxon_relative_abundance, top_n_taxa)
write.csv(taxon_relative_abundance, file.path(results_dir, taxon_relative_output), row.names = FALSE)

long_taxa <- taxon_table_to_long(taxon_relative_abundance[, c("taxon", sample_columns), drop = FALSE])
long_taxa <- merge(
  long_taxa,
  sample_metadata[, c("sample_id", "group")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

candidate_statistics <- make_biomarker_statistics(long_taxa, group_order, pseudocount)
candidate_statistics$overall_mean_percent <- candidate_statistics$overall_mean_abundance * 100
candidate_statistics$mean_enriched_group_percent <- candidate_statistics$mean_enriched_group * 100
candidate_statistics$mean_other_groups_percent <- candidate_statistics$mean_other_groups * 100
write.csv(candidate_statistics, file.path(results_dir, candidate_statistics_output), row.names = FALSE)

biomarkers <- select_biomarkers(
  candidate_statistics,
  min_prevalence = min_prevalence,
  min_mean_abundance = min_mean_relative_abundance,
  fdr_cutoff = kruskal_fdr_cutoff,
  effect_score_cutoff = effect_score_cutoff,
  max_biomarkers_to_plot = max_biomarkers_to_plot,
  minimum_biomarkers_to_plot = minimum_biomarkers_to_plot
)

if (nrow(biomarkers) == 0) {
  stop("No biomarkers available after filtering. Lower the filtering thresholds and rerun.", call. = FALSE)
}

biomarkers$enriched_group <- factor(biomarkers$enriched_group, levels = group_order)
biomarkers <- biomarkers[order(biomarkers$enriched_group, biomarkers$lefse_like_score), , drop = FALSE]
biomarkers$taxon <- factor(biomarkers$taxon, levels = biomarkers$taxon)
write.csv(biomarkers, file.path(results_dir, biomarker_table_output), row.names = FALSE)
write.csv(biomarkers, file.path(results_dir, barplot_table_output), row.names = FALSE)

selected_taxa <- as.character(biomarkers$taxon)
heatmap_data <- long_taxa[long_taxa$taxon %in% selected_taxa, , drop = FALSE]
heatmap_summary <- aggregate(
  relative_abundance ~ taxon + group,
  data = heatmap_data,
  FUN = mean
)
heatmap_summary$mean_percent <- heatmap_summary$relative_abundance * 100
heatmap_summary$group <- factor(heatmap_summary$group, levels = group_order)
heatmap_summary$taxon <- factor(heatmap_summary$taxon, levels = levels(biomarkers$taxon))
write.csv(heatmap_summary, file.path(results_dir, heatmap_table_output), row.names = FALSE)

bar_plot <- ggplot(biomarkers, aes(x = lefse_like_score, y = taxon, fill = enriched_group)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.18) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "LEfSe-style biomarkers",
    subtitle = paste0(target_taxonomic_level, " level; score = log10 enriched-group ratio"),
    x = "LEfSe-style effect score",
    y = NULL,
    fill = "Enriched group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black", face = "italic"),
    axis.text.x = element_text(color = "black"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

heatmap_plot <- ggplot(heatmap_summary, aes(x = group, y = taxon, fill = mean_percent)) +
  geom_tile(color = "white", linewidth = 0.45) +
  scale_fill_gradientn(colors = heatmap_palette, name = "Mean relative\nabundance (%)") +
  labs(
    title = "Biomarker abundance profile",
    subtitle = "Group mean relative abundance",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", angle = 35, hjust = 1),
    axis.text.y = element_text(color = "black", face = "italic"),
    legend.position = "right",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(figures_dir, barplot_pdf), bar_plot, width = barplot_width, height = barplot_height, units = "in")
ggsave(file.path(figures_dir, barplot_png), bar_plot, width = barplot_width, height = barplot_height, units = "in", dpi = png_dpi)
ggsave(file.path(figures_dir, heatmap_pdf), heatmap_plot, width = heatmap_width, height = heatmap_height, units = "in")
ggsave(file.path(figures_dir, heatmap_png), heatmap_plot, width = heatmap_width, height = heatmap_height, units = "in", dpi = png_dpi)

message("LEfSe-style biomarker demo completed.")
message("Target taxonomic level: ", target_taxonomic_level)
message("Tested taxa: ", nrow(candidate_statistics))
message("Plotted biomarkers: ", nrow(biomarkers))
message("Outputs written to results/ and figures/.")
