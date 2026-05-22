# Purpose:
#   Build a reproducible ternary plot showing how top microbial taxa are
#   distributed across three groups in the shared simulated metal-contaminated
#   soil microbiome toy dataset.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/taxonomy_table.csv
#
# Output files:
#   - results/relative_abundance_by_feature.csv
#   - results/group_mean_abundance_by_<target_taxonomic_level>.csv
#   - results/ternary_plotting_table_<target_taxonomic_level>.csv
#   - figures/<output_pdf_file>
#   - figures/<output_png_file>
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   abundance_table.csv: feature_id, one column per sample_id
#   taxonomy_table.csv: feature_id and the selected target_taxonomic_level
#
# User-editable settings:
#   Edit the settings block below to change input paths, axis groups,
#   taxonomy level, top_n, output names, figure size, and color palette.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")

target_taxonomic_level <- "Genus"

group_for_axis_1 <- "Control"
group_for_axis_2 <- "Mining"
group_for_axis_3 <- "Smelting"

top_n <- 12

relative_abundance_output <- "relative_abundance_by_feature.csv"
group_mean_output_prefix <- "group_mean_abundance_by"
ternary_table_output_prefix <- "ternary_plotting_table"

output_pdf_file <- "ternary_taxa_distribution.pdf"
output_png_file <- "ternary_taxa_distribution.png"

figure_width <- 7.2
figure_height <- 6.2
png_dpi <- 300

color_palette <- c(
  "#ff6e61", "#ffd561", "#6ccb79", "#4d97ff",
  "#ff70c8", "#ff9f1a", "#2ec2b3", "#9e4fde",
  "#ff575c", "#1981c2", "#91abad", "#816d5a"
)

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

if (!requireNamespace("ggtern", quietly = TRUE)) {
  stop(
    "The R package 'ggtern' is required for this module but is not installed.\n",
    "Please install it before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(ggtern))

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

input_files <- c(metadata_file, abundance_file, taxonomy_file)
check_files_exist(input_files)

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_table <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
taxonomy_table <- read.csv(taxonomy_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", "group"), "sample_metadata.csv")
check_required_columns(abundance_table, "feature_id", "abundance_table.csv")
check_required_columns(taxonomy_table, c("feature_id", target_taxonomic_level), "taxonomy_table.csv")

axis_groups <- c(group_for_axis_1, group_for_axis_2, group_for_axis_3)
missing_axis_groups <- setdiff(axis_groups, unique(sample_metadata$group))
if (length(missing_axis_groups) > 0) {
  stop(
    "Selected axis group(s) are not present in sample_metadata.csv: ",
    paste(missing_axis_groups, collapse = ", "),
    call. = FALSE
  )
}

check_sample_alignment(sample_metadata, abundance_table)
check_feature_alignment(abundance_table, taxonomy_table)

sample_columns <- setdiff(colnames(abundance_table), "feature_id")
sample_metadata <- sample_metadata[match(sample_columns, sample_metadata$sample_id), , drop = FALSE]

relative_abundance <- calculate_relative_abundance(abundance_table)
relative_abundance_path <- file.path(results_dir, relative_abundance_output)
write.csv(relative_abundance, relative_abundance_path, row.names = FALSE)

taxon_abundance <- aggregate_by_taxonomy(
  relative_abundance = relative_abundance,
  taxonomy = taxonomy_table,
  target_taxonomic_level = target_taxonomic_level
)

group_mean_abundance <- calculate_group_means(
  taxon_abundance = taxon_abundance,
  metadata = sample_metadata
)

safe_level_name <- gsub("[^A-Za-z0-9_]+", "_", target_taxonomic_level)
group_mean_output <- paste0(group_mean_output_prefix, "_", safe_level_name, ".csv")
group_mean_path <- file.path(results_dir, group_mean_output)
write.csv(group_mean_abundance, group_mean_path, row.names = FALSE)

ternary_table <- make_ternary_table(
  group_means = group_mean_abundance,
  axis_groups = axis_groups,
  top_n = top_n
)

ternary_table_output <- paste0(ternary_table_output_prefix, "_", safe_level_name, ".csv")
ternary_table_path <- file.path(results_dir, ternary_table_output)
write.csv(ternary_table, ternary_table_path, row.names = FALSE)

plot_palette <- rep(color_palette, length.out = nrow(ternary_table))
names(plot_palette) <- ternary_table$taxon

ternary_plot <- ggtern(
  data = ternary_table,
  aes(
    x = axis_1_proportion,
    y = axis_2_proportion,
    z = axis_3_proportion,
    color = taxon,
    size = overall_mean_percent
  )
) +
  geom_point(alpha = 0.88) +
  scale_color_manual(values = plot_palette) +
  scale_size_continuous(range = c(2.8, 8.5), name = "Mean abundance (%)") +
  labs(
    title = paste0("Top ", nrow(ternary_table), " taxa distribution"),
    subtitle = paste0("Taxonomic level: ", target_taxonomic_level),
    x = group_for_axis_1,
    y = group_for_axis_2,
    z = group_for_axis_3,
    color = target_taxonomic_level
  ) +
  theme_rgbw(base_size = 12) +
  theme_showarrows() +
  theme_clockwise() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.35),
    panel.grid.minor = element_blank()
  )

pdf_path <- file.path(figures_dir, output_pdf_file)
png_path <- file.path(figures_dir, output_png_file)

ggsave(pdf_path, ternary_plot, width = figure_width, height = figure_height, units = "in")
ggsave(png_path, ternary_plot, width = figure_width, height = figure_height, units = "in", dpi = png_dpi)

message("Ternary taxa distribution demo completed.")
message("Wrote: ", relative_abundance_path)
message("Wrote: ", group_mean_path)
message("Wrote: ", ternary_table_path)
message("Wrote: ", pdf_path)
message("Wrote: ", png_path)
