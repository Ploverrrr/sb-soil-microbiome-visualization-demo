# Purpose:
#   Build a reproducible FAPROTAX-style functional profile visualization from
#   the shared simulated soil microbiome functional annotation toy data.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/functional_relative_abundance_by_sample.csv
#   - results/functional_group_summary.csv
#   - results/functional_bubble_plotting_table.csv
#   - results/functional_barplot_table.csv
#   - figures/faprotax_function_bubble_profile.pdf
#   - figures/faprotax_function_bubble_profile.png
#   - figures/faprotax_group_barplot.pdf
#   - figures/faprotax_group_barplot.png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   functional_annotation_table.csv: sample_id, annotation_source,
#     function_category, pathway, count
#
# User-editable settings:
#   Edit the settings block below to change input paths, selected annotation
#   sources, function grouping column, top_n, group order/colors, output names,
#   and figure sizes.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

annotation_sources_to_include <- c("FAPROTAX", "SulfurCycle", "NitrogenCycle", "CarbonCycle")
function_group_column <- "pathway"
top_n_functions <- 10

group_order <- c("Control", "Tailing", "Mining", "Smelting")
group_colors <- c(
  Control = "#b39cd0",
  Tailing = "#ffc75f",
  Mining = "#ff6f91",
  Smelting = "#4d97ff"
)

relative_abundance_output <- "functional_relative_abundance_by_sample.csv"
group_summary_output <- "functional_group_summary.csv"
bubble_plotting_output <- "functional_bubble_plotting_table.csv"
barplot_table_output <- "functional_barplot_table.csv"

bubble_pdf <- "faprotax_function_bubble_profile.pdf"
bubble_png <- "faprotax_function_bubble_profile.png"
barplot_pdf <- "faprotax_group_barplot.pdf"
barplot_png <- "faprotax_group_barplot.png"

bubble_width <- 9
bubble_height <- 5.8
barplot_width <- 8.5
barplot_height <- 5.6
png_dpi <- 300

bubble_color_palette <- c("#00c9c8", "white", "#ed5e93")
bubble_size_range <- c(1.8, 8.5)

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

check_files_exist(c(metadata_file, functional_annotation_file))

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", "group"), "sample_metadata.csv")
check_required_columns(
  functional_annotation,
  c("sample_id", "annotation_source", "function_category", function_group_column, "count"),
  "functional_annotation_table.csv"
)
check_sample_ids(sample_metadata, functional_annotation)

missing_sources <- setdiff(annotation_sources_to_include, unique(functional_annotation$annotation_source))
if (length(missing_sources) > 0) {
  stop("Selected annotation source(s) not found: ", paste(missing_sources, collapse = ", "), call. = FALSE)
}

missing_groups <- setdiff(group_order, unique(sample_metadata$group))
if (length(missing_groups) > 0) {
  stop("Selected group_order value(s) not found in metadata: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}

filtered_functional <- functional_annotation[
  functional_annotation$annotation_source %in% annotation_sources_to_include,
  ,
  drop = FALSE
]

aggregated_counts <- aggregate_function_counts(filtered_functional, function_group_column)
relative_abundance <- make_relative_function_matrix(aggregated_counts)

if (top_n_functions < 2) stop("top_n_functions must be at least 2.", call. = FALSE)
if (top_n_functions > nrow(relative_abundance)) {
  message(
    "top_n_functions (", top_n_functions,
    ") is larger than the number of available functions (", nrow(relative_abundance),
    "); using all available functions."
  )
  top_n_functions <- nrow(relative_abundance)
}

top_functions <- select_top_functions(relative_abundance, top_n_functions)
write.csv(top_functions, file.path(results_dir, relative_abundance_output), row.names = FALSE)

long_top_functions <- long_function_table(top_functions)
long_top_functions <- merge(
  long_top_functions,
  sample_metadata[, c("sample_id", "group")],
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)
long_top_functions$group <- factor(long_top_functions$group, levels = group_order)
sample_order <- sample_metadata$sample_id[order(match(sample_metadata$group, group_order), sample_metadata$sample_id)]
if ("replicate" %in% colnames(sample_metadata)) {
  sample_order <- sample_metadata$sample_id[order(match(sample_metadata$group, group_order), sample_metadata$replicate)]
}
long_top_functions$sample_id <- factor(long_top_functions$sample_id, levels = sample_order)
long_top_functions$function_name <- factor(long_top_functions$function_name, levels = rev(top_functions$function_name))
long_top_functions$relative_abundance_percent <- long_top_functions$relative_abundance * 100

write.csv(long_top_functions, file.path(results_dir, bubble_plotting_output), row.names = FALSE)

group_summary <- summarize_by_group(long_top_functions, sample_metadata)
group_summary$group <- factor(group_summary$group, levels = group_order)
group_summary$function_name <- factor(group_summary$function_name, levels = rev(top_functions$function_name))
group_summary$mean_percent <- group_summary$mean_relative_abundance * 100
group_summary$se_percent <- group_summary$se_relative_abundance * 100
write.csv(group_summary, file.path(results_dir, group_summary_output), row.names = FALSE)
write.csv(group_summary, file.path(results_dir, barplot_table_output), row.names = FALSE)

bubble_plot <- ggplot(long_top_functions, aes(x = sample_id, y = function_name)) +
  geom_point(aes(size = relative_abundance_percent, color = relative_abundance_percent), alpha = 0.88) +
  scale_color_gradientn(colors = bubble_color_palette, name = "Relative abundance (%)") +
  scale_size_continuous(range = bubble_size_range, name = "Relative abundance (%)") +
  labs(
    title = "FAPROTAX-style functional profile",
    subtitle = paste0("Selected sources: ", paste(annotation_sources_to_include, collapse = ", ")),
    x = NULL,
    y = function_group_column
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  )

bar_plot <- ggplot(group_summary, aes(x = mean_percent, y = function_name, fill = group)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  geom_errorbar(
    aes(xmin = pmax(mean_percent - se_percent, 0), xmax = mean_percent + se_percent),
    position = position_dodge(width = 0.78),
    width = 0.28,
    linewidth = 0.35
  ) +
  scale_fill_manual(values = group_colors, drop = FALSE) +
  labs(
    title = "Group-level functional profile",
    subtitle = "Mean relative abundance +/- SE",
    x = "Relative abundance (%)",
    y = function_group_column,
    fill = "Group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )

ggsave(file.path(figures_dir, bubble_pdf), bubble_plot, width = bubble_width, height = bubble_height, units = "in")
ggsave(file.path(figures_dir, bubble_png), bubble_plot, width = bubble_width, height = bubble_height, units = "in", dpi = png_dpi)
ggsave(file.path(figures_dir, barplot_pdf), bar_plot, width = barplot_width, height = barplot_height, units = "in")
ggsave(file.path(figures_dir, barplot_png), bar_plot, width = barplot_width, height = barplot_height, units = "in", dpi = png_dpi)

message("FAPROTAX-style functional profile demo completed.")
message("Selected annotation sources: ", paste(annotation_sources_to_include, collapse = ", "))
message("Top functions: ", nrow(top_functions))
message("Outputs written to results/ and figures/.")
