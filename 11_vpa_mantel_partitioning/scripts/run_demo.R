# Purpose:
#   Build an independent VPA and Mantel partitioning demo from shared simulated
#   soil microbiome toy data. The workflow follows the original scripts:
#   Hellinger-transformed species/function matrices, vegan::varpart() for
#   Sb / co-metal / nutrient variable groups, partial RDA testing, and a
#   ggcor Mantel-link plot over an environmental Spearman correlation heatmap.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/environmental_variables.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/functional_annotation_table.csv
#
# Output files:
#   - results/vpa_fraction_table.csv
#   - results/partial_rda_tests.csv
#   - results/mantel_results.csv
#   - results/environment_spearman_correlation.csv
#   - results/species_response_matrix.csv
#   - results/function_response_matrix.csv
#   - figures/fig1_vpa_combined.pdf/png
#   - figures/fig2_mantel_env_correlation.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   environmental_variables.csv: sample_id plus selected environmental variables
#   abundance_table.csv: feature_id plus one numeric column per sample_id
#   functional_annotation_table.csv: sample_id, selected functional feature column,
#                                    selected numeric value column
#
# User-editable settings:
#   Edit the settings block below to change input paths, environmental variable
#   groups, feature selection, Mantel response blocks, plot colors, output file
#   names, figure sizes, and random seed.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
environmental_file <- file.path(shared_data_dir, "environmental_variables.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
functional_annotation_file <- file.path(shared_data_dir, "functional_annotation_table.csv")

group_column <- "group"
random_seed <- 123
permutations <- 999

top_n_taxa_for_vpa <- 60
top_n_functions_for_vpa <- 40
functional_feature_column <- "ko_id"
functional_value_column <- "abundance"

sb_variables <- c("Sb_total")
co_metal_variables <- c("Cu", "As")
nutrient_variables <- c("TOC", "SO4", "pH")
mantel_environmental_variables <- c("pH", "EC", "TOC", "TN", "SO4", "NO3", "Sb_total", "Sb_III", "Sb_V", "As", "Cu")

vpa_group_names <- c("Sb", "Cu & As", "Nutrients")
vpa_colors <- c("#E64B3599", "#4DBBD599", "#00A08799")

mantel_r_breaks <- c(-Inf, 0.2, 0.4, Inf)
mantel_r_labels <- c("<0.2", "0.2-0.4", ">=0.4")
mantel_p_breaks <- c(-Inf, 0.01, 0.05, Inf)
mantel_p_labels <- c("<0.01", "0.01-0.05", ">=0.05")
mantel_link_colors <- c("<0.01" = "#62a11b", "0.01-0.05" = "#68edcb", ">=0.05" = "snow3")
mantel_link_sizes <- c("<0.2" = 0.35, "0.2-0.4" = 0.7, ">=0.4" = 1.25)
mantel_heatmap_colors <- c(low = "#c6f093", mid = "#f8fff8", high = "#163f00")

vpa_fraction_output <- "vpa_fraction_table.csv"
partial_rda_output <- "partial_rda_tests.csv"
mantel_output <- "mantel_results.csv"
env_correlation_output <- "environment_spearman_correlation.csv"
species_matrix_output <- "species_response_matrix.csv"
function_matrix_output <- "function_response_matrix.csv"

vpa_pdf <- "fig1_vpa_combined.pdf"
vpa_png <- "fig1_vpa_combined.png"
mantel_pdf <- "fig2_mantel_env_correlation.pdf"
mantel_png <- "fig2_mantel_env_correlation.png"

vpa_width <- 10
vpa_height <- 5
mantel_width <- 8
mantel_height <- 6
png_dpi <- 300

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("vegan", "ggplot2", "ggcor", "dplyr")
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
  library(vegan)
  library(ggplot2)
  library(ggcor)
  library(dplyr)
})

set.seed(random_seed)

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(metadata_file, environmental_file, abundance_file, functional_annotation_file))

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
environmental_variables <- read.csv(environmental_file, stringsAsFactors = FALSE, check.names = FALSE)
functional_annotation <- read.csv(functional_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_matrix <- read_abundance_matrix(abundance_file)

check_required_columns(metadata, c("sample_id", group_column), "sample_metadata.csv")
metadata <- metadata[metadata$sample_id %in% colnames(abundance_matrix), , drop = FALSE]
metadata <- metadata[match(colnames(abundance_matrix), metadata$sample_id), , drop = FALSE]
if (anyNA(metadata$sample_id)) stop("sample_metadata.csv and abundance_table.csv sample IDs do not fully match.", call. = FALSE)

environmental_vars <- unique(c(sb_variables, co_metal_variables, nutrient_variables, mantel_environmental_variables))
env_raw <- align_environment(environmental_variables, metadata, environmental_vars)

species_matrix <- make_taxonomy_response_matrix(abundance_matrix, metadata, top_n_taxa_for_vpa)
function_matrix <- make_function_response_matrix(
  functional_annotation = functional_annotation,
  metadata = metadata,
  value_column = functional_value_column,
  feature_column = functional_feature_column,
  top_n_features = top_n_functions_for_vpa
)

write.csv(data.frame(sample_id = rownames(species_matrix), species_matrix, check.names = FALSE), file.path(results_dir, species_matrix_output), row.names = FALSE)
write.csv(data.frame(sample_id = rownames(function_matrix), function_matrix, check.names = FALSE), file.path(results_dir, function_matrix_output), row.names = FALSE)

species_hellinger <- vegan::decostand(species_matrix, method = "hellinger")
function_hellinger <- vegan::decostand(function_matrix, method = "hellinger")
env_scaled <- as.data.frame(scale(env_raw))

env_sb <- env_scaled[, sb_variables, drop = FALSE]
env_co_metals <- env_scaled[, co_metal_variables, drop = FALSE]
env_nutrients <- env_scaled[, nutrient_variables, drop = FALSE]

vpa_species <- vegan::varpart(species_hellinger, env_sb, env_co_metals, env_nutrients)
vpa_genes <- vegan::varpart(function_hellinger, env_sb, env_co_metals, env_nutrients)

vpa_fraction_table <- rbind(
  clean_vpa_fractions(vpa_species, "Species community"),
  clean_vpa_fractions(vpa_genes, "Functional genes")
)
write.csv(vpa_fraction_table, file.path(results_dir, vpa_fraction_output), row.names = FALSE)

partial_species <- run_partial_rda_test(
  response_matrix = species_matrix,
  focal_env = env_sb,
  covariates = cbind(env_co_metals, env_nutrients),
  permutations = permutations,
  seed = random_seed
)
partial_species$response <- "Species community"

partial_genes <- run_partial_rda_test(
  response_matrix = function_matrix,
  focal_env = env_sb,
  covariates = cbind(env_co_metals, env_nutrients),
  permutations = permutations,
  seed = random_seed
)
partial_genes$response <- "Functional genes"

partial_tests <- rbind(partial_species, partial_genes)
partial_tests$term <- rownames(partial_tests)
partial_tests <- partial_tests[, c("response", "term", setdiff(colnames(partial_tests), c("response", "term"))), drop = FALSE]
write.csv(partial_tests, file.path(results_dir, partial_rda_output), row.names = FALSE)

save_base_vpa_pair(
  vpa_species = vpa_species,
  vpa_genes = vpa_genes,
  pdf_file = file.path(figures_dir, vpa_pdf),
  png_file = file.path(figures_dir, vpa_png),
  width = vpa_width,
  height = vpa_height,
  colors = vpa_colors,
  group_names = vpa_group_names
)

mantel_env <- env_raw[, mantel_environmental_variables, drop = FALSE]
mantel_response <- cbind(
  species_matrix,
  function_matrix
)
species_indices <- seq_len(ncol(species_matrix))
function_indices <- seq(from = ncol(species_matrix) + 1, length.out = ncol(function_matrix))

mantel_results <- ggcor::mantel_test(
  mantel_response,
  mantel_env,
  mantel.fun = "mantel",
  spec.dist.method = "bray",
  env.dist.method = "euclidean",
  spec.select = list(
    Species_community = species_indices,
    Functional_genes = function_indices
  ),
  env.select = NULL
) %>%
  mutate(
    rd = cut(r, breaks = mantel_r_breaks, labels = mantel_r_labels),
    pd = cut(p.value, breaks = mantel_p_breaks, labels = mantel_p_labels)
  )

write.csv(as.data.frame(mantel_results), file.path(results_dir, mantel_output), row.names = FALSE)

env_correlation <- stats::cor(mantel_env, method = "spearman")
write.csv(env_correlation, file.path(results_dir, env_correlation_output), row.names = TRUE)

corr_df <- as.data.frame(as.table(env_correlation), stringsAsFactors = FALSE)
colnames(corr_df) <- c("var_y", "var_x", "spearman_r")
corr_df$x <- match(corr_df$var_x, mantel_environmental_variables)
corr_df$y <- match(corr_df$var_y, mantel_environmental_variables)
corr_df <- corr_df[corr_df$x > corr_df$y, , drop = FALSE]

axis_df <- data.frame(
  variable = mantel_environmental_variables,
  x = seq_along(mantel_environmental_variables),
  y = seq_along(mantel_environmental_variables),
  stringsAsFactors = FALSE
)

mantel_plot_data <- as.data.frame(mantel_results)
mantel_plot_data$pd <- factor(mantel_plot_data$pd, levels = mantel_p_labels)
mantel_plot_data$rd <- factor(mantel_plot_data$rd, levels = mantel_r_labels)
mantel_plot_data$x <- -0.45
mantel_plot_data$y <- ifelse(mantel_plot_data$spec == "Species_community", length(mantel_environmental_variables) + 1.65, length(mantel_environmental_variables) + 0.95)
mantel_plot_data$xend <- match(mantel_plot_data$env, mantel_environmental_variables)
mantel_plot_data$yend <- match(mantel_plot_data$env, mantel_environmental_variables)
mantel_plot_data$curvature <- ifelse(mantel_plot_data$spec == "Species_community", -0.25, -0.15)

response_label_df <- unique(mantel_plot_data[, c("spec", "x", "y"), drop = FALSE])
response_label_df$label <- gsub("_", " ", response_label_df$spec)
response_label_df$label_x <- response_label_df$x - 1.45

mantel_plot <- ggplot() +
  geom_tile(data = corr_df, aes(x = x, y = y, fill = spearman_r), color = "white", linewidth = 0.4) +
  geom_text(data = corr_df, aes(x = x, y = y, label = sprintf("%.2f", spearman_r)), size = 2.6, color = "black") +
  geom_text(data = axis_df, aes(x = x, y = -0.15, label = variable), angle = 45, hjust = 1, vjust = 1, size = 3.2, fontface = "bold") +
  geom_text(data = axis_df, aes(x = 0.05, y = y, label = variable), hjust = 1, size = 3.2, fontface = "bold") +
  geom_curve(
    data = mantel_plot_data,
    aes(x = x, y = y, xend = xend, yend = yend, colour = pd, linewidth = rd),
    curvature = -0.22,
    alpha = 0.62,
    lineend = "round"
  ) +
  geom_point(data = response_label_df, aes(x = x, y = y), shape = 21, fill = "#f8fff8", color = "#163f00", size = 3.6, stroke = 1) +
  geom_label(
    data = response_label_df,
    aes(x = label_x, y = y, label = label),
    hjust = 0,
    size = 3.4,
    fontface = "bold",
    fill = "white",
    color = "black",
    linewidth = 0,
    label.padding = unit(0.12, "lines")
  ) +
  scale_fill_gradient2(
    low = mantel_heatmap_colors["low"],
    mid = mantel_heatmap_colors["mid"],
    high = mantel_heatmap_colors["high"],
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_color_manual(values = mantel_link_colors, drop = FALSE) +
  scale_linewidth_manual(values = mantel_link_sizes, drop = FALSE) +
  coord_equal(
    xlim = c(-2.25, length(mantel_environmental_variables) + 0.8),
    ylim = c(-0.95, length(mantel_environmental_variables) + 2.15),
    clip = "off"
  ) +
  guides(
    linewidth = guide_legend(title = "Mantel's r", order = 2),
    colour = guide_legend(title = "Mantel's p", order = 3),
    fill = guide_colorbar(title = "Spearman's r", order = 4)
  ) +
  labs(x = NULL, y = NULL) +
  theme_void(base_size = 11) +
  theme(
    text = element_text(size = 11, face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    plot.margin = margin(18, 18, 34, 54)
  )

save_ggplot_pair(
  mantel_plot,
  file.path(figures_dir, mantel_pdf),
  file.path(figures_dir, mantel_png),
  mantel_width,
  mantel_height,
  png_dpi
)

message("VPA/Mantel partitioning demo completed.")
message("Species matrix: ", nrow(species_matrix), " samples x ", ncol(species_matrix), " features")
message("Function matrix: ", nrow(function_matrix), " samples x ", ncol(function_matrix), " features")
message("VPA groups: ", paste(vpa_group_names, collapse = " / "))
message("Mantel environmental variables: ", paste(mantel_environmental_variables, collapse = ", "))
message("Outputs written to results/ and figures/.")
