# Purpose:
#   Build a reproducible microbe-environment association network from shared
#   simulated metal-contaminated soil microbiome toy data.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/environmental_variables.csv
#   - ../data/toy_shared/abundance_table.csv
#   - ../data/toy_shared/taxonomy_table.csv
#
# Output files:
#   - results/relative_abundance_by_feature.csv
#   - results/taxon_abundance_by_sample.csv
#   - results/microbe_env_correlation_results.csv
#   - results/microbe_env_corr_matrix.csv
#   - results/network_edge_list.csv
#   - results/network_node_list.csv
#   - results/network_summary.csv
#   - figures/microbe_env_network.pdf
#   - figures/microbe_env_network.png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   environmental_variables.csv: sample_id and selected environmental variables
#   abundance_table.csv: feature_id, one column per sample_id
#   taxonomy_table.csv: feature_id and selected target_taxonomic_level
#
# User-editable settings:
#   Edit the settings block below to change input paths, taxonomy level,
#   selected environmental variables, thresholds, output names, figure size,
#   and network styling.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
environment_file <- file.path(shared_data_dir, "environmental_variables.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")
taxonomy_file <- file.path(shared_data_dir, "taxonomy_table.csv")

target_taxonomic_level <- "Genus"
top_n_taxa <- 14

environmental_variables <- c("pH", "EC", "TOC", "SO4", "NO3", "Sb_total", "As", "Cu", "Zn", "Cd")

correlation_method <- "spearman"
p_adjust_method <- "BH"
correlation_threshold <- 0.55
adjusted_p_threshold <- 0.05

relative_abundance_output <- "relative_abundance_by_feature.csv"
taxon_abundance_output <- "taxon_abundance_by_sample.csv"
correlation_output <- "microbe_env_correlation_results.csv"
correlation_matrix_output <- "microbe_env_corr_matrix.csv"
edge_list_output <- "network_edge_list.csv"
node_list_output <- "network_node_list.csv"
network_summary_output <- "network_summary.csv"

network_pdf <- "microbe_env_network.pdf"
network_png <- "microbe_env_network.png"

figure_width <- 8.5
figure_height <- 6.6
png_dpi <- 300

taxon_node_color <- "#4d97ff"
environment_node_color <- "#ff9f1a"
positive_edge_color <- "#ed5e93"
negative_edge_color <- "#00c9c8"
node_size_range <- c(4, 10)
edge_width_range <- c(0.4, 2.6)
layout_seed <- 20260522

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("igraph")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(igraph))

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

input_files <- c(metadata_file, environment_file, abundance_file, taxonomy_file)
check_files_exist(input_files)

sample_metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
environmental_data <- read.csv(environment_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_table <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
taxonomy_table <- read.csv(taxonomy_file, stringsAsFactors = FALSE, check.names = FALSE)

check_required_columns(sample_metadata, c("sample_id", "group"), "sample_metadata.csv")
check_required_columns(environmental_data, c("sample_id", environmental_variables), "environmental_variables.csv")
check_required_columns(abundance_table, "feature_id", "abundance_table.csv")
check_required_columns(taxonomy_table, c("feature_id", target_taxonomic_level), "taxonomy_table.csv")

non_numeric_env <- environmental_variables[!vapply(environmental_data[, environmental_variables, drop = FALSE], is.numeric, logical(1))]
if (length(non_numeric_env) > 0) {
  stop("Selected environmental variable(s) must be numeric: ", paste(non_numeric_env, collapse = ", "), call. = FALSE)
}

check_sample_alignment(sample_metadata, environmental_data, abundance_table)
check_feature_alignment(abundance_table, taxonomy_table)

if (top_n_taxa < 2) stop("top_n_taxa must be at least 2.", call. = FALSE)
if (correlation_threshold < 0 || correlation_threshold > 1) stop("correlation_threshold must be between 0 and 1.", call. = FALSE)
if (adjusted_p_threshold <= 0 || adjusted_p_threshold > 1) stop("adjusted_p_threshold must be in (0, 1].", call. = FALSE)

sample_columns <- setdiff(colnames(abundance_table), "feature_id")
sample_metadata <- sample_metadata[match(sample_columns, sample_metadata$sample_id), , drop = FALSE]
environmental_data <- environmental_data[match(sample_columns, environmental_data$sample_id), , drop = FALSE]

relative_abundance <- calculate_relative_abundance(abundance_table)
write.csv(relative_abundance, file.path(results_dir, relative_abundance_output), row.names = FALSE)

taxon_abundance <- aggregate_by_taxonomy(relative_abundance, taxonomy_table, target_taxonomic_level)
if (top_n_taxa > nrow(taxon_abundance)) {
  stop("top_n_taxa is larger than the number of available taxa: ", nrow(taxon_abundance), call. = FALSE)
}
selected_taxa <- select_top_taxa(taxon_abundance, top_n_taxa)
write.csv(selected_taxa, file.path(results_dir, taxon_abundance_output), row.names = FALSE)

correlation_results <- calculate_microbe_env_correlations(
  taxon_abundance = selected_taxa,
  environmental_data = environmental_data,
  environmental_variables = environmental_variables,
  correlation_method = correlation_method,
  p_adjust_method = p_adjust_method
)
write.csv(correlation_results, file.path(results_dir, correlation_output), row.names = FALSE)

filtered_matrix <- make_filtered_matrix(
  correlation_results = correlation_results,
  taxa = selected_taxa$taxon,
  environmental_variables = environmental_variables,
  correlation_threshold = correlation_threshold,
  adjusted_p_threshold = adjusted_p_threshold
)
write.csv(data.frame(taxon = rownames(filtered_matrix), filtered_matrix, check.names = FALSE), file.path(results_dir, correlation_matrix_output), row.names = FALSE)

edge_list <- correlation_results[
  abs(correlation_results$correlation) >= correlation_threshold &
    correlation_results$p_adjust <= adjusted_p_threshold,
  ,
  drop = FALSE
]

if (nrow(edge_list) == 0) {
  stop("No network edges passed the current thresholds. Lower correlation_threshold or adjusted_p_threshold.", call. = FALSE)
}

edge_list <- data.frame(
  source = edge_list$taxon,
  target = edge_list$environmental_variable,
  correlation = edge_list$correlation,
  abs_correlation = abs(edge_list$correlation),
  p_value = edge_list$p_value,
  p_adjust = edge_list$p_adjust,
  association = edge_list$association,
  stringsAsFactors = FALSE
)
write.csv(edge_list, file.path(results_dir, edge_list_output), row.names = FALSE)

taxon_nodes <- unique(edge_list$source)
environment_nodes <- unique(edge_list$target)
taxonomy_for_nodes <- taxonomy_table[match(taxon_nodes, taxonomy_table[[target_taxonomic_level]]), , drop = FALSE]

node_list <- rbind(
  data.frame(
    node = taxon_nodes,
    node_type = "taxon",
    display_label = taxon_nodes,
    Phylum = taxonomy_for_nodes$Phylum,
    Class = taxonomy_for_nodes$Class,
    Order = taxonomy_for_nodes$Order,
    Family = taxonomy_for_nodes$Family,
    Genus = taxonomy_for_nodes$Genus,
    stringsAsFactors = FALSE
  ),
  data.frame(
    node = environment_nodes,
    node_type = "environment",
    display_label = environment_nodes,
    Phylum = NA_character_,
    Class = NA_character_,
    Order = NA_character_,
    Family = NA_character_,
    Genus = NA_character_,
    stringsAsFactors = FALSE
  )
)

graph <- graph_from_data_frame(edge_list, directed = FALSE, vertices = node_list)
node_list$degree <- degree(graph)[node_list$node]
write.csv(node_list, file.path(results_dir, node_list_output), row.names = FALSE)

network_summary <- data.frame(
  target_taxonomic_level = target_taxonomic_level,
  top_n_taxa = top_n_taxa,
  environmental_variable_count = length(environmental_variables),
  correlation_method = correlation_method,
  correlation_threshold = correlation_threshold,
  adjusted_p_threshold = adjusted_p_threshold,
  node_count = vcount(graph),
  edge_count = ecount(graph),
  positive_edges = sum(edge_list$association == "positive"),
  negative_edges = sum(edge_list$association == "negative"),
  stringsAsFactors = FALSE
)
write.csv(network_summary, file.path(results_dir, network_summary_output), row.names = FALSE)

V(graph)$node_type <- node_list$node_type[match(V(graph)$name, node_list$node)]
V(graph)$degree <- degree(graph)
V(graph)$color <- ifelse(V(graph)$node_type == "taxon", taxon_node_color, environment_node_color)
V(graph)$size <- scale_to_range(V(graph)$degree, node_size_range)
V(graph)$label.cex <- ifelse(V(graph)$node_type == "environment", 0.85, 0.65)
V(graph)$label.color <- "black"

E(graph)$color <- ifelse(E(graph)$association == "positive", positive_edge_color, negative_edge_color)
E(graph)$width <- scale_to_range(E(graph)$abs_correlation, edge_width_range)

set.seed(layout_seed)
layout_matrix <- layout_with_fr(graph, weights = E(graph)$abs_correlation)

pdf(file.path(figures_dir, network_pdf), width = figure_width, height = figure_height)
par(mar = c(0, 0, 3, 0), family = "sans")
plot(
  graph,
  layout = layout_matrix,
  vertex.frame.color = "white",
  vertex.label = V(graph)$name,
  vertex.label.cex = V(graph)$label.cex,
  vertex.label.color = V(graph)$label.color,
  edge.curved = 0.08,
  main = "Microbe-environment association network"
)
legend(
  "bottomleft",
  legend = c("Taxon", "Environmental variable", "Positive association", "Negative association"),
  pch = c(21, 21, NA, NA),
  pt.bg = c(taxon_node_color, environment_node_color, NA, NA),
  col = c("black", "black", positive_edge_color, negative_edge_color),
  lty = c(NA, NA, 1, 1),
  bty = "n",
  cex = 0.8
)
dev.off()

png(file.path(figures_dir, network_png), width = figure_width, height = figure_height, units = "in", res = png_dpi)
par(mar = c(0, 0, 3, 0), family = "sans")
plot(
  graph,
  layout = layout_matrix,
  vertex.frame.color = "white",
  vertex.label = V(graph)$name,
  vertex.label.cex = V(graph)$label.cex,
  vertex.label.color = V(graph)$label.color,
  edge.curved = 0.08,
  main = "Microbe-environment association network"
)
legend(
  "bottomleft",
  legend = c("Taxon", "Environmental variable", "Positive association", "Negative association"),
  pch = c(21, 21, NA, NA),
  pt.bg = c(taxon_node_color, environment_node_color, NA, NA),
  col = c("black", "black", positive_edge_color, negative_edge_color),
  lty = c(NA, NA, 1, 1),
  bty = "n",
  cex = 0.8
)
dev.off()

message("Microbe-environment network demo completed.")
message("Selected taxa: ", nrow(selected_taxa))
message("Network nodes: ", vcount(graph))
message("Network edges: ", ecount(graph))
message("Outputs written to results/ and figures/.")
