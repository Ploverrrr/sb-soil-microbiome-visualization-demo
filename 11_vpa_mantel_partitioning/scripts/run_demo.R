# Purpose:
#   Build an independent VPA and Mantel partitioning demo from shared simulated
#   soil microbiome toy data. The workflow follows a reference workflow pattern:
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
mantel_link_sizes <- c("<0.2" = 0.5, "0.2-0.4" = 1, ">=0.4" = 2)
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

varechem <- env_raw[, mantel_environmental_variables, drop = FALSE]
varespec <- species_matrix

# Compatibility shim only: the legacy plotting pattern calls bare geom_square(), but
# this ggcor/ggplot2 combination needs width and height mapped for that layer.
geom_square <- function(...) {
  suppressWarnings(ggcor::geom_square(aes(width = 1, height = 1), ...))
}

mantel <- mantel_test(varespec, varechem, mantel.fun = 'mantel',
                      spec.dist.method = 'bray', env.dist.method = 'euclidean',
                      spec.select = list(callvulg = 1:2,B = 3:44), #spec.select = NULL
                      env.select = NULL) %>% # Define response blocks for Mantel testing.
  mutate(rd = cut(r,breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("<0.2","0.2-0.4",">=0.4")), # Mantel r categories for plotting.
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("<0.01","0.01-0.05",">=0.05")))  # Mantel p-value categories for plotting.

mantel_results <- mantel
write.csv(as.data.frame(mantel_results), file.path(results_dir, mantel_output), row.names = FALSE)

env_correlation <- stats::cor(varechem, method = "spearman")
write.csv(env_correlation, file.path(results_dir, env_correlation_output), row.names = TRUE)

p1 <- quickcor (varechem, method = "spearman", type = "upper")+ 
  # Draw an environmental correlation heatmap; method can also be pearson or kendall.
  geom_square()+  # Use square tiles for the correlation heatmap.
  scale_fill_gradient2(low = '#c6f093', mid = '#f8fff8', high = '#163f00')+ # Environmental correlation colors.
  anno_link(aes (colour = pd, size = rd), data = mantel) + # Mantel association links.
  scale_color_manual (values = c("#62a11b","#68edcb", "snow3")) + # Link colors by Mantel p-value.
  scale_size_manual (values = c(0.5, 1, 2))+ # Link widths by Mantel r category.
  guides (size = guide_legend(title ="Mantel's r", # Plot legends.
                              order = 2),
          colour = guide_legend (title = "Mantel's p",
                                 order = 3),
          fill = guide_colorbar (title = "Spearman's r", order = 4))

save_ggplot_pair(
  p1,
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
