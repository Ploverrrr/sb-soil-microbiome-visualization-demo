# Purpose:
#   Create a shared simulated toy dataset for all demo modules in this
#   environmental microbiome visualization portfolio.
#
# Output directory:
#   data/toy_shared/
#
# Output files:
#   - sample_metadata.csv
#   - environmental_variables.csv
#   - taxonomy_table.csv
#   - abundance_table.csv
#   - functional_annotation_table.csv
#
# Notes:
#   This script uses simulated/desensitized toy data only. It does not read
#   private original files and does not use real sample measurements.

set.seed(20260522)

output_dir <- file.path("data", "toy_shared")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

groups <- c("Control", "Tailing", "Mining", "Smelting")
replicates_per_group <- 6
sample_metadata <- do.call(
  rbind,
  lapply(groups, function(group_name) {
    replicate_id <- seq_len(replicates_per_group)
    data.frame(
      sample_id = paste0(substr(group_name, 1, 1), sprintf("%02d", replicate_id)),
      group = group_name,
      site_type = ifelse(group_name == "Control", "reference", "metal_contaminated"),
      contamination_level = switch(
        group_name,
        "Control" = "reference",
        "Tailing" = "moderate",
        "Mining" = "high",
        "Smelting" = "high"
      ),
      replicate = replicate_id,
      batch = rep(c("Batch_A", "Batch_B"), length.out = replicates_per_group),
      stringsAsFactors = FALSE
    )
  })
)
rownames(sample_metadata) <- NULL

sample_ids <- sample_metadata$sample_id
group_effect <- c(Control = 0, Tailing = 1, Mining = 2, Smelting = 3)
group_numeric <- unname(group_effect[sample_metadata$group])

bounded_normal <- function(mean, sd, min_value = -Inf, max_value = Inf, n = length(mean)) {
  values <- rnorm(n, mean = mean, sd = sd)
  pmin(pmax(values, min_value), max_value)
}

environmental_variables <- data.frame(
  sample_id = sample_ids,
  pH = bounded_normal(6.8 - 0.28 * group_numeric, 0.18, 4.5, 8.2),
  EC = bounded_normal(120 + 95 * group_numeric, 22, 50, 520),
  TOC = bounded_normal(34 - 3.2 * group_numeric + ifelse(sample_metadata$group == "Tailing", 3, 0), 3.4, 10, 45),
  TN = bounded_normal(2.8 - 0.18 * group_numeric, 0.22, 1.2, 4.0),
  TP = bounded_normal(0.72 + 0.06 * group_numeric, 0.07, 0.3, 1.2),
  SO4 = bounded_normal(58 + 52 * group_numeric + ifelse(sample_metadata$group == "Smelting", 35, 0), 14, 20, 300),
  NO3 = bounded_normal(24 - 2.5 * group_numeric, 4.2, 4, 45),
  Sb = bounded_normal(1.5 + c(0, 8, 18, 25)[group_numeric + 1], 2.1, 0.2, 45),
  Sb3 = bounded_normal(0.35 + c(0, 1.8, 4.3, 7.0)[group_numeric + 1], 0.7, 0.05, 15),
  Sb5 = bounded_normal(1.1 + c(0, 5.0, 11.0, 15.0)[group_numeric + 1], 1.6, 0.1, 28),
  As = bounded_normal(3 + c(0, 6, 12, 15)[group_numeric + 1], 1.8, 0.2, 30),
  Cu = bounded_normal(18 + c(0, 26, 52, 70)[group_numeric + 1], 8, 5, 120),
  Zn = bounded_normal(45 + c(0, 32, 75, 95)[group_numeric + 1], 12, 15, 180),
  Cd = bounded_normal(0.18 + c(0, 0.45, 0.9, 1.2)[group_numeric + 1], 0.14, 0.02, 2.5),
  Fe = bounded_normal(18 + c(0, 2, 6, 9)[group_numeric + 1], 2.1, 8, 35),
  Mn = bounded_normal(0.42 + c(0, 0.12, 0.28, 0.38)[group_numeric + 1], 0.08, 0.1, 1.2),
  stringsAsFactors = FALSE
)

feature_count <- 96
feature_ids <- paste0("ASV", sprintf("%03d", seq_len(feature_count)))

phyla <- c(
  "Proteobacteria", "Actinobacteriota", "Acidobacteriota", "Chloroflexi",
  "Bacteroidota", "Firmicutes", "Gemmatimonadota", "Nitrospirota"
)
classes <- c(
  "Alphaproteobacteria", "Gammaproteobacteria", "Actinobacteria",
  "Acidobacteriae", "Anaerolineae", "Bacteroidia", "Bacilli",
  "Gemmatimonadetes", "Nitrospiria"
)
orders <- c(
  "Rhizobiales", "Burkholderiales", "Streptomycetales", "Acidobacteriales",
  "Anaerolineales", "Flavobacteriales", "Bacillales", "Gemmatimonadales",
  "Nitrospirales"
)
families <- c(
  "Rhizobiaceae", "Comamonadaceae", "Streptomycetaceae", "Acidobacteriaceae",
  "Anaerolineaceae", "Flavobacteriaceae", "Bacillaceae",
  "Gemmatimonadaceae", "Nitrospiraceae"
)
genus_pool <- c(
  "Aciditolerans", "Metallibacter", "Sulfuritalea", "Nitrospira",
  "Rhizomicrobium", "Terrimonas", "Cupriavidus", "Pseudarthrobacter",
  "Gemmatimonas", "Sideroxydans", "Bacillus", "Flavobacterium",
  "Anaerolinea", "Streptomyces", "Bradyrhizobium", "Thiobacillus"
)

taxonomy_table <- data.frame(
  feature_id = feature_ids,
  Kingdom = "Bacteria",
  Phylum = rep(phyla, length.out = feature_count),
  Class = rep(classes, length.out = feature_count),
  Order = rep(orders, length.out = feature_count),
  Family = rep(families, length.out = feature_count),
  Genus = rep(genus_pool, length.out = feature_count),
  stringsAsFactors = FALSE
)
taxonomy_table$taxon_label <- paste0(taxonomy_table$Genus, "_", sprintf("%02d", ave(seq_len(feature_count), taxonomy_table$Genus, FUN = seq_along)))

taxon_pattern <- rep(c("control_enriched", "tailing_enriched", "mining_enriched", "smelting_enriched", "broad", "rare"), length.out = feature_count)
names(taxon_pattern) <- feature_ids

base_means <- runif(feature_count, min = 55, max = 520)
abundance_matrix <- matrix(0, nrow = feature_count, ncol = length(sample_ids))
rownames(abundance_matrix) <- feature_ids
colnames(abundance_matrix) <- sample_ids

for (feature_index in seq_len(feature_count)) {
  pattern <- taxon_pattern[feature_index]
  for (sample_index in seq_along(sample_ids)) {
    group_name <- sample_metadata$group[sample_index]
    multiplier <- switch(
      pattern,
      "control_enriched" = ifelse(group_name == "Control", 4.2, 0.75),
      "tailing_enriched" = ifelse(group_name == "Tailing", 4.0, 0.8),
      "mining_enriched" = ifelse(group_name == "Mining", 4.4, 0.7),
      "smelting_enriched" = ifelse(group_name == "Smelting", 4.6, 0.65),
      "broad" = 1.4,
      "rare" = 0.22,
      1
    )
    if (pattern == "broad" && group_name != "Control") {
      multiplier <- multiplier + 0.1 * group_numeric[sample_index]
    }
    lambda <- base_means[feature_index] * multiplier * runif(1, 0.82, 1.18)
    abundance_matrix[feature_index, sample_index] <- rpois(1, lambda = max(lambda, 1))
  }
}

abundance_table <- data.frame(
  feature_id = rownames(abundance_matrix),
  abundance_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

annotation_sources <- c("KEGG", "FAPROTAX", "BacMet", "SulfurCycle", "NitrogenCycle", "CarbonCycle")
function_catalog <- data.frame(
  annotation_source = annotation_sources,
  function_category = c(
    "Energy metabolism", "Element cycling", "Metal resistance",
    "Sulfur cycling", "Nitrogen cycling", "Carbon cycling"
  ),
  pathway = c(
    "Oxidative phosphorylation", "Aerobic chemoheterotrophy",
    "Metal efflux and detoxification", "Dissimilatory sulfur metabolism",
    "Nitrification and denitrification", "Organic carbon degradation"
  ),
  ko_id = c("K02111", "FAP001", "BMR001", "K11180", "K00370", "K01190"),
  stringsAsFactors = FALSE
)

gene_ids <- paste0("gene_", sprintf("%03d", seq_len(72)))
gene_feature_map <- data.frame(
  gene_id = gene_ids,
  feature_id = sample(feature_ids, length(gene_ids), replace = TRUE),
  catalog_index = rep(seq_len(nrow(function_catalog)), length.out = length(gene_ids)),
  stringsAsFactors = FALSE
)

functional_rows <- vector("list", length(gene_ids) * length(sample_ids))
row_index <- 1
for (gene_index in seq_along(gene_ids)) {
  catalog_row <- function_catalog[gene_feature_map$catalog_index[gene_index], ]
  linked_feature <- gene_feature_map$feature_id[gene_index]
  pattern <- taxon_pattern[linked_feature]
  base_count <- runif(1, 18, 95)
  for (sample_index in seq_along(sample_ids)) {
    group_name <- sample_metadata$group[sample_index]
    contamination <- group_numeric[sample_index]
    source <- catalog_row$annotation_source
    multiplier <- 1
    if (source %in% c("BacMet", "SulfurCycle")) {
      multiplier <- 1 + 0.65 * contamination
    } else if (source == "FAPROTAX" && pattern %in% c("mining_enriched", "smelting_enriched")) {
      multiplier <- ifelse(group_name %in% c("Mining", "Smelting"), 2.3, 0.9)
    } else if (source == "NitrogenCycle") {
      multiplier <- ifelse(group_name == "Control", 1.6, 1.0)
    } else if (source == "CarbonCycle") {
      multiplier <- ifelse(group_name == "Control", 1.8, 1.1)
    }
    count_value <- rpois(1, lambda = max(base_count * multiplier * runif(1, 0.75, 1.25), 1))
    functional_rows[[row_index]] <- data.frame(
      gene_id = gene_ids[gene_index],
      feature_id = linked_feature,
      sample_id = sample_ids[sample_index],
      annotation_source = source,
      function_category = catalog_row$function_category,
      pathway = catalog_row$pathway,
      ko_id = catalog_row$ko_id,
      count = count_value,
      abundance = NA_real_,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1
  }
}

functional_annotation_table <- do.call(rbind, functional_rows)
sample_totals <- ave(functional_annotation_table$count, functional_annotation_table$sample_id, FUN = sum)
functional_annotation_table$abundance <- round(functional_annotation_table$count / sample_totals, 8)

write.csv(sample_metadata, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)
write.csv(environmental_variables, file.path(output_dir, "environmental_variables.csv"), row.names = FALSE)
write.csv(taxonomy_table, file.path(output_dir, "taxonomy_table.csv"), row.names = FALSE)
write.csv(abundance_table, file.path(output_dir, "abundance_table.csv"), row.names = FALSE)
write.csv(functional_annotation_table, file.path(output_dir, "functional_annotation_table.csv"), row.names = FALSE)

message("Shared toy dataset created in: ", output_dir)
message("sample_metadata: ", nrow(sample_metadata), " rows x ", ncol(sample_metadata), " columns")
message("environmental_variables: ", nrow(environmental_variables), " rows x ", ncol(environmental_variables), " columns")
message("taxonomy_table: ", nrow(taxonomy_table), " rows x ", ncol(taxonomy_table), " columns")
message("abundance_table: ", nrow(abundance_table), " rows x ", ncol(abundance_table), " columns")
message("functional_annotation_table: ", nrow(functional_annotation_table), " rows x ", ncol(functional_annotation_table), " columns")
