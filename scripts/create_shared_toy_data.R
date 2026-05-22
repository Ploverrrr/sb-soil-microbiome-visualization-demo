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
        "Smelting" = "very_high"
      ),
      replicate = replicate_id,
      batch = rep(c("Batch_A", "Batch_B"), length.out = replicates_per_group),
      stringsAsFactors = FALSE
    )
  })
)
rownames(sample_metadata) <- NULL

sample_ids <- sample_metadata$sample_id
group_numeric <- c(Control = 0, Tailing = 1, Mining = 2, Smelting = 3)[sample_metadata$group]

bounded_normal <- function(mean, sd, min_value = -Inf, max_value = Inf) {
  values <- rnorm(length(mean), mean = mean, sd = sd)
  pmin(pmax(values, min_value), max_value)
}

group_mean <- function(control, tailing, mining, smelting) {
  values <- c(Control = control, Tailing = tailing, Mining = mining, Smelting = smelting)
  unname(values[sample_metadata$group])
}

sb_iii <- bounded_normal(group_mean(0.18, 1.85, 4.9, 7.3), 0.35, 0.03, 12)
sb_v <- bounded_normal(group_mean(0.95, 5.8, 12.5, 17.2), 0.85, 0.08, 26)
sb_total <- (sb_iii + sb_v) * bounded_normal(rep(1.015, length(sample_ids)), 0.025, 0.96, 1.08)

environmental_variables <- data.frame(
  sample_id = sample_ids,
  pH = bounded_normal(group_mean(6.82, 6.45, 6.05, 5.82), 0.16, 4.5, 8.2),
  EC = bounded_normal(group_mean(118, 205, 340, 390), 24, 60, 520),
  TOC = bounded_normal(group_mean(35.5, 33.0, 27.0, 24.5), 2.8, 12, 45),
  TN = bounded_normal(group_mean(2.85, 2.65, 2.25, 2.05), 0.18, 1.1, 4.0),
  TP = bounded_normal(group_mean(0.68, 0.78, 0.88, 0.94), 0.06, 0.3, 1.3),
  SO4 = bounded_normal(group_mean(55, 118, 190, 245), 18, 20, 330),
  NO3 = bounded_normal(group_mean(27, 23, 18, 16), 3.5, 4, 45),
  Sb_total = round(sb_total, 4),
  Sb_III = round(sb_iii, 4),
  Sb_V = round(sb_v, 4),
  As = bounded_normal(group_mean(3.2, 8.7, 15.8, 19.5), 1.5, 0.2, 32),
  Cu = bounded_normal(group_mean(18, 46, 86, 102), 7.5, 5, 140),
  Zn = bounded_normal(group_mean(45, 76, 128, 150), 11, 15, 210),
  Cd = bounded_normal(group_mean(0.16, 0.55, 1.05, 1.35), 0.12, 0.02, 2.8),
  Fe = bounded_normal(group_mean(18.0, 20.5, 24.0, 27.5), 1.9, 8, 36),
  Mn = bounded_normal(group_mean(0.40, 0.55, 0.74, 0.84), 0.07, 0.1, 1.3),
  stringsAsFactors = FALSE
)

taxonomy_lineages <- data.frame(
  Phylum = c(
    "Pseudomonadota", "Pseudomonadota", "Pseudomonadota", "Pseudomonadota",
    "Pseudomonadota", "Pseudomonadota", "Pseudomonadota", "Pseudomonadota",
    "Actinobacteriota", "Actinobacteriota", "Actinobacteriota", "Actinobacteriota",
    "Actinobacteriota", "Actinobacteriota",
    "Acidobacteriota", "Acidobacteriota", "Acidobacteriota", "Acidobacteriota",
    "Acidobacteriota",
    "Bacteroidota", "Bacteroidota", "Bacteroidota", "Bacteroidota", "Bacteroidota",
    "Bacillota", "Bacillota", "Bacillota", "Bacillota", "Bacillota",
    "Chloroflexi", "Chloroflexi", "Chloroflexi", "Chloroflexi",
    "Gemmatimonadota", "Gemmatimonadota", "Nitrospirota", "Methylomirabilota",
    "Desulfobacterota"
  ),
  Class = c(
    "Alphaproteobacteria", "Alphaproteobacteria", "Alphaproteobacteria", "Gammaproteobacteria",
    "Gammaproteobacteria", "Gammaproteobacteria", "Betaproteobacteria", "Acidithiobacillia",
    "Actinomycetia", "Actinomycetia", "Actinomycetia", "Actinomycetia",
    "Thermoleophilia", "Thermoleophilia",
    "Acidobacteriae", "Acidobacteriae", "Blastocatellia", "Blastocatellia",
    "Vicinamibacteria",
    "Bacteroidia", "Bacteroidia", "Bacteroidia", "Chitinophagia", "Chitinophagia",
    "Bacilli", "Bacilli", "Bacilli", "Clostridia", "Clostridia",
    "Anaerolineae", "Anaerolineae", "Ktedonobacteria", "Thermomicrobia",
    "Gemmatimonadetes", "Gemmatimonadetes", "Nitrospiria", "Methylomirabilia",
    "Desulfobacteria"
  ),
  Order = c(
    "Rhizobiales", "Sphingomonadales", "Rhodospirillales", "Burkholderiales",
    "Pseudomonadales", "Xanthomonadales", "Nitrosomonadales", "Acidithiobacillales",
    "Streptomycetales", "Micrococcales", "Frankiales", "Propionibacteriales",
    "Solirubrobacterales", "Gaiellales",
    "Acidobacteriales", "Acidobacteriales", "Blastocatellales", "Pyrinomonadales",
    "Vicinamibacterales",
    "Flavobacteriales", "Sphingobacteriales", "Cytophagales", "Chitinophagales", "Chitinophagales",
    "Bacillales", "Bacillales", "Lactobacillales", "Clostridiales", "Lachnospirales",
    "Anaerolineales", "Anaerolineales", "Ktedonobacterales", "Thermomicrobiales",
    "Gemmatimonadales", "Longimicrobiales", "Nitrospirales", "Methylomirabilales",
    "Desulfobacterales"
  ),
  Family = c(
    "Rhizobiaceae", "Sphingomonadaceae", "Acetobacteraceae", "Comamonadaceae",
    "Pseudomonadaceae", "Xanthomonadaceae", "Nitrosomonadaceae", "Acidithiobacillaceae",
    "Streptomycetaceae", "Micrococcaceae", "Frankiaceae", "Nocardioidaceae",
    "Solirubrobacteraceae", "Gaiellaceae",
    "Acidobacteriaceae", "Terriglobaceae", "Blastocatellaceae", "Pyrinomonadaceae",
    "Vicinamibacteraceae",
    "Flavobacteriaceae", "Sphingobacteriaceae", "Cytophagaceae", "Chitinophagaceae", "Chitinophagaceae",
    "Bacillaceae", "Paenibacillaceae", "Lactobacillaceae", "Clostridiaceae", "Lachnospiraceae",
    "Anaerolineaceae", "Caldilineaceae", "Ktedonobacteraceae", "Thermomicrobiaceae",
    "Gemmatimonadaceae", "Longimicrobiaceae", "Nitrospiraceae", "Methylomirabilaceae",
    "Desulfobacteraceae"
  ),
  Genus = c(
    "Bradyrhizobium", "Sphingomonas", "Acidocella", "Cupriavidus",
    "Pseudomonas", "Lysobacter", "Nitrosomonas", "Acidithiobacillus",
    "Streptomyces", "Arthrobacter", "Frankia", "Nocardioides",
    "Solirubrobacter", "Gaiella",
    "Acidipila", "Terriglobus", "Blastocatella", "Pyrinomonas",
    "Vicinamibacter",
    "Flavobacterium", "Pedobacter", "Cytophaga", "Terrimonas", "Segetibacter",
    "Bacillus", "Paenibacillus", "Lactobacillus", "Clostridium", "Anaerostipes",
    "Anaerolinea", "Caldilinea", "Ktedonobacter", "Thermomicrobium",
    "Gemmatimonas", "Longimicrobium", "Nitrospira", "Methylomirabilis",
    "Desulfobacter"
  ),
  stringsAsFactors = FALSE
)

features_per_lineage <- 4
feature_count <- nrow(taxonomy_lineages) * features_per_lineage
feature_ids <- paste0("ASV", sprintf("%03d", seq_len(feature_count)))
lineage_index <- rep(seq_len(nrow(taxonomy_lineages)), each = features_per_lineage)
taxonomy_table <- cbind(
  data.frame(
    feature_id = feature_ids,
    Kingdom = "Bacteria",
    stringsAsFactors = FALSE
  ),
  taxonomy_lineages[lineage_index, ]
)
taxonomy_table$taxon_label <- paste0(
  taxonomy_table$Genus,
  "_",
  sprintf("%02d", ave(seq_len(feature_count), taxonomy_table$Genus, FUN = seq_along))
)
rownames(taxonomy_table) <- NULL

lineage_pattern <- rep(
  c("control_enriched", "tailing_enriched", "mining_enriched", "smelting_enriched", "broad", "rare"),
  length.out = nrow(taxonomy_lineages)
)
taxon_pattern <- lineage_pattern[lineage_index]
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
      "tailing_enriched" = ifelse(group_name == "Tailing", 4.0, ifelse(group_name == "Control", 0.95, 0.78)),
      "mining_enriched" = ifelse(group_name == "Mining", 4.4, ifelse(group_name == "Smelting", 1.2, 0.7)),
      "smelting_enriched" = ifelse(group_name == "Smelting", 4.6, ifelse(group_name == "Mining", 1.15, 0.65)),
      "broad" = 1.35 + 0.08 * unname(group_numeric[sample_index]),
      "rare" = 0.22,
      1
    )
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

function_catalog <- data.frame(
  annotation_source = c(
    "KEGG", "KEGG", "KEGG", "FAPROTAX",
    "BacMet", "BacMet", "SulfurCycle", "SulfurCycle",
    "SulfurCycle", "NitrogenCycle", "NitrogenCycle", "CarbonCycle"
  ),
  function_category = c(
    "Energy metabolism", "Membrane transport", "Stress response", "Element cycling",
    "Metal resistance", "Metal resistance", "Sulfur cycling", "Sulfur cycling",
    "Sulfur cycling", "Nitrogen cycling", "Nitrogen cycling", "Carbon cycling"
  ),
  pathway = c(
    "Oxidative phosphorylation", "ABC transporters", "Two-component stress response",
    "Aerobic chemoheterotrophy", "Copper resistance", "Antimony efflux",
    "Dissimilatory sulfite reduction", "Thiosulfate oxidation", "Sulfate assimilation",
    "Nitrification", "Denitrification", "Organic carbon degradation"
  ),
  ko_id = c(
    "K02111", "K02003", "K07636", "FAP001",
    "BMR_COP", "BMR_ANT", "K11180", "K17218",
    "K00958", "K10944", "K00370", "K01190"
  ),
  trend = c(
    "broad", "contamination_mild", "contamination_mild", "carbon_control",
    "metal_strong", "metal_strong", "sulfur_strong", "sulfur_smelting",
    "sulfur_mild", "nitrogen_control", "nitrogen_mining", "carbon_control"
  ),
  stringsAsFactors = FALSE
)

gene_ids <- paste0("gene_", sprintf("%03d", seq_len(96)))
gene_feature_map <- data.frame(
  gene_id = gene_ids,
  feature_id = sample(feature_ids, length(gene_ids), replace = TRUE),
  catalog_index = rep(seq_len(nrow(function_catalog)), length.out = length(gene_ids)),
  stringsAsFactors = FALSE
)

trend_multiplier <- function(trend, group_name, contamination_score) {
  switch(
    trend,
    "broad" = 1.15 + runif(1, -0.08, 0.08),
    "contamination_mild" = 1 + 0.25 * contamination_score,
    "metal_strong" = 0.75 + c(Control = 0.0, Tailing = 0.9, Mining = 1.9, Smelting = 2.35)[group_name],
    "sulfur_strong" = 0.85 + c(Control = 0.0, Tailing = 0.75, Mining = 1.55, Smelting = 2.05)[group_name],
    "sulfur_smelting" = ifelse(group_name == "Smelting", 3.1, ifelse(group_name == "Mining", 1.85, 0.95)),
    "sulfur_mild" = 1 + c(Control = 0.0, Tailing = 0.35, Mining = 0.75, Smelting = 0.95)[group_name],
    "nitrogen_control" = ifelse(group_name == "Control", 1.85, ifelse(group_name == "Tailing", 1.35, 0.9)),
    "nitrogen_mining" = ifelse(group_name %in% c("Tailing", "Mining"), 1.65, 1.0),
    "carbon_control" = ifelse(group_name == "Control", 1.95, ifelse(group_name == "Tailing", 1.45, 0.95)),
    1
  )
}

functional_rows <- vector("list", length(gene_ids) * length(sample_ids))
row_index <- 1
for (gene_index in seq_along(gene_ids)) {
  catalog_row <- function_catalog[gene_feature_map$catalog_index[gene_index], ]
  linked_feature <- gene_feature_map$feature_id[gene_index]
  base_count <- runif(1, 18, 95)
  for (sample_index in seq_along(sample_ids)) {
    group_name <- sample_metadata$group[sample_index]
    contamination_score <- unname(group_numeric[sample_index])
    multiplier <- trend_multiplier(catalog_row$trend, group_name, contamination_score)
    count_value <- rpois(1, lambda = max(base_count * multiplier * runif(1, 0.78, 1.22), 1))
    functional_rows[[row_index]] <- data.frame(
      gene_id = gene_ids[gene_index],
      feature_id = linked_feature,
      sample_id = sample_ids[sample_index],
      annotation_source = catalog_row$annotation_source,
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
