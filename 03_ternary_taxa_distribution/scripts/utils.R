check_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop(
      "Missing input file(s):\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}

check_required_columns <- function(data, required_columns, table_name) {
  missing <- setdiff(required_columns, colnames(data))
  if (length(missing) > 0) {
    stop(
      table_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

check_sample_alignment <- function(metadata, abundance) {
  metadata_ids <- metadata$sample_id
  abundance_ids <- colnames(abundance)[colnames(abundance) != "feature_id"]

  missing_in_metadata <- setdiff(abundance_ids, metadata_ids)
  missing_in_abundance <- setdiff(metadata_ids, abundance_ids)

  if (length(missing_in_metadata) > 0 || length(missing_in_abundance) > 0) {
    message <- c()
    if (length(missing_in_metadata) > 0) {
      message <- c(message, paste0("Sample columns missing from metadata: ", paste(missing_in_metadata, collapse = ", ")))
    }
    if (length(missing_in_abundance) > 0) {
      message <- c(message, paste0("Metadata samples missing from abundance table: ", paste(missing_in_abundance, collapse = ", ")))
    }
    stop(paste(message, collapse = "\n"), call. = FALSE)
  }

  invisible(TRUE)
}

check_feature_alignment <- function(abundance, taxonomy) {
  abundance_features <- abundance$feature_id
  taxonomy_features <- taxonomy$feature_id

  missing_in_taxonomy <- setdiff(abundance_features, taxonomy_features)
  if (length(missing_in_taxonomy) > 0) {
    stop(
      "Some abundance feature_id values are missing from taxonomy_table.csv: ",
      paste(head(missing_in_taxonomy, 10), collapse = ", "),
      ifelse(length(missing_in_taxonomy) > 10, " ...", ""),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

calculate_relative_abundance <- function(abundance_table) {
  feature_ids <- abundance_table$feature_id
  count_matrix <- as.matrix(abundance_table[, colnames(abundance_table) != "feature_id", drop = FALSE])
  storage.mode(count_matrix) <- "numeric"

  sample_totals <- colSums(count_matrix, na.rm = TRUE)
  if (any(sample_totals <= 0)) {
    stop("All abundance sample columns must have positive totals.", call. = FALSE)
  }

  relative_matrix <- sweep(count_matrix, 2, sample_totals, FUN = "/")
  data.frame(
    feature_id = feature_ids,
    relative_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

aggregate_by_taxonomy <- function(relative_abundance, taxonomy, target_taxonomic_level) {
  taxonomy_subset <- taxonomy[, c("feature_id", target_taxonomic_level), drop = FALSE]
  colnames(taxonomy_subset)[2] <- "taxon"
  taxonomy_subset$taxon[is.na(taxonomy_subset$taxon) | taxonomy_subset$taxon == ""] <- "Unclassified"

  merged <- merge(taxonomy_subset, relative_abundance, by = "feature_id", all.y = TRUE, sort = FALSE)
  sample_columns <- setdiff(colnames(relative_abundance), "feature_id")
  aggregated_matrix <- rowsum(as.matrix(merged[, sample_columns, drop = FALSE]), group = merged$taxon, reorder = FALSE)

  data.frame(
    taxon = rownames(aggregated_matrix),
    aggregated_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

calculate_group_means <- function(taxon_abundance, metadata) {
  sample_columns <- setdiff(colnames(taxon_abundance), "taxon")
  groups <- unique(metadata$group)

  result <- data.frame(taxon = taxon_abundance$taxon, stringsAsFactors = FALSE)
  for (group_name in groups) {
    group_samples <- metadata$sample_id[metadata$group == group_name]
    group_samples <- intersect(group_samples, sample_columns)
    result[[group_name]] <- rowMeans(taxon_abundance[, group_samples, drop = FALSE], na.rm = TRUE)
  }

  result$overall_mean_abundance <- rowMeans(taxon_abundance[, sample_columns, drop = FALSE], na.rm = TRUE)
  result
}

make_ternary_table <- function(group_means, axis_groups, top_n) {
  required <- c("taxon", axis_groups, "overall_mean_abundance")
  missing <- setdiff(required, colnames(group_means))
  if (length(missing) > 0) {
    stop("Group mean table is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  ranked <- group_means[order(group_means$overall_mean_abundance, decreasing = TRUE), , drop = FALSE]
  top_taxa <- head(ranked, top_n)
  axis_values <- as.matrix(top_taxa[, axis_groups, drop = FALSE])
  storage.mode(axis_values) <- "numeric"
  axis_totals <- rowSums(axis_values, na.rm = TRUE)

  keep <- axis_totals > 0
  top_taxa <- top_taxa[keep, , drop = FALSE]
  axis_values <- axis_values[keep, , drop = FALSE]
  axis_totals <- axis_totals[keep]

  proportions <- sweep(axis_values, 1, axis_totals, FUN = "/")
  colnames(proportions) <- c("axis_1_proportion", "axis_2_proportion", "axis_3_proportion")

  data.frame(
    taxon = top_taxa$taxon,
    axis_1_group = axis_groups[1],
    axis_2_group = axis_groups[2],
    axis_3_group = axis_groups[3],
    axis_1_mean_abundance = axis_values[, 1],
    axis_2_mean_abundance = axis_values[, 2],
    axis_3_mean_abundance = axis_values[, 3],
    proportions,
    overall_mean_abundance = top_taxa$overall_mean_abundance,
    overall_mean_percent = top_taxa$overall_mean_abundance * 100,
    stringsAsFactors = FALSE
  )
}
