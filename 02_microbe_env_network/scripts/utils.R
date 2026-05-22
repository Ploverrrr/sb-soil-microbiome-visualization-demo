check_files_exist <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    stop("Missing input file(s):\n", paste(missing, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

check_required_columns <- function(data, required_columns, table_name) {
  missing <- setdiff(required_columns, colnames(data))
  if (length(missing) > 0) {
    stop(table_name, " is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

check_sample_alignment <- function(metadata, environmental, abundance) {
  metadata_ids <- metadata$sample_id
  environmental_ids <- environmental$sample_id
  abundance_ids <- setdiff(colnames(abundance), "feature_id")

  problems <- c()
  if (!setequal(metadata_ids, environmental_ids)) {
    problems <- c(problems, "sample_metadata.csv and environmental_variables.csv have different sample_id sets.")
  }
  if (!setequal(metadata_ids, abundance_ids)) {
    problems <- c(problems, "sample_metadata.csv sample_id values do not match abundance_table.csv sample columns.")
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

check_feature_alignment <- function(abundance, taxonomy) {
  missing_in_taxonomy <- setdiff(abundance$feature_id, taxonomy$feature_id)
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
  count_matrix <- as.matrix(abundance_table[, setdiff(colnames(abundance_table), "feature_id"), drop = FALSE])
  storage.mode(count_matrix) <- "numeric"
  sample_totals <- colSums(count_matrix, na.rm = TRUE)
  if (any(sample_totals <= 0)) stop("All abundance sample columns must have positive totals.", call. = FALSE)
  relative_matrix <- sweep(count_matrix, 2, sample_totals, FUN = "/")
  data.frame(feature_id = feature_ids, relative_matrix, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_by_taxonomy <- function(relative_abundance, taxonomy, target_taxonomic_level) {
  taxonomy_subset <- taxonomy[, c("feature_id", target_taxonomic_level), drop = FALSE]
  colnames(taxonomy_subset)[2] <- "taxon"
  taxonomy_subset$taxon[is.na(taxonomy_subset$taxon) | taxonomy_subset$taxon == ""] <- "Unclassified"

  merged <- merge(taxonomy_subset, relative_abundance, by = "feature_id", all.y = TRUE, sort = FALSE)
  sample_columns <- setdiff(colnames(relative_abundance), "feature_id")
  aggregated <- rowsum(as.matrix(merged[, sample_columns, drop = FALSE]), group = merged$taxon, reorder = FALSE)

  data.frame(taxon = rownames(aggregated), aggregated, check.names = FALSE, stringsAsFactors = FALSE)
}

select_top_taxa <- function(taxon_abundance, top_n_taxa) {
  sample_columns <- setdiff(colnames(taxon_abundance), "taxon")
  taxon_abundance$overall_mean_abundance <- rowMeans(taxon_abundance[, sample_columns, drop = FALSE], na.rm = TRUE)
  ranked <- taxon_abundance[order(taxon_abundance$overall_mean_abundance, decreasing = TRUE), , drop = FALSE]
  head(ranked, top_n_taxa)
}

calculate_microbe_env_correlations <- function(taxon_abundance, environmental_data, environmental_variables, correlation_method, p_adjust_method) {
  sample_columns <- setdiff(colnames(taxon_abundance), c("taxon", "overall_mean_abundance"))
  taxon_matrix <- as.matrix(taxon_abundance[, sample_columns, drop = FALSE])
  rownames(taxon_matrix) <- taxon_abundance$taxon

  rows <- vector("list", nrow(taxon_matrix) * length(environmental_variables))
  index <- 1
  for (taxon_name in rownames(taxon_matrix)) {
    for (env_name in environmental_variables) {
      microbe_values <- as.numeric(taxon_matrix[taxon_name, environmental_data$sample_id])
      env_values <- environmental_data[[env_name]]
      test <- suppressWarnings(cor.test(microbe_values, env_values, method = correlation_method, exact = FALSE))
      rows[[index]] <- data.frame(
        taxon = taxon_name,
        environmental_variable = env_name,
        correlation = unname(test$estimate),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }

  result <- do.call(rbind, rows)
  result$p_adjust <- p.adjust(result$p_value, method = p_adjust_method)
  result$association <- ifelse(result$correlation >= 0, "positive", "negative")
  result
}

make_filtered_matrix <- function(correlation_results, taxa, environmental_variables, correlation_threshold, adjusted_p_threshold) {
  filtered <- correlation_results[
    abs(correlation_results$correlation) >= correlation_threshold &
      correlation_results$p_adjust <= adjusted_p_threshold,
    ,
    drop = FALSE
  ]

  mat <- matrix(0, nrow = length(taxa), ncol = length(environmental_variables), dimnames = list(taxa, environmental_variables))
  if (nrow(filtered) > 0) {
    mat[cbind(filtered$taxon, filtered$environmental_variable)] <- filtered$correlation
  }
  mat
}

scale_to_range <- function(x, output_range) {
  if (length(x) == 0) return(numeric(0))
  if (diff(range(x, na.rm = TRUE)) == 0) return(rep(mean(output_range), length(x)))
  output_range[1] + (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * diff(output_range)
}
