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
    stop(
      table_name,
      " is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
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

  if (length(problems) > 0) {
    stop(paste(problems, collapse = "\n"), call. = FALSE)
  }

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
  if (any(sample_totals <= 0)) {
    stop("All abundance sample columns must have positive totals.", call. = FALSE)
  }

  relative_matrix <- sweep(count_matrix, 2, sample_totals, FUN = "/")
  data.frame(feature_id = feature_ids, relative_matrix, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_taxonomy_abundance <- function(relative_abundance, taxonomy, target_taxonomic_level) {
  taxonomy_subset <- taxonomy[, c("feature_id", target_taxonomic_level), drop = FALSE]
  colnames(taxonomy_subset)[2] <- "feature"
  taxonomy_subset$feature[is.na(taxonomy_subset$feature) | taxonomy_subset$feature == ""] <- "Unclassified"

  merged <- merge(taxonomy_subset, relative_abundance, by = "feature_id", all.y = TRUE, sort = FALSE)
  sample_columns <- setdiff(colnames(relative_abundance), "feature_id")
  aggregated <- rowsum(as.matrix(merged[, sample_columns, drop = FALSE]), group = merged$feature, reorder = FALSE)

  data.frame(feature = rownames(aggregated), aggregated, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_function_abundance <- function(functional_annotation, feature_column) {
  check_required_columns(
    functional_annotation,
    c("sample_id", feature_column, "count"),
    "functional_annotation_table.csv"
  )

  functional_annotation[[feature_column]][is.na(functional_annotation[[feature_column]]) | functional_annotation[[feature_column]] == ""] <- "Unannotated"
  grouped <- aggregate(
    count ~ sample_id + feature,
    data = data.frame(
      sample_id = functional_annotation$sample_id,
      feature = functional_annotation[[feature_column]],
      count = functional_annotation$count,
      stringsAsFactors = FALSE
    ),
    FUN = sum
  )

  sample_ids <- unique(grouped$sample_id)
  features <- unique(grouped$feature)
  wide <- matrix(0, nrow = length(features), ncol = length(sample_ids), dimnames = list(features, sample_ids))
  wide[cbind(grouped$feature, grouped$sample_id)] <- grouped$count

  sample_totals <- colSums(wide)
  relative <- sweep(wide, 2, sample_totals, FUN = "/")
  data.frame(feature = rownames(relative), relative, check.names = FALSE, stringsAsFactors = FALSE)
}

select_top_features <- function(feature_abundance, top_n_features) {
  sample_columns <- setdiff(colnames(feature_abundance), "feature")
  feature_abundance$overall_mean_abundance <- rowMeans(feature_abundance[, sample_columns, drop = FALSE], na.rm = TRUE)
  ranked <- feature_abundance[order(feature_abundance$overall_mean_abundance, decreasing = TRUE), , drop = FALSE]
  head(ranked, top_n_features)
}

make_correlation_table <- function(feature_abundance, environmental_data, environmental_variables, correlation_method, p_adjust_method) {
  sample_columns <- setdiff(colnames(feature_abundance), c("feature", "overall_mean_abundance"))
  feature_matrix <- as.matrix(feature_abundance[, sample_columns, drop = FALSE])
  rownames(feature_matrix) <- feature_abundance$feature

  rows <- vector("list", nrow(feature_matrix) * length(environmental_variables))
  index <- 1
  for (feature_name in rownames(feature_matrix)) {
    for (env_name in environmental_variables) {
      feature_values <- as.numeric(feature_matrix[feature_name, environmental_data$sample_id])
      env_values <- environmental_data[[env_name]]
      test <- suppressWarnings(cor.test(feature_values, env_values, method = correlation_method, exact = FALSE))
      rows[[index]] <- data.frame(
        feature = feature_name,
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
  result$significance <- ifelse(
    result$p_adjust < 0.001, "***",
    ifelse(result$p_adjust < 0.01, "**", ifelse(result$p_adjust < 0.05, "*", ""))
  )
  result
}

scale_to_range <- function(x, output_range) {
  if (length(x) == 0) return(numeric(0))
  if (diff(range(x, na.rm = TRUE)) == 0) return(rep(mean(output_range), length(x)))
  output_range[1] + (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE)) * diff(output_range)
}
