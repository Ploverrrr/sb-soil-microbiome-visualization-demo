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

check_sample_ids <- function(metadata, functional_annotation) {
  missing_in_metadata <- setdiff(unique(functional_annotation$sample_id), metadata$sample_id)
  missing_in_function <- setdiff(metadata$sample_id, unique(functional_annotation$sample_id))

  problems <- c()
  if (length(missing_in_metadata) > 0) {
    problems <- c(problems, paste0("Functional sample_id values missing from metadata: ", paste(missing_in_metadata, collapse = ", ")))
  }
  if (length(missing_in_function) > 0) {
    problems <- c(problems, paste0("Metadata sample_id values missing from functional table: ", paste(missing_in_function, collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

aggregate_function_counts <- function(functional_annotation, function_group_column) {
  data <- data.frame(
    sample_id = functional_annotation$sample_id,
    function_name = functional_annotation[[function_group_column]],
    count = functional_annotation$count,
    stringsAsFactors = FALSE
  )
  data$function_name[is.na(data$function_name) | data$function_name == ""] <- "Unannotated"
  aggregate(count ~ sample_id + function_name, data = data, FUN = sum)
}

make_relative_function_matrix <- function(aggregated_counts) {
  sample_ids <- unique(aggregated_counts$sample_id)
  function_names <- unique(aggregated_counts$function_name)
  count_matrix <- matrix(0, nrow = length(function_names), ncol = length(sample_ids), dimnames = list(function_names, sample_ids))
  count_matrix[cbind(aggregated_counts$function_name, aggregated_counts$sample_id)] <- aggregated_counts$count
  sample_totals <- colSums(count_matrix)
  if (any(sample_totals <= 0)) stop("All samples must have positive functional count totals.", call. = FALSE)
  relative_matrix <- sweep(count_matrix, 2, sample_totals, FUN = "/")
  data.frame(function_name = rownames(relative_matrix), relative_matrix, check.names = FALSE, stringsAsFactors = FALSE)
}

select_top_functions <- function(relative_abundance, top_n_functions) {
  sample_columns <- setdiff(colnames(relative_abundance), "function_name")
  relative_abundance$overall_mean_abundance <- rowMeans(relative_abundance[, sample_columns, drop = FALSE], na.rm = TRUE)
  ranked <- relative_abundance[order(relative_abundance$overall_mean_abundance, decreasing = TRUE), , drop = FALSE]
  head(ranked, top_n_functions)
}

long_function_table <- function(relative_abundance) {
  sample_columns <- setdiff(colnames(relative_abundance), c("function_name", "overall_mean_abundance"))
  rows <- vector("list", length(sample_columns))
  for (i in seq_along(sample_columns)) {
    sample_id <- sample_columns[i]
    rows[[i]] <- data.frame(
      function_name = relative_abundance$function_name,
      sample_id = sample_id,
      relative_abundance = relative_abundance[[sample_id]],
      overall_mean_abundance = relative_abundance$overall_mean_abundance,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

summarize_by_group <- function(long_data, metadata) {
  if ("group" %in% colnames(long_data)) {
    merged <- long_data
  } else {
    merged <- merge(long_data, metadata[, c("sample_id", "group")], by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  mean_table <- aggregate(relative_abundance ~ function_name + group, data = merged, FUN = mean)
  sd_table <- aggregate(relative_abundance ~ function_name + group, data = merged, FUN = sd)
  n_table <- aggregate(relative_abundance ~ function_name + group, data = merged, FUN = length)
  colnames(mean_table)[3] <- "mean_relative_abundance"
  colnames(sd_table)[3] <- "sd_relative_abundance"
  colnames(n_table)[3] <- "n_samples"
  summary <- merge(mean_table, sd_table, by = c("function_name", "group"), sort = FALSE)
  summary <- merge(summary, n_table, by = c("function_name", "group"), sort = FALSE)
  summary$se_relative_abundance <- summary$sd_relative_abundance / sqrt(summary$n_samples)
  summary
}
