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

read_abundance_matrix <- function(abundance_file) {
  abundance <- read.csv(abundance_file, stringsAsFactors = FALSE, check.names = FALSE)
  check_required_columns(abundance, "feature_id", "abundance_table.csv")

  matrix_data <- as.matrix(abundance[, setdiff(colnames(abundance), "feature_id"), drop = FALSE])
  storage.mode(matrix_data) <- "numeric"
  rownames(matrix_data) <- abundance$feature_id

  if (anyNA(matrix_data)) stop("abundance_table.csv contains non-numeric or missing values.", call. = FALSE)
  if (any(matrix_data < 0)) stop("abundance_table.csv contains negative values.", call. = FALSE)

  matrix_data
}

make_taxonomy_response_matrix <- function(abundance_matrix, metadata, top_n_features) {
  sample_ids <- metadata$sample_id
  missing_samples <- setdiff(sample_ids, colnames(abundance_matrix))
  if (length(missing_samples) > 0) {
    stop("Sample(s) missing from abundance_table.csv: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  abundance_matrix <- abundance_matrix[, sample_ids, drop = FALSE]
  feature_totals <- rowSums(abundance_matrix)
  top_features <- names(sort(feature_totals, decreasing = TRUE))[seq_len(min(top_n_features, length(feature_totals)))]
  response <- t(abundance_matrix[top_features, , drop = FALSE])
  response[, colSums(response) > 0, drop = FALSE]
}

make_function_response_matrix <- function(functional_annotation, metadata, value_column, feature_column, top_n_features) {
  check_required_columns(
    functional_annotation,
    c("sample_id", feature_column, value_column),
    "functional_annotation_table.csv"
  )

  sample_ids <- metadata$sample_id
  subset_data <- functional_annotation[functional_annotation$sample_id %in% sample_ids, , drop = FALSE]
  aggregated <- aggregate(
    subset_data[[value_column]],
    by = list(sample_id = subset_data$sample_id, feature = subset_data[[feature_column]]),
    FUN = sum
  )
  colnames(aggregated)[3] <- value_column

  wide <- reshape(aggregated, idvar = "sample_id", timevar = "feature", direction = "wide")
  colnames(wide) <- sub(paste0("^", value_column, "[.]"), "", colnames(wide))
  rownames(wide) <- wide$sample_id
  wide <- wide[sample_ids, setdiff(colnames(wide), "sample_id"), drop = FALSE]
  wide[is.na(wide)] <- 0

  feature_totals <- colSums(wide)
  top_features <- names(sort(feature_totals, decreasing = TRUE))[seq_len(min(top_n_features, length(feature_totals)))]
  as.matrix(wide[, top_features, drop = FALSE])
}

align_environment <- function(environmental_variables, metadata, environmental_vars) {
  check_required_columns(environmental_variables, c("sample_id", environmental_vars), "environmental_variables.csv")
  missing_samples <- setdiff(metadata$sample_id, environmental_variables$sample_id)
  if (length(missing_samples) > 0) {
    stop("Sample(s) missing from environmental_variables.csv: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }

  env <- environmental_variables[match(metadata$sample_id, environmental_variables$sample_id), c("sample_id", environmental_vars), drop = FALSE]
  rownames(env) <- env$sample_id
  env <- env[, environmental_vars, drop = FALSE]
  env[] <- lapply(env, function(x) as.numeric(as.character(x)))
  if (anyNA(env)) stop("Configured environmental variables contain missing or non-numeric values.", call. = FALSE)
  env
}

clean_vpa_fractions <- function(vpa_obj, response_name) {
  fractions <- vpa_obj$part$indfract
  data.frame(
    response = response_name,
    fraction = rownames(fractions),
    adjusted_r2 = fractions$Adj.R.square,
    cleaned_adjusted_r2 = pmax(fractions$Adj.R.square, 0),
    stringsAsFactors = FALSE
  )
}

run_partial_rda_test <- function(response_matrix, focal_env, covariates, permutations, seed) {
  set.seed(seed)
  test_data <- data.frame(focal = focal_env[, 1], covariates, check.names = FALSE)
  response_hellinger <- vegan::decostand(response_matrix, method = "hellinger")
  model <- vegan::rda(response_hellinger ~ focal + Condition(as.matrix(test_data[, setdiff(colnames(test_data), "focal"), drop = FALSE])), data = test_data)
  result <- vegan::anova.cca(model, permutations = permutations)
  as.data.frame(result)
}

save_base_vpa_pair <- function(vpa_species, vpa_genes, pdf_file, png_file, width, height, colors, group_names) {
  pdf(pdf_file, width = width, height = height, family = "serif")
  par(mfrow = c(1, 2), mar = c(2, 2, 4, 2), family = "serif")
  plot(vpa_species, bg = colors, Xnames = group_names, cutoff = 0, digits = 2, cex = 1.25)
  title(main = "(A) Species Community", cex.main = 1.6, font.main = 2)
  plot(vpa_genes, bg = colors, Xnames = group_names, cutoff = 0, digits = 2, cex = 1.25)
  title(main = "(B) Functional Genes", cex.main = 1.6, font.main = 2)
  dev.off()

  png(png_file, width = width, height = height, units = "in", res = 300, bg = "white")
  par(mfrow = c(1, 2), mar = c(2, 2, 4, 2), family = "serif")
  plot(vpa_species, bg = colors, Xnames = group_names, cutoff = 0, digits = 2, cex = 1.25)
  title(main = "(A) Species Community", cex.main = 1.6, font.main = 2)
  plot(vpa_genes, bg = colors, Xnames = group_names, cutoff = 0, digits = 2, cex = 1.25)
  title(main = "(B) Functional Genes", cex.main = 1.6, font.main = 2)
  dev.off()

  invisible(TRUE)
}

save_ggplot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(TRUE)
}
