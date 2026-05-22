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

check_sample_ids <- function(metadata, abundance) {
  sample_columns <- setdiff(colnames(abundance), "feature_id")
  missing_in_metadata <- setdiff(sample_columns, metadata$sample_id)
  missing_in_abundance <- setdiff(metadata$sample_id, sample_columns)

  problems <- c()
  if (length(missing_in_metadata) > 0) {
    problems <- c(problems, paste0("Abundance sample columns missing from metadata: ", paste(missing_in_metadata, collapse = ", ")))
  }
  if (length(missing_in_abundance) > 0) {
    problems <- c(problems, paste0("Metadata sample_id values missing from abundance table: ", paste(missing_in_abundance, collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

check_feature_ids <- function(abundance, taxonomy) {
  missing_in_taxonomy <- setdiff(abundance$feature_id, taxonomy$feature_id)
  missing_in_abundance <- setdiff(taxonomy$feature_id, abundance$feature_id)

  problems <- c()
  if (length(missing_in_taxonomy) > 0) {
    problems <- c(problems, paste0("feature_id values missing from taxonomy table: ", paste(head(missing_in_taxonomy, 10), collapse = ", ")))
  }
  if (length(missing_in_abundance) > 0) {
    problems <- c(problems, paste0("feature_id values missing from abundance table: ", paste(head(missing_in_abundance, 10), collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

calculate_relative_abundance <- function(abundance) {
  sample_columns <- setdiff(colnames(abundance), "feature_id")
  count_matrix <- as.matrix(abundance[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"
  sample_totals <- colSums(count_matrix, na.rm = TRUE)
  if (any(sample_totals <= 0)) stop("All samples must have positive abundance totals.", call. = FALSE)
  relative_matrix <- sweep(count_matrix, 2, sample_totals, FUN = "/")
  data.frame(feature_id = abundance$feature_id, relative_matrix, check.names = FALSE, stringsAsFactors = FALSE)
}

aggregate_to_taxon <- function(relative_abundance, taxonomy, target_taxonomic_level) {
  merged <- merge(
    taxonomy[, c("feature_id", target_taxonomic_level)],
    relative_abundance,
    by = "feature_id",
    all.y = TRUE,
    sort = FALSE
  )
  merged$taxon <- merged[[target_taxonomic_level]]
  merged$taxon[is.na(merged$taxon) | merged$taxon == ""] <- "Unclassified"
  sample_columns <- setdiff(colnames(relative_abundance), "feature_id")
  aggregate(merged[, sample_columns, drop = FALSE], by = list(taxon = merged$taxon), FUN = sum)
}

taxon_table_to_long <- function(taxon_relative_abundance) {
  sample_columns <- setdiff(colnames(taxon_relative_abundance), "taxon")
  rows <- vector("list", length(sample_columns))
  for (i in seq_along(sample_columns)) {
    sample_id <- sample_columns[i]
    rows[[i]] <- data.frame(
      taxon = taxon_relative_abundance$taxon,
      sample_id = sample_id,
      relative_abundance = taxon_relative_abundance[[sample_id]],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

standard_error <- function(x) {
  if (length(x) <= 1) return(0)
  stats::sd(x) / sqrt(length(x))
}

safe_kruskal_p <- function(values, groups) {
  if (length(unique(groups)) < 2 || length(unique(values)) < 2) return(1)
  tryCatch(stats::kruskal.test(values ~ groups)$p.value, error = function(e) 1)
}

safe_wilcox_p <- function(values, is_enriched_group) {
  if (length(unique(is_enriched_group)) < 2 || length(unique(values)) < 2) return(1)
  tryCatch(stats::wilcox.test(values ~ is_enriched_group, exact = FALSE)$p.value, error = function(e) 1)
}

make_biomarker_statistics <- function(long_data, group_order, pseudocount) {
  taxa <- unique(long_data$taxon)
  rows <- vector("list", length(taxa))

  for (i in seq_along(taxa)) {
    taxon_name <- taxa[i]
    taxon_data <- long_data[long_data$taxon == taxon_name, , drop = FALSE]
    taxon_data$group <- factor(taxon_data$group, levels = group_order)

    group_means <- tapply(taxon_data$relative_abundance, taxon_data$group, mean, na.rm = TRUE)
    group_means <- group_means[group_order]
    group_means[is.na(group_means)] <- 0

    enriched_group <- names(group_means)[which.max(group_means)]
    mean_enriched <- unname(group_means[enriched_group])
    mean_other <- mean(group_means[names(group_means) != enriched_group])
    effect_score <- log10((mean_enriched + pseudocount) / (mean_other + pseudocount))

    rows[[i]] <- data.frame(
      taxon = taxon_name,
      enriched_group = enriched_group,
      prevalence = mean(taxon_data$relative_abundance > 0),
      overall_mean_abundance = mean(taxon_data$relative_abundance),
      mean_enriched_group = mean_enriched,
      mean_other_groups = mean_other,
      lefse_like_score = effect_score,
      kruskal_p = safe_kruskal_p(taxon_data$relative_abundance, taxon_data$group),
      wilcoxon_one_vs_rest_p = safe_wilcox_p(
        taxon_data$relative_abundance,
        taxon_data$group == enriched_group
      ),
      stringsAsFactors = FALSE
    )
  }

  statistics <- do.call(rbind, rows)
  statistics$kruskal_fdr <- p.adjust(statistics$kruskal_p, method = "BH")
  statistics[order(statistics$kruskal_fdr, -statistics$lefse_like_score), , drop = FALSE]
}

select_biomarkers <- function(statistics, min_prevalence, min_mean_abundance, fdr_cutoff, effect_score_cutoff, max_biomarkers_to_plot, minimum_biomarkers_to_plot) {
  filtered <- statistics[
    statistics$prevalence >= min_prevalence &
      statistics$overall_mean_abundance >= min_mean_abundance &
      statistics$kruskal_fdr <= fdr_cutoff &
      statistics$lefse_like_score >= effect_score_cutoff,
    ,
    drop = FALSE
  ]
  filtered$selection_reason <- "passes_filter"

  if (nrow(filtered) < minimum_biomarkers_to_plot) {
    relaxed <- statistics[
      statistics$prevalence >= min_prevalence &
        statistics$overall_mean_abundance >= min_mean_abundance,
      ,
      drop = FALSE
    ]
    relaxed <- relaxed[order(relaxed$kruskal_fdr, -relaxed$lefse_like_score), , drop = FALSE]
    relaxed <- head(relaxed, max_biomarkers_to_plot)
    relaxed$selection_reason <- "relaxed_for_demo_plot"
    filtered <- relaxed
  } else {
    filtered <- filtered[order(filtered$enriched_group, -filtered$lefse_like_score), , drop = FALSE]
    filtered <- head(filtered, max_biomarkers_to_plot)
  }

  filtered
}
