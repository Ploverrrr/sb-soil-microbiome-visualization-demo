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

make_microeco_tables <- function(sample_metadata, abundance_table, taxonomy_table, group_column, group_order) {
  sample_columns <- setdiff(colnames(abundance_table), "feature_id")

  otu_table <- abundance_table[, sample_columns, drop = FALSE]
  rownames(otu_table) <- abundance_table$feature_id

  sample_table <- sample_metadata
  rownames(sample_table) <- sample_table$sample_id
  sample_table <- sample_table[sample_columns, , drop = FALSE]
  sample_table$Group <- factor(sample_table[[group_column]], levels = group_order)

  rank_columns <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  prefix <- c(
    Kingdom = "k__",
    Phylum = "p__",
    Class = "c__",
    Order = "o__",
    Family = "f__",
    Genus = "g__"
  )

  tax_table <- taxonomy_table[, rank_columns, drop = FALSE]
  for (rank in rank_columns) {
    tax_table[[rank]] <- paste0(prefix[[rank]], tax_table[[rank]])
  }
  rownames(tax_table) <- taxonomy_table$feature_id
  tax_table <- tax_table[rownames(otu_table), , drop = FALSE]

  list(
    otu_table = otu_table,
    sample_table = sample_table,
    tax_table = tax_table
  )
}

strip_tax_prefixes <- function(x) {
  cleaned <- gsub("(^|[|])([kpcofg]__)", "\\1", x)
  gsub("\\|", " | ", cleaned)
}

save_plot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
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
  aggregated <- aggregate(merged[, sample_columns, drop = FALSE], by = list(taxon = merged$taxon), FUN = sum)
  aggregated$rank <- target_taxonomic_level
  aggregated$taxon_name <- aggregated$taxon
  aggregated$taxon <- paste0(substr(target_taxonomic_level, 1, 1), "__", aggregated$taxon)
  aggregated[, c("taxon", "rank", "taxon_name", sample_columns), drop = FALSE]
}

aggregate_to_all_taxonomic_levels <- function(relative_abundance, taxonomy, taxonomic_levels) {
  sample_columns <- setdiff(colnames(relative_abundance), "feature_id")
  rows <- vector("list", length(taxonomic_levels))

  for (i in seq_along(taxonomic_levels)) {
    rank <- taxonomic_levels[i]
    rank_table <- aggregate_to_taxon(relative_abundance, taxonomy, rank)
    rows[[i]] <- rank_table[, c("taxon", "rank", "taxon_name", sample_columns), drop = FALSE]
  }

  combined <- do.call(rbind, rows)
  combined[order(match(combined$rank, taxonomic_levels), combined$taxon_name), , drop = FALSE]
}

cleanup_default_rplots <- function(path = "Rplots.pdf") {
  if (file.exists(path)) unlink(path)
  invisible(TRUE)
}
