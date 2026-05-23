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

check_sample_ids <- function(sample_metadata, abundance_table) {
  sample_columns <- setdiff(colnames(abundance_table), "feature_id")
  missing_in_metadata <- setdiff(sample_columns, sample_metadata$sample_id)
  missing_in_abundance <- setdiff(sample_metadata$sample_id, sample_columns)

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

check_feature_ids <- function(abundance_table, taxonomy_table) {
  missing_in_taxonomy <- setdiff(abundance_table$feature_id, taxonomy_table$feature_id)
  missing_in_abundance <- setdiff(taxonomy_table$feature_id, abundance_table$feature_id)

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

aggregate_counts_to_taxon <- function(abundance_table, taxonomy_table, target_taxonomic_level) {
  sample_columns <- setdiff(colnames(abundance_table), "feature_id")
  merged <- merge(
    taxonomy_table[, c("feature_id", target_taxonomic_level), drop = FALSE],
    abundance_table,
    by = "feature_id",
    all.y = TRUE,
    sort = FALSE
  )
  merged$taxon <- merged[[target_taxonomic_level]]
  merged$taxon[is.na(merged$taxon) | merged$taxon == ""] <- "Unclassified"

  aggregated <- aggregate(
    merged[, sample_columns, drop = FALSE],
    by = list(taxon = merged$taxon),
    FUN = sum
  )
  aggregated[, sample_columns] <- round(aggregated[, sample_columns, drop = FALSE], 0)
  aggregated
}

run_deseq2_contrast <- function(count_table, sample_metadata, control_group, treatment_group, group_column, fit_type) {
  sample_subset <- sample_metadata[sample_metadata[[group_column]] %in% c(control_group, treatment_group), , drop = FALSE]
  sample_subset <- sample_subset[order(match(sample_subset[[group_column]], c(control_group, treatment_group)), sample_subset$sample_id), , drop = FALSE]
  count_matrix <- as.matrix(count_table[, sample_subset$sample_id, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- count_table$taxon

  col_data <- data.frame(
    condition = factor(
      ifelse(sample_subset[[group_column]] == control_group, "control", "Contaminant"),
      levels = c("control", "Contaminant")
    ),
    row.names = sample_subset$sample_id
  )

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~condition
  )
  dds <- DESeq2::DESeq(dds, fitType = fit_type, quiet = TRUE)
  result <- DESeq2::results(dds, contrast = c("condition", "Contaminant", "control"))
  result <- data.frame(result, stringsAsFactors = FALSE, check.names = FALSE)
  result$taxon <- rownames(result)
  result <- result[!is.na(result$pvalue), , drop = FALSE]
  result
}

classify_change <- function(diff_result, log2fc_threshold, pvalue_threshold) {
  diff_result$change <- ifelse(
    diff_result$pvalue < pvalue_threshold & abs(diff_result$log2FoldChange) >= log2fc_threshold,
    ifelse(diff_result$log2FoldChange > 0, "Up", "Down"),
    "Stable"
  )
  diff_result
}

select_top_labeled_taxa <- function(diff_result, top_label_n) {
  up <- diff_result[diff_result$change == "Up", , drop = FALSE]
  up <- up[order(up$log2FoldChange, decreasing = TRUE), , drop = FALSE]
  down <- diff_result[diff_result$change == "Down", , drop = FALSE]
  down <- down[order(down$log2FoldChange), , drop = FALSE]
  rbind(head(down, top_label_n), head(up, top_label_n))
}

select_heatmap_taxa <- function(diff_result, top_n_each_direction) {
  changed <- diff_result[diff_result$change != "Stable", , drop = FALSE]
  if (nrow(changed) == 0) {
    changed <- diff_result[order(diff_result$pvalue), , drop = FALSE]
    return(head(changed$taxon, top_n_each_direction * 2))
  }

  up <- changed[changed$change == "Up", , drop = FALSE]
  up <- up[order(up$log2FoldChange, decreasing = TRUE), , drop = FALSE]
  down <- changed[changed$change == "Down", , drop = FALSE]
  down <- down[order(down$log2FoldChange), , drop = FALSE]
  unique(c(head(up$taxon, top_n_each_direction), head(down$taxon, top_n_each_direction)))
}

make_zscore_matrix <- function(count_table, selected_taxa, sample_order) {
  selected_counts <- count_table[count_table$taxon %in% selected_taxa, c("taxon", sample_order), drop = FALSE]
  selected_counts <- selected_counts[match(selected_taxa, selected_counts$taxon), , drop = FALSE]
  matrix_data <- as.matrix(selected_counts[, sample_order, drop = FALSE])
  storage.mode(matrix_data) <- "numeric"
  rownames(matrix_data) <- selected_counts$taxon
  log_matrix <- log2(matrix_data + 1)
  scaled <- t(scale(t(log_matrix)))
  scaled[is.na(scaled)] <- 0
  scaled
}

save_ggplot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(TRUE)
}

save_complex_heatmap_pair <- function(heatmap_object, pdf_file, png_file, width, height, dpi) {
  grDevices::pdf(pdf_file, width = width, height = height, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(heatmap_object)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)

  grDevices::png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  ComplexHeatmap::draw(heatmap_object)
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  invisible(TRUE)
}

save_circular_heatmap_pair <- function(matrix_data, pdf_file, png_file, width, height, dpi, color_function, row_name_color) {
  plot_function <- function() {
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(gap.after = c(30))
    circlize::circos.heatmap(
      matrix_data,
      col = color_function,
      dend.side = "inside",
      rownames.side = "outside",
      rownames.col = row_name_color,
      rownames.cex = 0.5,
      rownames.font = 0.5,
      bg.border = "black",
      show.sector.labels = FALSE,
      track.height = 0.4,
      dend.track.height = 0.2,
      cluster = TRUE
    )
    ComplexHeatmap::draw(ComplexHeatmap::Legend(
      title = "Expression",
      col_fun = color_function,
      direction = "vertical",
      title_position = "topcenter"
    ))
  }

  grDevices::pdf(pdf_file, width = width, height = height, bg = "white")
  plot_function()
  grDevices::dev.off()

  grDevices::png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white")
  plot_function()
  grDevices::dev.off()
  invisible(TRUE)
}
