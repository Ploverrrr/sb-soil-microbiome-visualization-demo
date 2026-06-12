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

check_sample_ids <- function(sample_metadata, functional_annotation) {
  missing_in_metadata <- setdiff(unique(functional_annotation$sample_id), sample_metadata$sample_id)
  missing_in_function <- setdiff(sample_metadata$sample_id, unique(functional_annotation$sample_id))

  problems <- c()
  if (length(missing_in_metadata) > 0) {
    problems <- c(problems, paste0("Functional annotation sample_id values missing from metadata: ", paste(missing_in_metadata, collapse = ", ")))
  }
  if (length(missing_in_function) > 0) {
    problems <- c(problems, paste0("Metadata sample_id values missing from functional annotation table: ", paste(missing_in_function, collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

aggregate_ko_counts <- function(functional_annotation, ko_id_pattern) {
  filtered <- functional_annotation[grepl(ko_id_pattern, functional_annotation$ko_id), , drop = FALSE]
  if (nrow(filtered) == 0) {
    stop("No KO IDs matched ko_id_pattern: ", ko_id_pattern, call. = FALSE)
  }

  aggregated <- aggregate(count ~ ko_id + sample_id, data = filtered, FUN = sum)
  wide <- reshape(
    aggregated,
    idvar = "ko_id",
    timevar = "sample_id",
    direction = "wide"
  )
  colnames(wide) <- sub("^count[.]", "", colnames(wide))
  wide[is.na(wide)] <- 0
  wide[, setdiff(colnames(wide), "ko_id")] <- round(wide[, setdiff(colnames(wide), "ko_id"), drop = FALSE], 0)
  wide
}

run_deseq2_ko_contrast <- function(ko_count_table, sample_metadata, control_group, treatment_group, group_column, fit_type) {
  sample_subset <- sample_metadata[sample_metadata[[group_column]] %in% c(control_group, treatment_group), , drop = FALSE]
  sample_subset <- sample_subset[order(match(sample_subset[[group_column]], c(control_group, treatment_group)), sample_subset$sample_id), , drop = FALSE]
  count_matrix <- as.matrix(ko_count_table[, sample_subset$sample_id, drop = FALSE])
  storage.mode(count_matrix) <- "integer"
  rownames(count_matrix) <- ko_count_table$ko_id

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
  result$ko_id <- rownames(result)
  result <- result[!is.na(result$pvalue), , drop = FALSE]
  result
}

classify_differential_ko <- function(diff_result, log2fc_threshold, pvalue_threshold) {
  diff_result$change <- ifelse(
    diff_result$pvalue < pvalue_threshold & abs(diff_result$log2FoldChange) >= log2fc_threshold,
    ifelse(diff_result$log2FoldChange > 0, "Up", "Down"),
    "Stable"
  )
  diff_result[order(diff_result$pvalue, -abs(diff_result$log2FoldChange)), , drop = FALSE]
}

select_foreground_kos <- function(diff_result, min_ko_count, log2fc_threshold, pvalue_threshold) {
  selected <- diff_result[
    diff_result$pvalue < pvalue_threshold &
      abs(diff_result$log2FoldChange) >= log2fc_threshold,
    ,
    drop = FALSE
  ]
  if (nrow(selected) < min_ko_count) {
    selected <- diff_result[order(diff_result$pvalue, -abs(diff_result$log2FoldChange)), , drop = FALSE]
    selected <- head(selected, min_ko_count)
    selected$selection_reason <- "top_ranked_for_toy_demo"
  } else {
    selected$selection_reason <- "passes_threshold"
  }
  selected
}

make_toy_kegg_term2gene <- function() {
  data.frame(
    term = c(
      rep("ko01100 Metabolic pathways", 6),
      rep("ko01200 Carbon metabolism", 3),
      rep("ko00190 Oxidative phosphorylation", 1),
      rep("ko02010 ABC transporters", 1),
      rep("ko02020 Two-component system", 1),
      rep("ko00920 Sulfur metabolism", 3),
      rep("ko00910 Nitrogen metabolism", 2),
      rep("ko01120 Microbial metabolism in diverse environments", 7),
      rep("ko02025 Biofilm formation - stress response", 2)
    ),
    gene = c(
      "K02111", "K01190", "K11180", "K17218", "K00958", "K00370",
      "K01190", "K11180", "K00958",
      "K02111",
      "K02003",
      "K07636",
      "K11180", "K17218", "K00958",
      "K10944", "K00370",
      "K02111", "K02003", "K07636", "K11180", "K17218", "K00958", "K01190",
      "K02003", "K07636"
    ),
    stringsAsFactors = FALSE
  )
}

make_toy_module_term2gene <- function() {
  data.frame(
    term = c(
      rep("M00175 Nitrogen metabolism module", 2),
      rep("M00596 Sulfur oxidation module", 2),
      rep("M00176 Assimilatory sulfate reduction module", 2),
      rep("M00157 Energy metabolism module", 2),
      rep("M00690 Stress adaptation module", 2),
      rep("M00001 Central carbon module", 2)
    ),
    gene = c(
      "K10944", "K00370",
      "K17218", "K00958",
      "K11180", "K00958",
      "K02111", "K02003",
      "K02003", "K07636",
      "K01190", "K02111"
    ),
    stringsAsFactors = FALSE
  )
}

run_offline_enrichment <- function(ko_list, universe, term2gene, pvalue_cutoff, p_adjust_method, qvalue_cutoff, min_gs_size, max_gs_size) {
  clusterProfiler::enricher(
    gene = ko_list,
    universe = universe,
    TERM2GENE = term2gene,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    qvalueCutoff = qvalue_cutoff,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size
  )
}

run_clusterprofiler_kegg <- function(ko_list, pvalue_cutoff, p_adjust_method, qvalue_cutoff, min_gs_size, max_gs_size) {
  ko_result <- clusterProfiler::enrichKEGG(
    gene = ko_list,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    qvalueCutoff = qvalue_cutoff
  )

  module_result <- clusterProfiler::enrichMKEGG(
    gene = ko_list,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    qvalueCutoff = qvalue_cutoff
  )

  list(ko = ko_result, module = module_result)
}

run_enrichment_backend <- function(enrichment_backend, ko_list, universe, pvalue_cutoff, p_adjust_method, qvalue_cutoff, min_gs_size, max_gs_size) {
  if (identical(enrichment_backend, "clusterprofiler_kegg")) {
    message("Running online clusterProfiler KEGG backend: enrichKEGG() + enrichMKEGG().")
    return(run_clusterprofiler_kegg(
      ko_list = ko_list,
      pvalue_cutoff = pvalue_cutoff,
      p_adjust_method = p_adjust_method,
      qvalue_cutoff = qvalue_cutoff,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size
    ))
  }

  if (identical(enrichment_backend, "toy_offline")) {
    message("Running reproducible offline toy backend: clusterProfiler::enricher() with KEGG-like TERM2GENE.")
    return(list(
      ko = run_offline_enrichment(
        ko_list = ko_list,
        universe = universe,
        term2gene = make_toy_kegg_term2gene(),
        pvalue_cutoff = pvalue_cutoff,
        p_adjust_method = p_adjust_method,
        qvalue_cutoff = qvalue_cutoff,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size
      ),
      module = run_offline_enrichment(
        ko_list = ko_list,
        universe = universe,
        term2gene = make_toy_module_term2gene(),
        pvalue_cutoff = pvalue_cutoff,
        p_adjust_method = p_adjust_method,
        qvalue_cutoff = qvalue_cutoff,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size
      )
    ))
  }

  message("Trying online clusterProfiler KEGG backend first; falling back to offline toy backend only if KEGG access fails.")
  tryCatch(
    run_enrichment_backend(
      "clusterprofiler_kegg",
      ko_list = ko_list,
      universe = universe,
      pvalue_cutoff = pvalue_cutoff,
      p_adjust_method = p_adjust_method,
      qvalue_cutoff = qvalue_cutoff,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size
    ),
    error = function(error) {
      message("Online KEGG backend failed: ", conditionMessage(error))
      message("Falling back to the offline toy backend so the demo remains reproducible.")
      run_enrichment_backend(
        "toy_offline",
        ko_list = ko_list,
        universe = universe,
        pvalue_cutoff = pvalue_cutoff,
        p_adjust_method = p_adjust_method,
        qvalue_cutoff = qvalue_cutoff,
        min_gs_size = min_gs_size,
        max_gs_size = max_gs_size
      )
    }
  )
}

make_enrichment_plots <- function(enrichment_result, show_category, color_by) {
  list(
    bar = barplot(enrichment_result, showCategory = show_category, color = color_by),
    dot = dotplot(enrichment_result, showCategory = show_category, color = color_by)
  )
}

parse_gene_ratio <- function(x) {
  vapply(strsplit(x, "/", fixed = TRUE), function(parts) as.numeric(parts[1]) / as.numeric(parts[2]), numeric(1))
}

make_enrichment_bubble_plot <- function(enrichment_table, title, color_palette, size_range) {
  if (nrow(enrichment_table) == 0) {
    stop("Cannot draw enrichment bubble plot from an empty table.", call. = FALSE)
  }
  enrichment_table$GeneRatioDecimal <- parse_gene_ratio(enrichment_table$GeneRatio)
  enrichment_table$Description <- factor(enrichment_table$Description, levels = rev(enrichment_table$Description))

  ggplot2::ggplot(
    enrichment_table,
    ggplot2::aes(x = GeneRatioDecimal, y = Description, size = Count, color = p.adjust)
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::scale_color_gradient(low = color_palette[1], high = color_palette[2], name = "Adjusted P") +
    ggplot2::scale_size_continuous(range = size_range, name = "Count") +
    ggplot2::labs(title = title, x = "GeneRatio", y = "KEGG term") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

save_ggplot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(TRUE)
}
