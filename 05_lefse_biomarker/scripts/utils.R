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

taxon_table_to_long <- function(taxon_relative_abundance) {
  sample_columns <- setdiff(colnames(taxon_relative_abundance), c("taxon", "rank", "taxon_name", "overall_mean_abundance"))
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

select_biomarkers <- function(statistics, min_prevalence, min_mean_abundance, fdr_cutoff, effect_score_cutoff, max_biomarkers_to_plot, minimum_biomarkers_to_plot, group_order) {
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
    filtered <- filtered[order(filtered$kruskal_fdr, -filtered$lefse_like_score), , drop = FALSE]
  }

  per_group_target <- max(1, floor(max_biomarkers_to_plot / length(group_order)))
  balanced <- do.call(
    rbind,
    lapply(group_order, function(group_name) {
      group_rows <- filtered[filtered$enriched_group == group_name, , drop = FALSE]
      group_rows <- group_rows[order(group_rows$kruskal_fdr, -group_rows$lefse_like_score), , drop = FALSE]
      head(group_rows, per_group_target)
    })
  )

  remaining_slots <- max_biomarkers_to_plot - nrow(balanced)
  if (remaining_slots > 0) {
    remaining <- filtered[!filtered$taxon %in% balanced$taxon, , drop = FALSE]
    remaining <- remaining[order(remaining$kruskal_fdr, -remaining$lefse_like_score), , drop = FALSE]
    balanced <- rbind(balanced, head(remaining, remaining_slots))
  }

  balanced <- balanced[order(match(balanced$enriched_group, group_order), -balanced$lefse_like_score), , drop = FALSE]
  balanced
}

p_to_stars <- function(p_value) {
  ifelse(
    p_value <= 0.001, "***",
    ifelse(p_value <= 0.01, "**", ifelse(p_value <= 0.05, "*", ""))
  )
}

summarize_biomarker_abundance <- function(long_data, biomarkers, group_order) {
  selected <- long_data[long_data$taxon %in% biomarkers$taxon, , drop = FALSE]
  selected$group <- factor(selected$group, levels = group_order)

  mean_table <- aggregate(relative_abundance ~ taxon + group, data = selected, FUN = mean)
  se_table <- aggregate(relative_abundance ~ taxon + group, data = selected, FUN = standard_error)
  colnames(mean_table)[3] <- "mean_relative_abundance"
  colnames(se_table)[3] <- "se_relative_abundance"

  summary <- merge(mean_table, se_table, by = c("taxon", "group"), sort = FALSE)
  summary <- merge(
    summary,
    biomarkers[, c("taxon", "enriched_group", "kruskal_fdr", "lefse_like_score")],
    by = "taxon",
    all.x = TRUE,
    sort = FALSE
  )
  summary$mean_percent <- summary$mean_relative_abundance * 100
  summary$se_percent <- summary$se_relative_abundance * 100
  summary$significance <- p_to_stars(summary$kruskal_fdr)
  summary
}

build_cladogram_tables <- function(taxonomy, all_rank_abundance, biomarkers, taxonomic_levels, max_labels) {
  lineage_columns <- c("Kingdom", taxonomic_levels)
  lineage <- unique(taxonomy[, lineage_columns, drop = FALSE])
  lineage <- lineage[do.call(order, lineage[, lineage_columns, drop = FALSE]), , drop = FALSE]

  genus_order <- unique(lineage$Genus)
  genus_angles <- seq(pi / 2, pi / 2 - 2 * pi, length.out = length(genus_order) + 1)[seq_along(genus_order)]
  names(genus_angles) <- genus_order

  nodes <- data.frame(
    node_key = "Kingdom|Bacteria",
    taxon = "k__Bacteria",
    rank = "Kingdom",
    taxon_name = "Bacteria",
    parent_key = NA_character_,
    stringsAsFactors = FALSE
  )

  for (rank in taxonomic_levels) {
    parent_rank <- if (rank == taxonomic_levels[1]) "Kingdom" else taxonomic_levels[match(rank, taxonomic_levels) - 1]
    rank_nodes <- unique(lineage[, c(parent_rank, rank), drop = FALSE])
    colnames(rank_nodes) <- c("parent_name", "taxon_name")
    rank_nodes$node_key <- paste(rank, rank_nodes$taxon_name, sep = "|")
    rank_nodes$taxon <- paste0(substr(rank, 1, 1), "__", rank_nodes$taxon_name)
    rank_nodes$rank <- rank
    rank_nodes$parent_key <- paste(parent_rank, rank_nodes$parent_name, sep = "|")
    nodes <- rbind(nodes, rank_nodes[, c("node_key", "taxon", "rank", "taxon_name", "parent_key")])
  }

  node_angle <- numeric(nrow(nodes))
  for (i in seq_len(nrow(nodes))) {
    node <- nodes[i, ]
    if (node$rank == "Kingdom") {
      descendant_genera <- genus_order
    } else {
      subset_lineage <- lineage[lineage[[node$rank]] == node$taxon_name, , drop = FALSE]
      descendant_genera <- unique(subset_lineage$Genus)
    }
    node_angle[i] <- mean(genus_angles[descendant_genera], na.rm = TRUE)
  }

  depth_lookup <- c(Kingdom = 0, setNames(seq_along(taxonomic_levels), taxonomic_levels))
  nodes$depth <- unname(depth_lookup[nodes$rank])
  nodes$radius <- nodes$depth
  nodes$angle <- node_angle
  nodes$x <- nodes$depth * cos(nodes$angle)
  nodes$y <- nodes$depth * sin(nodes$angle)

  abundance_lookup <- all_rank_abundance[, c("taxon", "overall_mean_abundance"), drop = FALSE]
  nodes <- merge(nodes, abundance_lookup, by = "taxon", all.x = TRUE, sort = FALSE)
  nodes$overall_mean_abundance[is.na(nodes$overall_mean_abundance)] <- 0
  nodes$overall_mean_percent <- nodes$overall_mean_abundance * 100

  biomarker_lookup <- biomarkers[, c("taxon", "enriched_group", "lefse_like_score"), drop = FALSE]
  nodes <- merge(nodes, biomarker_lookup, by = "taxon", all.x = TRUE, sort = FALSE)
  nodes$enriched_group <- as.character(nodes$enriched_group)
  nodes$plot_group <- ifelse(is.na(nodes$enriched_group), "Not significant", nodes$enriched_group)

  parent_coordinates <- nodes[, c("node_key", "x", "y"), drop = FALSE]
  colnames(parent_coordinates) <- c("parent_key", "x_parent", "y_parent")
  edges <- merge(nodes[!is.na(nodes$parent_key), c("node_key", "parent_key", "x", "y"), drop = FALSE], parent_coordinates, by = "parent_key", all.x = TRUE, sort = FALSE)

  label_nodes <- nodes[!is.na(nodes$enriched_group), , drop = FALSE]
  label_nodes <- label_nodes[order(-label_nodes$lefse_like_score), , drop = FALSE]
  label_nodes <- head(label_nodes, max_labels)
  label_nodes$label <- c(letters, LETTERS)[seq_len(nrow(label_nodes))]
  label_nodes$label_x <- label_nodes$x * 1.03
  label_nodes$label_y <- label_nodes$y * 1.03

  ring_angles <- seq(0, 2 * pi, length.out = 361)
  rings <- do.call(
    rbind,
    lapply(seq_along(taxonomic_levels), function(radius_value) {
      data.frame(
        rank = taxonomic_levels[radius_value],
        radius = radius_value,
        angle = ring_angles,
        x = radius_value * cos(ring_angles),
        y = radius_value * sin(ring_angles),
        stringsAsFactors = FALSE
      )
    })
  )

  rank_labels <- data.frame(
    rank = c("Kingdom", taxonomic_levels),
    radius = c(0, seq_along(taxonomic_levels)),
    x = c(0, seq_along(taxonomic_levels)),
    y = 0,
    label = c("Kingdom", taxonomic_levels),
    stringsAsFactors = FALSE
  )

  list(nodes = nodes, edges = edges, labels = label_nodes, rings = rings, rank_labels = rank_labels)
}

save_two_panel_plot <- function(filename, left_plot, right_plot, width, height, left_width = 0.42, dpi = 300) {
  extension <- tolower(tools::file_ext(filename))
  if (extension == "pdf") {
    grDevices::pdf(filename, width = width, height = height)
  } else if (extension == "png") {
    grDevices::png(filename, width = width, height = height, units = "in", res = dpi)
  } else {
    stop("Unsupported output extension for combined plot: ", extension, call. = FALSE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  grid::grid.newpage()
  layout <- grid::grid.layout(
    nrow = 1,
    ncol = 2,
    widths = grid::unit(c(left_width, 1 - left_width), "null")
  )
  grid::pushViewport(grid::viewport(layout = layout))
  print(left_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(right_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  invisible(TRUE)
}
