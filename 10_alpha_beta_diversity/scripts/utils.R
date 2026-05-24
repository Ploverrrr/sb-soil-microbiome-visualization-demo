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

  feature_ids <- abundance$feature_id
  abundance_matrix <- as.matrix(abundance[, setdiff(colnames(abundance), "feature_id"), drop = FALSE])
  storage.mode(abundance_matrix) <- "numeric"
  rownames(abundance_matrix) <- feature_ids

  if (anyNA(abundance_matrix)) {
    stop("abundance_table.csv contains non-numeric or missing abundance values.", call. = FALSE)
  }
  if (any(abundance_matrix < 0)) {
    stop("abundance_table.csv contains negative abundance values.", call. = FALSE)
  }

  abundance_matrix
}

align_community_and_metadata <- function(abundance_matrix, metadata, group_column) {
  sample_ids <- colnames(abundance_matrix)
  missing_in_metadata <- setdiff(sample_ids, metadata$sample_id)
  missing_in_abundance <- setdiff(metadata$sample_id, sample_ids)

  problems <- c()
  if (length(missing_in_metadata) > 0) {
    problems <- c(problems, paste0("Samples in abundance_table.csv missing from sample_metadata.csv: ", paste(missing_in_metadata, collapse = ", ")))
  }
  if (length(missing_in_abundance) > 0) {
    problems <- c(problems, paste0("Samples in sample_metadata.csv missing from abundance_table.csv: ", paste(missing_in_abundance, collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)

  metadata <- metadata[match(sample_ids, metadata$sample_id), , drop = FALSE]
  check_required_columns(metadata, c("sample_id", group_column), "sample_metadata.csv")

  community <- t(abundance_matrix)
  rownames(community) <- sample_ids
  list(community = community, metadata = metadata)
}

calculate_alpha_indices <- function(community_counts, base = exp(1)) {
  observed_species <- vegan::specnumber(community_counts)
  richness <- rowSums(community_counts > 0)
  estimate <- vegan::estimateR(community_counts)
  chao1 <- as.numeric(estimate["S.chao1", ])
  ace <- as.numeric(estimate["S.ACE", ])
  shannon <- vegan::diversity(community_counts, index = "shannon", base = base)
  simpson <- vegan::diversity(community_counts, index = "simpson")
  pielou <- shannon / log(observed_species, base = base)
  pielou[!is.finite(pielou)] <- NA_real_
  goods_coverage <- 1 - rowSums(community_counts == 1) / rowSums(community_counts)

  data.frame(
    sample_id = rownames(community_counts),
    observed_species = as.numeric(observed_species),
    richness = as.numeric(richness),
    chao1 = chao1,
    ace = ace,
    shannon = as.numeric(shannon),
    simpson = as.numeric(simpson),
    pielou = as.numeric(pielou),
    goods_coverage = as.numeric(goods_coverage),
    stringsAsFactors = FALSE
  )
}

pairwise_test_table <- function(data, value_columns, group_column, group_order, method = "t.test") {
  comparisons <- utils::combn(group_order, 2, simplify = FALSE)
  rows <- list()
  index <- 1

  for (value_column in value_columns) {
    for (comparison in comparisons) {
      group_a <- comparison[1]
      group_b <- comparison[2]
      values_a <- data[data[[group_column]] == group_a, value_column]
      values_b <- data[data[[group_column]] == group_b, value_column]
      p_value <- tryCatch(
        {
          test <- switch(
            method,
            "wilcox.test" = stats::wilcox.test(values_a, values_b),
            stats::t.test(values_a, values_b)
          )
          test$p.value
        },
        error = function(e) NA_real_
      )
      rows[[index]] <- data.frame(
        metric = value_column,
        group_1 = group_a,
        group_2 = group_b,
        method = method,
        p = p_value,
        p_signif = ifelse(is.na(p_value), "not_tested", p_to_stars(p_value)),
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }

  do.call(rbind, rows)
}

p_to_stars <- function(p_value) {
  ifelse(
    p_value < 0.001, "***",
    ifelse(p_value < 0.01, "**", ifelse(p_value < 0.05, "*", "ns"))
  )
}

make_pairwise_distance_table <- function(distance_matrix, metadata, group_column, group_order) {
  sample_ids <- rownames(distance_matrix)
  rows <- list()
  index <- 1

  for (i in seq_len(length(sample_ids) - 1)) {
    for (j in seq((i + 1), length(sample_ids))) {
      sample_1 <- sample_ids[i]
      sample_2 <- sample_ids[j]
      group_1 <- metadata[[group_column]][match(sample_1, metadata$sample_id)]
      group_2 <- metadata[[group_column]][match(sample_2, metadata$sample_id)]
      ordered_groups <- group_order[group_order %in% c(group_1, group_2)]
      comparison <- if (group_1 == group_2) group_1 else paste(ordered_groups, collapse = " vs ")
      rows[[index]] <- data.frame(
        sample_1 = sample_1,
        sample_2 = sample_2,
        group_1 = group_1,
        group_2 = group_2,
        comparison = comparison,
        bray_curtis = distance_matrix[sample_1, sample_2],
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }

  do.call(rbind, rows)
}

make_focal_distance_table <- function(distance_matrix, metadata, group_column) {
  sample_ids <- rownames(distance_matrix)
  rows <- list()
  index <- 1

  for (sample_id in sample_ids) {
    other_samples <- setdiff(sample_ids, sample_id)
    focal_group <- metadata[[group_column]][match(sample_id, metadata$sample_id)]
    for (other_sample in other_samples) {
      rows[[index]] <- data.frame(
        sample_id = sample_id,
        comparison_sample_id = other_sample,
        group = focal_group,
        bray_curtis = distance_matrix[sample_id, other_sample],
        stringsAsFactors = FALSE
      )
      index <- index + 1
    }
  }

  do.call(rbind, rows)
}

save_ggplot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(TRUE)
}
