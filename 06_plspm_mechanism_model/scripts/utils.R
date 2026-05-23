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

check_sample_ids <- function(sample_metadata, environmental_variables, abundance_table, functional_annotation) {
  metadata_ids <- sample_metadata$sample_id
  id_tables <- list(
    environmental_variables = environmental_variables$sample_id,
    abundance_table = setdiff(colnames(abundance_table), "feature_id"),
    functional_annotation_table = unique(functional_annotation$sample_id)
  )
  problems <- c()
  for (table_name in names(id_tables)) {
    ids <- id_tables[[table_name]]
    missing_from_metadata <- setdiff(ids, metadata_ids)
    missing_from_table <- setdiff(metadata_ids, ids)
    if (length(missing_from_metadata) > 0) {
      problems <- c(problems, paste0(table_name, " sample_id values missing from metadata: ", paste(missing_from_metadata, collapse = ", ")))
    }
    if (length(missing_from_table) > 0) {
      problems <- c(problems, paste0("metadata sample_id values missing from ", table_name, ": ", paste(missing_from_table, collapse = ", ")))
    }
  }

  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

calculate_alpha_diversity <- function(abundance_table) {
  sample_columns <- setdiff(colnames(abundance_table), "feature_id")
  count_matrix <- as.matrix(abundance_table[, sample_columns, drop = FALSE])
  storage.mode(count_matrix) <- "numeric"

  rows <- vector("list", length(sample_columns))
  for (i in seq_along(sample_columns)) {
    sample_id <- sample_columns[i]
    counts <- count_matrix[, sample_id]
    total <- sum(counts, na.rm = TRUE)
    if (total <= 0) stop("All samples must have positive abundance totals.", call. = FALSE)

    proportion_all <- counts / total
    pseudo_rarefied_counts <- round(proportion_all * 2000)
    present <- pseudo_rarefied_counts[pseudo_rarefied_counts > 0]
    richness <- length(present)
    f1 <- sum(pseudo_rarefied_counts == 1, na.rm = TRUE)
    f2 <- sum(pseudo_rarefied_counts == 2, na.rm = TRUE)
    chao1 <- if (f2 > 0) richness + (f1 * f1) / (2 * f2) else richness + (f1 * (f1 - 1)) / 2
    proportion <- counts[counts > 0] / total
    shannon <- -sum(proportion * log(proportion))
    simpson <- 1 - sum(proportion^2)
    pielou <- if (richness > 1) shannon / log(richness) else 0

    rows[[i]] <- data.frame(
      sample_id = sample_id,
      Chao1 = chao1,
      Shannon = shannon,
      Simpson = simpson,
      Pielou = pielou,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

summarize_function_indicators <- function(functional_annotation) {
  category_counts <- aggregate(
    count ~ sample_id + function_category,
    data = functional_annotation,
    FUN = sum
  )
  total_counts <- aggregate(count ~ sample_id, data = category_counts, FUN = sum)
  colnames(total_counts)[2] <- "total_function_count"
  category_counts <- merge(category_counts, total_counts, by = "sample_id", sort = FALSE)
  category_counts$relative_abundance <- category_counts$count / category_counts$total_function_count

  wanted <- c(
    MetalRes = "Metal resistance",
    Sulfur = "Sulfur cycling",
    Nitrogen = "Nitrogen cycling",
    Carbon = "Carbon cycling"
  )

  sample_ids <- unique(functional_annotation$sample_id)
  output <- data.frame(sample_id = sample_ids, stringsAsFactors = FALSE)
  for (indicator in names(wanted)) {
    category <- wanted[[indicator]]
    subset_data <- category_counts[category_counts$function_category == category, c("sample_id", "relative_abundance"), drop = FALSE]
    colnames(subset_data)[2] <- indicator
    output <- merge(output, subset_data, by = "sample_id", all.x = TRUE, sort = FALSE)
  }
  output[is.na(output)] <- 0
  output
}

build_plspm_input <- function(sample_metadata, environmental_variables, abundance_table, functional_annotation) {
  alpha_diversity <- calculate_alpha_diversity(abundance_table)
  function_indicators <- summarize_function_indicators(functional_annotation)

  input <- merge(sample_metadata[, c("sample_id", "group"), drop = FALSE], environmental_variables, by = "sample_id", sort = FALSE)
  input <- merge(input, function_indicators, by = "sample_id", sort = FALSE)
  input <- merge(input, alpha_diversity, by = "sample_id", sort = FALSE)
  input
}

scale_indicator_table <- function(model_input, indicator_columns) {
  numeric_table <- model_input[, indicator_columns, drop = FALSE]
  for (column in indicator_columns) {
    if (!is.numeric(numeric_table[[column]])) {
      stop("PLS-PM indicator is not numeric: ", column, call. = FALSE)
    }
    if (stats::var(numeric_table[[column]], na.rm = TRUE) == 0) {
      stop("PLS-PM indicator has zero variance: ", column, call. = FALSE)
    }
  }
  scaled <- as.data.frame(scale(numeric_table), stringsAsFactors = FALSE)
  if (any(!is.finite(as.matrix(scaled)))) {
    stop("Scaled PLS-PM input contains non-finite values.", call. = FALSE)
  }
  scaled
}

make_block_table <- function(blocks) {
  do.call(
    rbind,
    lapply(names(blocks), function(block_name) {
      data.frame(
        latent_variable = block_name,
        manifest_variable = blocks[[block_name]],
        stringsAsFactors = FALSE
      )
    })
  )
}

flatten_inner_model <- function(inner_model) {
  rows <- vector("list", length(inner_model))
  for (i in seq_along(inner_model)) {
    target <- names(inner_model)[i]
    table <- as.data.frame(inner_model[[i]])
    table$target_latent_variable <- target
    table$predictor_latent_variable <- rownames(table)
    rownames(table) <- NULL
    rows[[i]] <- table[, c("target_latent_variable", "predictor_latent_variable", setdiff(colnames(table), c("target_latent_variable", "predictor_latent_variable"))), drop = FALSE]
  }
  do.call(rbind, rows)
}

save_plspm_base_plot <- function(pdf_file, png_file, width, height, dpi, plot_function) {
  grDevices::pdf(pdf_file, width = width, height = height, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(plot_function())
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)

  grDevices::png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(plot_function())
  grDevices::dev.off()
  on.exit(NULL, add = FALSE)
  invisible(TRUE)
}

save_effects_plot <- function(effects_table, pdf_file, png_file, width, height, dpi, positive_color, negative_color) {
  plot_function <- function() {
    ordered <- effects_table[order(effects_table$total), , drop = FALSE]
    values <- ordered$total
    labels <- ordered$relationships
    colors <- ifelse(values >= 0, positive_color, negative_color)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mar = c(5, 10, 3, 2))
    graphics::barplot(
      values,
      names.arg = labels,
      horiz = TRUE,
      las = 1,
      col = colors,
      border = "grey40",
      xlab = "Total effect",
      main = "PLS-PM total effects"
    )
    graphics::abline(v = 0, col = "grey35", lwd = 1)
  }

  save_plspm_base_plot(pdf_file, png_file, width, height, dpi, plot_function)
}
