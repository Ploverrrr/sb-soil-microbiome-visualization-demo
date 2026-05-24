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

check_sample_ids <- function(environmental_variables, functional_annotation) {
  missing_in_environment <- setdiff(unique(functional_annotation$sample_id), environmental_variables$sample_id)
  missing_in_function <- setdiff(environmental_variables$sample_id, unique(functional_annotation$sample_id))

  problems <- c()
  if (length(missing_in_environment) > 0) {
    problems <- c(problems, paste0("Functional annotation sample_id values missing from environmental table: ", paste(missing_in_environment, collapse = ", ")))
  }
  if (length(missing_in_function) > 0) {
    problems <- c(problems, paste0("Environmental sample_id values missing from functional annotation table: ", paste(missing_in_function, collapse = ", ")))
  }
  if (length(problems) > 0) stop(paste(problems, collapse = "\n"), call. = FALSE)
  invisible(TRUE)
}

make_sulfur_gene_matrix <- function(functional_annotation, target_ko_names, value_column) {
  target_kos <- names(target_ko_names)
  subset_data <- functional_annotation[functional_annotation$ko_id %in% target_kos, , drop = FALSE]
  if (nrow(subset_data) == 0) {
    stop("No target sulfur KO IDs were found in functional_annotation_table.csv.", call. = FALSE)
  }

  aggregated <- aggregate(
    subset_data[[value_column]],
    by = list(ko_id = subset_data$ko_id, sample_id = subset_data$sample_id),
    FUN = sum
  )
  colnames(aggregated)[3] <- value_column

  wide <- reshape(
    aggregated,
    idvar = "sample_id",
    timevar = "ko_id",
    direction = "wide"
  )
  colnames(wide) <- sub(paste0("^", value_column, "[.]"), "", colnames(wide))

  for (ko in target_kos) {
    if (!ko %in% colnames(wide)) wide[[ko]] <- NA_real_
  }
  gene_matrix <- wide[, c("sample_id", target_kos), drop = FALSE]
  colnames(gene_matrix) <- c("sample_id", unname(target_ko_names[target_kos]))

  if ("dsrA" %in% colnames(gene_matrix) || "dsrB" %in% colnames(gene_matrix)) {
    dsr_cols <- intersect(c("dsrA", "dsrB"), colnames(gene_matrix))
    gene_matrix$dsrAB <- rowMeans(gene_matrix[, dsr_cols, drop = FALSE], na.rm = TRUE)
  }
  sox_cols <- intersect(c("soxB", "soxC", "soxD"), colnames(gene_matrix))
  if (length(sox_cols) > 0) {
    gene_matrix$soxBCD <- rowMeans(gene_matrix[, sox_cols, drop = FALSE], na.rm = TRUE)
  }
  if ("sat" %in% colnames(gene_matrix)) {
    gene_matrix$sat <- gene_matrix$sat
  }

  gene_matrix
}

make_pearson_table <- function(data, gene_vars, contaminant_vars) {
  rows <- list()
  index <- 1
  for (gene in gene_vars) {
    for (contaminant in contaminant_vars) {
      complete <- stats::complete.cases(data[, c(gene, contaminant), drop = FALSE])
      test <- stats::cor.test(data[[gene]][complete], data[[contaminant]][complete], method = "pearson")
      rows[[index]] <- data.frame(
        Gene = gene,
        Contaminant = contaminant,
        r = unname(test$estimate),
        p = test$p.value,
        sig = p_to_stars(test$p.value),
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

cor_pmat <- function(data, method = "pearson") {
  p_matrix <- matrix(NA_real_, ncol(data), ncol(data))
  colnames(p_matrix) <- colnames(data)
  rownames(p_matrix) <- colnames(data)
  for (i in seq_len(ncol(data))) {
    for (j in seq_len(ncol(data))) {
      complete <- stats::complete.cases(data[, c(i, j), drop = FALSE])
      p_matrix[i, j] <- stats::cor.test(data[[i]][complete], data[[j]][complete], method = method)$p.value
    }
  }
  p_matrix
}

select_scatter_pairs <- function(pearson_table, max_pairs) {
  significant <- pearson_table[pearson_table$p < 0.05, , drop = FALSE]
  significant <- significant[order(significant$p, -abs(significant$r)), , drop = FALSE]
  if (nrow(significant) == 0) {
    significant <- pearson_table[order(pearson_table$p, -abs(pearson_table$r)), , drop = FALSE]
  }

  diverse <- do.call(
    rbind,
    lapply(unique(significant$Gene), function(gene) {
      head(significant[significant$Gene == gene, , drop = FALSE], 1)
    })
  )
  diverse <- diverse[order(diverse$p, -abs(diverse$r)), , drop = FALSE]
  selected <- head(diverse, max_pairs)

  if (nrow(selected) < max_pairs) {
    remaining <- significant[
      !paste(significant$Gene, significant$Contaminant) %in% paste(selected$Gene, selected$Contaminant),
      ,
      drop = FALSE
    ]
    selected <- rbind(selected, head(remaining, max_pairs - nrow(selected)))
  }
  selected
}

save_ggplot_pair <- function(plot, pdf_file, png_file, width, height, dpi) {
  ggplot2::ggsave(pdf_file, plot, width = width, height = height, units = "in", bg = "white")
  ggplot2::ggsave(png_file, plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(TRUE)
}

fit_mlr_models <- function(data, predictor_vars, response_vars) {
  rows <- list()
  plots <- list()
  index <- 1

  for (response in response_vars) {
    model_data <- data[, c(predictor_vars, response), drop = FALSE]
    model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]
    formula <- stats::as.formula(paste(response, "~", paste(predictor_vars, collapse = " + ")))
    fit <- stats::lm(formula, data = model_data)
    fit_summary <- summary(fit)
    f_p <- stats::pf(
      fit_summary$fstatistic[1],
      fit_summary$fstatistic[2],
      fit_summary$fstatistic[3],
      lower.tail = FALSE
    )

    coef_table <- as.data.frame(fit_summary$coefficients)
    coef_table$Term <- rownames(coef_table)
    coef_table$Response <- response
    coef_table$R2 <- fit_summary$r.squared
    coef_table$R2_adj <- fit_summary$adj.r.squared
    coef_table$model_p <- f_p
    rows[[index]] <- coef_table
    index <- index + 1

    fitted_values <- stats::fitted(fit)
    actual_values <- model_data[[response]]
    plots[[paste0(response, "_fit")]] <- data.frame(
      response = response,
      fitted = fitted_values,
      actual = actual_values,
      R2 = fit_summary$r.squared,
      R2_adj = fit_summary$adj.r.squared,
      model_p = f_p,
      stringsAsFactors = FALSE
    )

    scaled_data <- as.data.frame(scale(model_data))
    scaled_fit <- stats::lm(formula, data = scaled_data)
    ci <- stats::confint(scaled_fit)
    beta_terms <- predictor_vars
    plots[[paste0(response, "_beta")]] <- data.frame(
      response = response,
      Gene = beta_terms,
      beta = stats::coef(scaled_fit)[beta_terms],
      lower = ci[beta_terms, 1],
      upper = ci[beta_terms, 2],
      p = coef_table[match(beta_terms, coef_table$Term), "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
  }

  mlr_table <- do.call(rbind, rows)
  rownames(mlr_table) <- NULL
  colnames(mlr_table)[1:4] <- c("Estimate", "Std.Error", "t.value", "p.value")
  list(table = mlr_table, plot_data = plots)
}
