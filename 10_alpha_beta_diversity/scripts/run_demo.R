# Purpose:
#   Build an independent alpha/beta diversity and ordination demo from shared
#   simulated soil microbiome toy data. The workflow follows the original
#   alpha/beta scripts: optional rarefaction, vegan alpha indices, ggpubr
#   boxplots, Bray-Curtis PCoA/PERMANOVA, NMDS, and beta-distance violin/boxplots.
#
# Input files:
#   - ../data/toy_shared/sample_metadata.csv
#   - ../data/toy_shared/abundance_table.csv
#
# Output files:
#   - results/alpha_diversity_indices.csv
#   - results/alpha_diversity_pairwise_tests.csv
#   - results/pcoa_coordinates.csv
#   - results/nmds_coordinates.csv
#   - results/permanova_results.csv
#   - results/beta_pairwise_distances.csv
#   - results/beta_distance_by_sample_group.csv
#   - figures/fig1_alpha_diversity_boxplots.pdf/png
#   - figures/fig2_pcoa_bray_curtis.pdf/png
#   - figures/fig3_nmds_bray_curtis.pdf/png
#   - figures/fig4_beta_distance_boxplot.pdf/png
#
# Required columns:
#   sample_metadata.csv: sample_id, group
#   abundance_table.csv: feature_id plus one numeric column per sample_id
#
# User-editable settings:
#   Edit the settings block below to change input paths, grouping column,
#   rarefaction, alpha metrics, beta distance method, group order, palettes,
#   output file names, figure sizes, and random seed.

# ----------------------------- User settings -----------------------------

shared_data_dir <- "../data/toy_shared"

metadata_file <- file.path(shared_data_dir, "sample_metadata.csv")
abundance_file <- file.path(shared_data_dir, "abundance_table.csv")

group_column <- "group"
group_order <- c("Control", "Tailing", "Mining", "Smelting")
group_palette <- c(
  Control = "#a3cfcd",
  Tailing = "#7d9fc0",
  Mining = "#8482ad",
  Smelting = "#8d658b"
)

random_seed <- 123
rarefy_counts <- TRUE
alpha_min_count_threshold <- 20
alpha_metrics_to_plot <- c("chao1", "shannon")
alpha_test_method <- "t.test"

beta_distance_method <- "bray"
use_relative_abundance_for_beta <- TRUE
nmds_trymax <- 100
ellipse_level_pcoa <- 0.90
ellipse_level_nmds <- 0.95
permutations <- 999

alpha_output <- "alpha_diversity_indices.csv"
alpha_test_output <- "alpha_diversity_pairwise_tests.csv"
pcoa_output <- "pcoa_coordinates.csv"
nmds_output <- "nmds_coordinates.csv"
permanova_output <- "permanova_results.csv"
beta_pairwise_output <- "beta_pairwise_distances.csv"
beta_focal_output <- "beta_distance_by_sample_group.csv"

alpha_pdf <- "fig1_alpha_diversity_boxplots.pdf"
alpha_png <- "fig1_alpha_diversity_boxplots.png"
pcoa_pdf <- "fig2_pcoa_bray_curtis.pdf"
pcoa_png <- "fig2_pcoa_bray_curtis.png"
nmds_pdf <- "fig3_nmds_bray_curtis.pdf"
nmds_png <- "fig3_nmds_bray_curtis.png"
beta_box_pdf <- "fig4_beta_distance_boxplot.pdf"
beta_box_png <- "fig4_beta_distance_boxplot.png"

alpha_width <- 7
alpha_height <- 3.8
ordination_width <- 5.5
ordination_height <- 4.8
beta_box_width <- 5.2
beta_box_height <- 4.3
png_dpi <- 300

# ----------------------------- Main workflow -----------------------------

source(file.path("scripts", "utils.R"))

required_packages <- c("vegan", "ggplot2", "ggpubr", "patchwork", "gghalves", "ggsignif")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this demo, then rerun: Rscript scripts/run_demo.R",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(gghalves)
  library(ggsignif)
})

set.seed(random_seed)

results_dir <- "results"
figures_dir <- "figures"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE)

check_files_exist(c(metadata_file, abundance_file))

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
abundance_matrix <- read_abundance_matrix(abundance_file)

check_required_columns(metadata, c("sample_id", group_column), "sample_metadata.csv")
aligned <- align_community_and_metadata(abundance_matrix, metadata, group_column)
community_counts <- aligned$community
metadata <- aligned$metadata

missing_groups <- setdiff(group_order, unique(metadata[[group_column]]))
if (length(missing_groups) > 0) {
  stop("The following configured groups are missing from sample_metadata.csv: ", paste(missing_groups, collapse = ", "), call. = FALSE)
}
metadata[[group_column]] <- factor(metadata[[group_column]], levels = group_order)

community_counts_for_alpha <- community_counts
if (alpha_min_count_threshold > 0) {
  community_counts_for_alpha[community_counts_for_alpha < alpha_min_count_threshold] <- 0
}

if (any(rowSums(community_counts_for_alpha) == 0)) {
  stop("At least one sample has zero total abundance after alpha_min_count_threshold filtering.", call. = FALSE)
}

if (rarefy_counts) {
  rarefaction_depth <- min(rowSums(community_counts_for_alpha))
  community_for_alpha <- vegan::rrarefy(community_counts_for_alpha, sample = rarefaction_depth)
} else {
  rarefaction_depth <- NA_real_
  community_for_alpha <- community_counts_for_alpha
}

alpha_indices <- calculate_alpha_indices(community_for_alpha)
alpha_indices <- merge(alpha_indices, metadata[, c("sample_id", group_column), drop = FALSE], by = "sample_id", sort = FALSE)
alpha_indices[[group_column]] <- factor(alpha_indices[[group_column]], levels = group_order)
write.csv(alpha_indices, file.path(results_dir, alpha_output), row.names = FALSE)

alpha_tests <- pairwise_test_table(
  data = alpha_indices,
  value_columns = alpha_metrics_to_plot,
  group_column = group_column,
  group_order = group_order,
  method = alpha_test_method
)
write.csv(alpha_tests, file.path(results_dir, alpha_test_output), row.names = FALSE)

alpha_comparisons <- utils::combn(group_order, 2, simplify = FALSE)
alpha_plots <- lapply(alpha_metrics_to_plot, function(metric) {
  testable_comparisons <- alpha_tests[
    alpha_tests$metric == metric & !is.na(alpha_tests$p) & alpha_tests$p < 0.05,
    c("group_1", "group_2"),
    drop = FALSE
  ]
  testable_comparisons <- lapply(seq_len(nrow(testable_comparisons)), function(i) unname(as.character(testable_comparisons[i, ])))

  plot <- ggpubr::ggboxplot(
    alpha_indices,
    x = group_column,
    y = metric,
    fill = group_column,
    palette = unname(group_palette[group_order]),
    width = 0.5,
    add = "jitter",
    add.params = list(size = 1.8, alpha = 0.65)
  )

  if (length(testable_comparisons) > 0) {
    plot <- plot +
      ggpubr::stat_compare_means(
      comparisons = testable_comparisons,
      label = "p.signif",
      method = alpha_test_method,
      size = 3
      )
  }

  plot +
    labs(x = "", y = metric, fill = "") +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 11, face = "bold", color = "black"),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
})
alpha_plot <- patchwork::wrap_plots(alpha_plots, nrow = 1)
save_ggplot_pair(alpha_plot, file.path(figures_dir, alpha_pdf), file.path(figures_dir, alpha_png), alpha_width, alpha_height, png_dpi)

if (use_relative_abundance_for_beta) {
  community_for_beta <- vegan::decostand(community_counts, method = "total")
} else {
  community_for_beta <- community_counts
}

beta_distance <- vegan::vegdist(community_for_beta, method = beta_distance_method, binary = FALSE)
beta_distance_matrix <- as.matrix(beta_distance)

pcoa <- stats::cmdscale(beta_distance, k = 2, eig = TRUE)
positive_eigen_sum <- sum(pcoa$eig[pcoa$eig > 0])
eig_percent <- round(pcoa$eig[1:2] / positive_eigen_sum * 100, 3)
pcoa_coordinates <- data.frame(
  sample_id = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2],
  stringsAsFactors = FALSE
)
pcoa_coordinates <- merge(pcoa_coordinates, metadata[, c("sample_id", group_column), drop = FALSE], by = "sample_id", sort = FALSE)
pcoa_coordinates[[group_column]] <- factor(pcoa_coordinates[[group_column]], levels = group_order)
write.csv(pcoa_coordinates, file.path(results_dir, pcoa_output), row.names = FALSE)

metadata$group_for_permanova <- metadata[[group_column]]
permanova <- vegan::adonis2(beta_distance ~ group_for_permanova, data = metadata, permutations = permutations)
permanova_table <- as.data.frame(permanova)
permanova_table$term <- rownames(permanova_table)
permanova_table <- permanova_table[, c("term", setdiff(colnames(permanova_table), "term")), drop = FALSE]
permanova_table$term <- sub("^group_for_permanova$", group_column, permanova_table$term)
write.csv(permanova_table, file.path(results_dir, permanova_output), row.names = FALSE)

adonis_caption <- paste0(
  "PERMANOVA R2 = ",
  round(permanova_table$R2[permanova_table$term %in% c(group_column, "Model")][1], 3),
  "; P = ",
  signif(permanova_table$`Pr(>F)`[permanova_table$term %in% c(group_column, "Model")][1], 3)
)

pcoa_plot <- ggplot(pcoa_coordinates, aes(x = PCoA1, y = PCoA2, color = .data[[group_column]])) +
  geom_point(aes(color = .data[[group_column]]), size = 4.5) +
  stat_ellipse(
    aes(fill = .data[[group_column]]),
    geom = "polygon",
    level = ellipse_level_pcoa,
    linetype = 2,
    linewidth = 0.5,
    alpha = 0.25,
    show.legend = TRUE
  ) +
  geom_hline(yintercept = 0, colour = "#BEBEBE", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "#BEBEBE", linetype = "dashed") +
  scale_color_manual(values = group_palette[group_order]) +
  scale_fill_manual(values = group_palette[group_order]) +
  guides(
    fill = "none",
    color = guide_legend(override.aes = list(fill = NA, linetype = 0, size = 4))
  ) +
  labs(
    x = paste0("PCoA1 (", eig_percent[1], "%)"),
    y = paste0("PCoA2 (", eig_percent[2], "%)"),
    caption = adonis_caption
  ) +
  theme(
    legend.position = c(0.82, 0.85),
    legend.title = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.caption = element_text(hjust = 0, size = 9)
  )
save_ggplot_pair(pcoa_plot, file.path(figures_dir, pcoa_pdf), file.path(figures_dir, pcoa_png), ordination_width, ordination_height, png_dpi)

set.seed(random_seed)
nmds <- vegan::metaMDS(
  community_for_beta,
  distance = beta_distance_method,
  k = 2,
  trymax = nmds_trymax,
  autotransform = FALSE,
  trace = FALSE
)
nmds_coordinates <- data.frame(
  sample_id = rownames(nmds$points),
  MDS1 = nmds$points[, 1],
  MDS2 = nmds$points[, 2],
  stress = nmds$stress,
  stringsAsFactors = FALSE
)
nmds_coordinates <- merge(nmds_coordinates, metadata[, c("sample_id", group_column), drop = FALSE], by = "sample_id", sort = FALSE)
nmds_coordinates[[group_column]] <- factor(nmds_coordinates[[group_column]], levels = group_order)
write.csv(nmds_coordinates, file.path(results_dir, nmds_output), row.names = FALSE)

nmds_plot <- ggplot(nmds_coordinates, aes(x = MDS1, y = MDS2)) +
  geom_point(aes(color = .data[[group_column]]), shape = 19, size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8, color = "grey50") +
  stat_ellipse(
    aes(fill = .data[[group_column]]),
    geom = "polygon",
    level = ellipse_level_nmds,
    linetype = 2,
    linewidth = 0.5,
    alpha = 0.22
  ) +
  scale_color_manual(values = group_palette[group_order]) +
  scale_fill_manual(values = group_palette[group_order]) +
  guides(
    fill = "none",
    color = guide_legend(override.aes = list(fill = NA, linetype = 0, size = 4))
  ) +
  ggtitle(paste0("Stress = ", round(nmds$stress, 3))) +
  theme_bw(base_size = 12) +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
save_ggplot_pair(nmds_plot, file.path(figures_dir, nmds_pdf), file.path(figures_dir, nmds_png), ordination_width, ordination_height, png_dpi)

pairwise_distances <- make_pairwise_distance_table(beta_distance_matrix, metadata, group_column, group_order)
write.csv(pairwise_distances, file.path(results_dir, beta_pairwise_output), row.names = FALSE)

focal_distances <- make_focal_distance_table(beta_distance_matrix, metadata, group_column)
focal_distances$group <- factor(focal_distances$group, levels = group_order)
write.csv(focal_distances, file.path(results_dir, beta_focal_output), row.names = FALSE)

beta_comparisons <- list(c("Control", "Tailing"), c("Control", "Mining"), c("Control", "Smelting"))
beta_y_max <- max(focal_distances$bray_curtis, na.rm = TRUE)
beta_y_positions <- beta_y_max + seq(0.04, 0.12, length.out = length(beta_comparisons))

beta_box_plot <- ggplot(focal_distances, aes(x = group, y = bray_curtis, fill = group)) +
  gghalves::geom_half_violin(
    position = position_nudge(x = 0.25),
    side = "r",
    width = 0.8,
    color = NA,
    alpha = 0.78
  ) +
  geom_boxplot(width = 0.38, linewidth = 0.9, outlier.color = NA, alpha = 0.95) +
  geom_jitter(aes(fill = group), shape = 21, size = 2.1, width = 0.14, alpha = 0.72, color = "grey30") +
  ggsignif::geom_signif(
    comparisons = beta_comparisons,
    map_signif_level = TRUE,
    test = "wilcox.test",
    y_position = beta_y_positions,
    size = 0.8,
    color = "black",
    textsize = 3.5
  ) +
  scale_fill_manual(values = group_palette[group_order]) +
  scale_y_continuous(limits = c(0, max(beta_y_positions) + 0.06), expand = c(0, 0)) +
  labs(x = NULL, y = "Bray-Curtis distance to other samples") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1),
    axis.text.x = element_text(color = "black", size = 12, angle = 25, hjust = 1),
    axis.text.y = element_text(color = "black", size = 12),
    legend.position = "none",
    axis.ticks = element_line(color = "black", linewidth = 0.8)
  )
save_ggplot_pair(beta_box_plot, file.path(figures_dir, beta_box_pdf), file.path(figures_dir, beta_box_png), beta_box_width, beta_box_height, png_dpi)

message("Alpha/beta diversity demo completed.")
message("Rarefaction depth: ", ifelse(is.na(rarefaction_depth), "not applied", rarefaction_depth))
message("Beta distance method: ", beta_distance_method)
message("PERMANOVA: ", adonis_caption)
message("NMDS stress: ", round(nmds$stress, 3))
message("Outputs written to results/ and figures/.")
