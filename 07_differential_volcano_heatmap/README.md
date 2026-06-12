# 07_differential_volcano_heatmap

This module demonstrates a reproducible differential abundance workflow for a simulated metal-contaminated soil microbiome study.

The demo follows a reference workflow pattern: DESeq2 differential analysis, `ggplot2` volcano plots with `ggrepel` labels, `patchwork` combined volcano panels, and `ComplexHeatmap`/`circlize` heatmap visualization. It does not use a precomputed differential result table or a manually edited heatmap matrix as input.

## What This Module Shows

- Aggregation of feature-level pseudo-counts to a selected taxonomic level.
- DESeq2 comparisons of each contaminated group against the control group.
- Up/Down/Stable classification using log2 fold-change and p-value thresholds.
- Demo volcano plot layout with cyan/grey/pink colors, threshold lines, and labeled top taxa.
- Rectangular `ComplexHeatmap` output and circular `circos.heatmap()` output for selected differential taxa.

## Shared Toy Inputs

The module reads:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/taxonomy_table.csv
```

The toy data are simulated/desensitized and are only intended to demonstrate the workflow. They do not represent real study results.

## How To Run

From this module directory:

```bash
cd 07_differential_volcano_heatmap
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## Main Outputs

Results:

- `results/taxon_count_matrix.csv`
- `results/deseq2_all_contrasts.csv`
- `results/deseq2_tailing_vs_control.csv`
- `results/deseq2_mining_vs_control.csv`
- `results/deseq2_smelting_vs_control.csv`
- `results/heatmap_zscore_matrix.csv`
- `results/heatmap_selected_taxa.csv`

Figures:

- `figures/volcano_tailing_vs_control.pdf`
- `figures/volcano_tailing_vs_control.png`
- `figures/volcano_mining_vs_control.pdf`
- `figures/volcano_mining_vs_control.png`
- `figures/volcano_smelting_vs_control.pdf`
- `figures/volcano_smelting_vs_control.png`
- `figures/deseq2_volcano_combined.pdf`
- `figures/deseq2_volcano_combined.png`
- `figures/differential_heatmap.pdf`
- `figures/differential_heatmap.png`
- `figures/differential_circular_heatmap.pdf`
- `figures/differential_circular_heatmap.png`

## User-Editable Settings

Edit the settings block at the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `control_group`
- `treatment_groups`
- `target_taxonomic_level`
- `deseq_fit_type`
- `log2fc_threshold`
- `pvalue_threshold`
- `top_label_n`
- `volcano_y_axis_cap`
- `heatmap_focus_group`
- `heatmap_top_n_each_direction`
- output file names
- figure width and height
- volcano and heatmap colors

## Replacing With Your Own Data

Use the same shared input structure:

- `sample_metadata.csv` must contain `sample_id` and `group`.
- `abundance_table.csv` must contain `feature_id` plus one column per sample.
- `taxonomy_table.csv` must contain `feature_id` and the selected `target_taxonomic_level`.

Keep sample IDs identical across metadata and abundance table column names. Keep feature IDs identical across abundance and taxonomy tables.

## Why Not Use A Finished Differential Table?

A precomputed `DESeq2_result.csv` or hand-edited heatmap matrix would only demonstrate plotting. This module recalculates the taxon count matrix, DESeq2 results, selected labels, and heatmap matrix from raw-like toy input so the workflow remains transparent and reproducible.

## Reference Workflow Features

The reference workflow uses `DESeq2`, volcano plots with `geom_point()`, `scale_color_manual()`, dashed threshold lines, `theme_bw()`, `ggrepel` labels for top features, `patchwork` for combined volcano panels, and `ComplexHeatmap`/`circlize` for heatmap visualization. This module keeps those plotting methods and visual grammar while using simulated shared toy data.
