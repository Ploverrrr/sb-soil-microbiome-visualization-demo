# 05_lefse_biomarker

This module demonstrates a reproducible LEfSe biomarker visualization workflow for a simulated metal-contaminated soil microbiome study.

The demo starts from shared raw-like toy input tables, builds a `microeco::microtable`, runs `microeco::trans_diff$new(method = "lefse")`, and uses the native `microeco` plotting functions for the LDA barplot, cladogram, and KW abundance plot. It does not use a precomputed LEfSe result table or a manually edited plotting table as input.

## What This Module Shows

- Taxonomic count table processing from feature-level pseudo-counts.
- Relative abundance calculation by sample for transparent intermediate output.
- `microeco` LEfSe analysis with `taxa_level = "all"`, matching the original script logic.
- Native `plot_diff_bar()` LDA score visualization.
- Native `plot_diff_cladogram()` taxonomic cladogram visualization.
- Native `plot_diff_abund()` KW abundance visualization using the same selected taxa as the LDA plot.
- Group-colored combined figures and an auxiliary abundance heatmap derived from the recalculated model output.

This is a portfolio-friendly LEfSe module using simulated data. The cladogram is not a hand-drawn approximation; it is generated through the same `microeco` plotting family used by the original reference script.

## Shared Toy Inputs

The module reads these shared simulated tables:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/taxonomy_table.csv
```

The toy data are simulated/desensitized and are only intended to demonstrate the workflow. They do not represent real study results.

## How To Run

From this module directory:

```bash
cd 05_lefse_biomarker
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## Main Outputs

Results:

- `results/taxon_relative_abundance_by_sample.csv`
- `results/lefse_candidate_statistics.csv`
- `results/lefse_biomarker_table.csv`
- `results/lefse_barplot_plotting_table.csv`
- `results/kw_abundance_plotting_table.csv`
- `results/cladogram_node_table.csv`
- `results/cladogram_edge_table.csv`
- `results/cladogram_label_table.csv`
- `results/biomarker_group_heatmap_table.csv`

Figures:

- `figures/lefse_biomarker_barplot.pdf`
- `figures/lefse_biomarker_barplot.png`
- `figures/lefse_cladogram.pdf`
- `figures/lefse_cladogram.png`
- `figures/lefse_barplot_cladogram_combined.pdf`
- `figures/lefse_barplot_cladogram_combined.png`
- `figures/lefse_kw_abundance_plot.pdf`
- `figures/lefse_kw_abundance_plot.png`
- `figures/lefse_lda_kw_combined.pdf`
- `figures/lefse_lda_kw_combined.png`
- `figures/biomarker_group_heatmap.pdf`
- `figures/biomarker_group_heatmap.png`

## User-Editable Settings

Edit the settings block at the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `target_taxonomic_level`
- `taxonomic_levels_for_lefse`
- `group_order`
- `lefse_alpha`
- `p_adjust_method`
- `lefse_norm`
- `lda_threshold`
- `plot_feature_count`
- `cladogram_use_taxa_num`
- `cladogram_use_feature_num`
- `clade_label_level`
- `cladogram_filter_taxa`
- output file names
- figure width and height
- group colors and heatmap palette

## Replacing With Your Own Data

Use the same input structure:

- `sample_metadata.csv` must contain `sample_id` and `group`.
- `abundance_table.csv` must contain `feature_id` plus one column per sample.
- `taxonomy_table.csv` must contain `feature_id` and the selected taxonomic ranks. The default all-rank workflow uses `Phylum`, `Class`, `Order`, `Family`, and `Genus`.

Keep sample IDs identical across metadata and abundance table column names. Keep feature IDs identical across abundance and taxonomy tables.

The script adds `k__`, `p__`, `c__`, `o__`, `f__`, and `g__` prefixes internally when creating the `microeco` object, so users do not need to prefix the shared CSV manually.

## Why Not Use A Hand-Edited LEfSe Table?

A precomputed `LDA_score.csv`, `res_diff.csv`, or manually edited plotting table would only demonstrate plotting. This module recalculates the LEfSe result with `microeco::trans_diff$new()` and then exports model-derived plotting tables so the workflow remains transparent and reproducible.

## Original-Script Features Preserved

The original reference script produced an LDA score barplot, a taxonomic cladogram, a combined LDA+cladogram figure, a KW abundance figure, and a combined LDA+KW figure. This module mirrors that functional structure with simulated data and the original `microeco` plotting workflow: `plot_diff_bar()`, `plot_diff_cladogram()`, `plot_diff_abund()`, `aplot::insert_left()`, and `gridExtra::arrangeGrob()`.
