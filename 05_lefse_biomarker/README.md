# 05_lefse_biomarker

This module demonstrates a reproducible LEfSe-style biomarker visualization workflow for a simulated metal-contaminated soil microbiome study.

The demo starts from shared raw-like toy input tables, calculates taxon relative abundance, performs group-wise biomarker screening, estimates a simple effect-size score, and produces publication-style biomarker figures. It does not use a precomputed LEfSe result table or a manually edited plotting table as input.

## What This Module Shows

- Taxonomic count table processing from feature-level pseudo-counts.
- Relative abundance calculation by sample.
- Aggregation to all major taxonomic levels by default, matching the original `taxa_level = "all"` LEfSe script logic.
- Group-wise biomarker screening using Kruskal-Wallis tests and one-vs-rest Wilcoxon tests.
- A LEfSe-style log-ratio effect score for ranked biomarker plotting.
- Group-colored biomarker barplot, taxonomic cladogram, KW abundance plot, combined figures, and auxiliary abundance heatmap.

This is a portfolio-friendly, simplified LEfSe-style workflow. It preserves the biomarker-discovery logic and figure style, but it is not a full replacement for the original LEfSe implementation.

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
- `top_n_taxa`
- `max_biomarkers_to_plot`
- `max_cladogram_labels`
- `min_prevalence`
- `min_mean_relative_abundance`
- `kruskal_fdr_cutoff`
- `effect_score_cutoff`
- output file names
- figure width and height
- group colors and heatmap palette

## Replacing With Your Own Data

Use the same input structure:

- `sample_metadata.csv` must contain `sample_id` and `group`.
- `abundance_table.csv` must contain `feature_id` plus one column per sample.
- `taxonomy_table.csv` must contain `feature_id` and the selected taxonomic ranks. The default all-rank workflow uses `Phylum`, `Class`, `Order`, `Family`, and `Genus`.

Keep sample IDs identical across metadata and abundance table column names. Keep feature IDs identical across abundance and taxonomy tables.

## Why Not Use A Hand-Edited LEfSe Table?

A precomputed `LDA_score.csv` or manually edited plotting table would only demonstrate plotting. This module recalculates abundance, statistics, effect scores, and plotting tables inside the demo so the workflow remains transparent and reproducible.

## Original-Script Features Preserved

The original reference script produced an LDA score barplot, a taxonomic cladogram, a combined LDA+cladogram figure, a KW abundance figure, and a combined LDA+KW figure. This module now mirrors that functional structure with simulated data and dependency-light R code.
