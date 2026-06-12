# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are designed to mimic a metal-contaminated soil microbiome study.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match sample columns in `abundance_table.csv`. |
| `group` | Experimental or site group used for biomarker testing. Default groups are `Control`, `Tailing`, `Mining`, and `Smelting`. |

Optional columns such as `site_type` and `replicate` can be present and are preserved in the shared dataset.

## `abundance_table.csv`

Required structure:

| Column | Description |
| --- | --- |
| `feature_id` | Unique microbial feature or ASV/OTU identifier. |
| sample columns | One column per `sample_id`, containing raw counts or pseudo-count abundance. |

The script converts each sample column to relative abundance internally.

## `taxonomy_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `feature_id` | Feature identifier matching `abundance_table.csv`. |
| `Kingdom` | Taxonomic kingdom. |
| `Phylum` | Taxonomic phylum. |
| `Class` | Taxonomic class. |
| `Order` | Taxonomic order. |
| `Family` | Taxonomic family. |
| `Genus` | Taxonomic genus. |
| `taxon_label` | Human-readable feature-level label. |

`target_taxonomic_level` in `scripts/run_demo.R` selects the LEfSe taxonomy level. The default is `all`, matching the reference `microeco::trans_diff$new(taxa_level = "all")` workflow pattern.

The shared CSV keeps clean taxonomy names without prefixes. The script adds `k__`, `p__`, `c__`, `o__`, `f__`, and `g__` prefixes internally before creating the `microeco::microtable`, because `microeco::plot_diff_cladogram()` expects prefixed taxonomy labels for native cladogram rendering.

## Output Tables

| Output | Description |
| --- | --- |
| `taxon_relative_abundance_by_sample.csv` | Selected taxonomic level by sample relative abundance table. |
| `lefse_candidate_statistics.csv` | Complete `microeco` LEfSe `res_diff` table exported after recalculation. |
| `lefse_biomarker_table.csv` | Top recalculated biomarkers used for demo figures. |
| `lefse_barplot_plotting_table.csv` | Native `plot_diff_bar()` plotting data. |
| `kw_abundance_plotting_table.csv` | Native `plot_diff_abund()` plotting data. |
| `cladogram_node_table.csv` | Native `plot_diff_cladogram()` data exported from the returned plot object. |
| `cladogram_edge_table.csv` | Parent-child cladogram fields exported from the native plot object where available. |
| `cladogram_label_table.csv` | Letter annotation fields exported from the native cladogram plot object where available. |
| `biomarker_group_heatmap_table.csv` | Model-derived group mean abundance table for the auxiliary heatmap. |
