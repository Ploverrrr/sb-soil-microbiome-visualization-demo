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

`target_taxonomic_level` in `scripts/run_demo.R` selects which taxonomy column is used for aggregation. The default is `all`, which aggregates `Phylum`, `Class`, `Order`, `Family`, and `Genus` to better match the original LEfSe script.

## Output Tables

| Output | Description |
| --- | --- |
| `taxon_relative_abundance_by_sample.csv` | Selected taxonomic level by sample relative abundance table. |
| `lefse_candidate_statistics.csv` | All tested taxa with prevalence, group means, Kruskal-Wallis FDR, Wilcoxon p value, and effect score. |
| `lefse_biomarker_table.csv` | Filtered biomarkers used for final figures. |
| `lefse_barplot_plotting_table.csv` | Plot-ready biomarker table for the LEfSe-style barplot. |
| `kw_abundance_plotting_table.csv` | Plot-ready group mean abundance table for the KW-style abundance plot. |
| `cladogram_node_table.csv` | Node coordinates, abundance, and enrichment groups for the simplified taxonomic cladogram. |
| `cladogram_edge_table.csv` | Parent-child edge coordinates for the simplified taxonomic cladogram. |
| `biomarker_group_heatmap_table.csv` | Plot-ready group mean abundance table for the heatmap. |
