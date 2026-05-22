# Data Schema

This module reads shared toy data from:

```text
../data/toy_shared/
```

The standard run mode is:

```bash
cd 02_microbe_env_network
Rscript scripts/run_demo.R
```

## `sample_metadata.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. Must match sample columns in the abundance table and rows in the environmental table. |
| `group` | Sample group. Used for study context and optional filtering in future extensions. |

## `environmental_variables.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. |
| user-selected environmental variables | Numeric variables used as environmental nodes in the network. |

Default environmental variables include pH, nutrients, sulfate, and metal contaminants such as `Sb_total`, `As`, `Cu`, `Zn`, and `Cd`.

## `abundance_table.csv`

Required structure:

| Column | Description |
|---|---|
| `feature_id` | Microbial feature identifier. |
| sample columns | One column per sample ID. Values are pseudo-count abundances. |

The script converts counts to sample-wise relative abundance before network analysis.

## `taxonomy_table.csv`

Required columns:

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. Must match `abundance_table.csv`. |
| `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus` | Taxonomic levels available for aggregation. |
| `taxon_label` | Optional display label. |

The default aggregation level is `Genus`.

## Output Tables

| Output | Description |
|---|---|
| `relative_abundance_by_feature.csv` | Feature-level relative abundance by sample. |
| `taxon_abundance_by_sample.csv` | Taxon-level abundance matrix after aggregation. |
| `microbe_env_correlation_results.csv` | Long-format correlation, p value, adjusted p value, and sign table. |
| `microbe_env_corr_matrix.csv` | Wide matrix of filtered correlations. |
| `network_edge_list.csv` | Edge list for the microbe-environment network. |
| `network_node_list.csv` | Node attributes for taxa and environmental variables. |
| `network_summary.csv` | Small summary table with node and edge counts. |
