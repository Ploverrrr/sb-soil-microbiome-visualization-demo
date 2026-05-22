# Data Schema

This module reads the shared toy dataset from:

```text
../data/toy_shared/
```

The standard run mode is:

```bash
cd 01_rf_correlation_heatmap
Rscript scripts/run_demo.R
```

## `sample_metadata.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. Must match abundance sample columns and environmental rows. |
| `group` | Site group. Used for optional interpretation and downstream extension. |

## `environmental_variables.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. |
| user-selected environmental variables | Numeric variables used for correlation and RF modeling. |

The default environmental variables are:

```text
pH, EC, TOC, TN, TP, SO4, NO3, Sb_total, Sb_III, Sb_V, As, Cu, Zn, Cd, Fe, Mn
```

The default RF response variable is `Sb_total`.

## `abundance_table.csv`

Required structure for `feature_source = "taxonomy"`:

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. |
| sample columns | One column per `sample_id`; values are pseudo-count abundances. |

The script converts pseudo-counts to sample-wise relative abundance.

## `taxonomy_table.csv`

Required columns for `feature_source = "taxonomy"`:

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. Must match `abundance_table.csv`. |
| `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus` | Taxonomic levels available for aggregation. |
| `taxon_label` | Optional display label. |

The `target_taxonomic_level` setting controls which taxonomy column is used for aggregation.

## `functional_annotation_table.csv`

Required columns for `feature_source = "function"`:

| Column | Description |
|---|---|
| `gene_id` | Simulated gene identifier. |
| `feature_id` | Linked microbial feature identifier. |
| `sample_id` | Sample identifier. |
| `annotation_source` | Functional annotation source. |
| `function_category` | Broad function category. |
| `pathway` | Pathway or function label used as the default functional feature. |
| `ko_id` | KO-like identifier. |
| `count` | Simulated functional pseudo-count. |
| `abundance` | Optional precomputed relative abundance; this module recalculates from `count`. |

## Output Tables

The module writes:

| Output | Description |
|---|---|
| `relative_abundance_by_feature.csv` | Feature-level relative abundance by sample for taxonomy mode. |
| `selected_feature_abundance.csv` | Top feature abundance matrix used for correlation and RF. |
| `correlation_results.csv` | Feature-environment correlation coefficients, p values, adjusted p values, and significance labels. |
| `rf_importance.csv` | Random forest feature importance and scaled importance. |
| `rf_model_performance.csv` | Out-of-bag pseudo R2 and final MSE from the RF model. |
| `heatmap_plotting_table.csv` | Plot-ready correlation table with RF bubble sizes. |
