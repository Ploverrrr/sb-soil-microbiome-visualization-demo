# Data Schema

This module reads three shared toy input tables from:

```text
../data/toy_shared/
```

The standard run mode is:

```bash
cd 03_ternary_taxa_distribution
Rscript scripts/run_demo.R
```

## `sample_metadata.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. Must match the sample columns in `abundance_table.csv`. |
| `group` | Site or treatment group. The default ternary axes use `Control`, `Mining`, and `Smelting`. |

Additional columns such as `site_type`, `contamination_level`, `replicate`, and `batch` may be present and are ignored by this module.

## `abundance_table.csv`

Required structure:

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. Must match `taxonomy_table.csv`. |
| sample columns | One column per sample ID. Values are pseudo-count abundances. |

The script converts pseudo-counts to sample-wise relative abundance before aggregation.

## `taxonomy_table.csv`

Required columns:

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. Must match `abundance_table.csv`. |
| `Kingdom` | Taxonomic kingdom. |
| `Phylum` | Taxonomic phylum. |
| `Class` | Taxonomic class. |
| `Order` | Taxonomic order. |
| `Family` | Taxonomic family. |
| `Genus` | Taxonomic genus. |
| `taxon_label` | Optional display label. |

The `target_taxonomic_level` setting in `scripts/run_demo.R` can be set to any taxonomy column present in this table, such as `Phylum`, `Family`, or `Genus`.

## Output Tables

The module writes intermediate tables to `results/`:

| Output | Description |
|---|---|
| `relative_abundance_by_feature.csv` | Feature-level relative abundance by sample. |
| `group_mean_abundance_by_<level>.csv` | Mean relative abundance by taxon and group. |
| `ternary_plotting_table_<level>.csv` | Top taxa and ternary proportions used for plotting. |
