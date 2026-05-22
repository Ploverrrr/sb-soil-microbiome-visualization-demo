# Shared Toy Data Schema

This repository uses one shared toy dataset as the common input foundation for all numbered demo modules. The dataset simulates a small metal-contaminated soil microbiome study with four groups:

- `Control`
- `Tailing`
- `Mining`
- `Smelting`

Each group has six simulated samples, for a total of 24 samples.

These data are simulated, toy, and desensitized. They are designed to demonstrate reproducible visualization workflows and do not represent real research results, original sample values, or private intermediate analysis tables.

## Location

The shared dataset is stored in:

```text
data/toy_shared/
```

Expected files:

```text
sample_metadata.csv
environmental_variables.csv
taxonomy_table.csv
abundance_table.csv
functional_annotation_table.csv
```

## `sample_metadata.csv`

Sample-level grouping and study design metadata.

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. Must match sample columns in `abundance_table.csv` and `sample_id` values in other tables. |
| `group` | Experimental/site group: `Control`, `Tailing`, `Mining`, or `Smelting`. |
| `site_type` | Broad site class, such as `reference` or `metal_contaminated`. |
| `contamination_level` | Simulated contamination category. |
| `replicate` | Replicate number within group. |
| `batch` | Simulated sequencing or processing batch. |

Supported modules: all modules.

## `environmental_variables.csv`

Simulated soil chemistry and contaminant measurements.

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. |
| `pH` | Simulated soil pH. |
| `EC` | Simulated electrical conductivity. |
| `TOC` | Simulated total organic carbon. |
| `TN` | Simulated total nitrogen. |
| `TP` | Simulated total phosphorus. |
| `SO4` | Simulated sulfate. |
| `NO3` | Simulated nitrate. |
| `Sb` | Simulated total antimony. |
| `Sb3` | Simulated trivalent antimony. |
| `Sb5` | Simulated pentavalent antimony. |
| `As` | Simulated arsenic. |
| `Cu` | Simulated copper. |
| `Zn` | Simulated zinc. |
| `Cd` | Simulated cadmium. |
| `Fe` | Simulated iron. |
| `Mn` | Simulated manganese. |

Supported modules: `01_rf_correlation_heatmap`, `02_microbe_env_network`, `06_plspm_mechanism_model`, `09_sulfur_gene_contaminant_association`, `10_alpha_beta_diversity`, and `11_vpa_mantel_partitioning`.

## `taxonomy_table.csv`

Feature-level taxonomic annotation table.

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. Must match `feature_id` in `abundance_table.csv`. |
| `Kingdom` | Taxonomic kingdom. |
| `Phylum` | Taxonomic phylum. |
| `Class` | Taxonomic class. |
| `Order` | Taxonomic order. |
| `Family` | Taxonomic family. |
| `Genus` | Taxonomic genus. |
| `taxon_label` | Human-readable toy taxon label for plotting. |

Supported modules: `02_microbe_env_network`, `03_ternary_taxa_distribution`, `04_faprotax_functional_profile`, `05_lefse_biomarker`, `07_differential_volcano_heatmap`, and `10_alpha_beta_diversity`.

## `abundance_table.csv`

Feature-by-sample pseudo-count table.

| Column | Description |
|---|---|
| `feature_id` | Feature identifier. |
| sample columns | One column per `sample_id`; values are simulated pseudo-count abundances. |

The simulated abundance table includes mild but clear group trends, including taxa enriched in `Control`, `Tailing`, `Mining`, and `Smelting`, plus broadly distributed and rare taxa.

Supported modules: `01_rf_correlation_heatmap`, `02_microbe_env_network`, `03_ternary_taxa_distribution`, `04_faprotax_functional_profile`, `05_lefse_biomarker`, `07_differential_volcano_heatmap`, `10_alpha_beta_diversity`, and `11_vpa_mantel_partitioning`.

## `functional_annotation_table.csv`

Long-format functional annotation and pseudo-count table.

| Column | Description |
|---|---|
| `gene_id` | Simulated gene identifier. |
| `feature_id` | Linked microbial feature identifier. |
| `sample_id` | Sample identifier. |
| `annotation_source` | Simulated annotation source, such as `KEGG`, `FAPROTAX`, `BacMet`, `SulfurCycle`, `NitrogenCycle`, or `CarbonCycle`. |
| `function_category` | Broad functional category. |
| `pathway` | Pathway or function name. |
| `ko_id` | KO-like or database-like identifier. |
| `count` | Simulated functional pseudo-count. |
| `abundance` | Sample-wise relative abundance of the simulated functional count. |

Supported modules: `04_faprotax_functional_profile`, `06_plspm_mechanism_model`, `07_differential_volcano_heatmap`, `08_kegg_enrichment`, `09_sulfur_gene_contaminant_association`, and `11_vpa_mantel_partitioning`.

## Replacing the Toy Data

Users can replace the shared toy data with their own data if they preserve the same column names and relationships:

1. `sample_id` must be consistent across metadata, environmental variables, abundance columns, and functional annotation rows.
2. `feature_id` must be consistent across abundance, taxonomy, and functional annotation tables.
3. Abundance tables should be feature-by-sample, with `feature_id` in the first column.
4. Functional annotation should remain in long format.

The numbered modules should read from `data/toy_shared/` by default, but their scripts may expose user-editable path settings for custom data.
