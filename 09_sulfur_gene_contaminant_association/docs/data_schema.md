# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are designed to mimic a metal-contaminated soil microbiome study.

## `environmental_variables.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match functional annotation sample IDs. |
| `Sb_III` | Simulated trivalent antimony. Replaces the original script's `Sb3` name. |
| `Sb_V` | Simulated pentavalent antimony. Replaces the original script's `Sb5` name. |
| `Sb_total` | Simulated total antimony. |
| `SO4` | Simulated sulfate. |
| `As` | Simulated arsenic. |
| `Cu` | Simulated copper. |

The selected contaminant columns can be edited with `contaminant_vars` in `scripts/run_demo.R`.

## `functional_annotation_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Sample identifier. |
| `ko_id` | KO-like identifier used to identify sulfur genes. |
| `abundance` | Default value column used for relative sulfur gene abundance. |
| `count` | Optional alternative value column. |

The default target KO mapping is:

| KO ID | Gene label |
| --- | --- |
| `K11180` | `dsrA` |
| `K11181` | `dsrB` |
| `K17218` | `soxB` |
| `K17222` | `soxC` |
| `K17223` | `soxD` |
| `K00958` | `sat` |

## Derived Variables

| Variable | Description |
| --- | --- |
| `dsrAB` | Mean of available `dsrA` and `dsrB` columns. |
| `soxBCD` | Mean of available `soxB`, `soxC`, and `soxD` columns. |
| `sat` | Sulfate assimilation marker when `K00958` is available. |

If some target KOs are absent in the toy data, the script keeps available indicators and records the resulting variables in output tables.

## Output Tables

| Output | Description |
| --- | --- |
| `sulfur_gene_abundance_by_sample.csv` | Sample-level target sulfur gene abundance table. |
| `sulfur_gene_contaminant_analysis_table.csv` | Merged sulfur gene and environmental variable table used for analysis. |
| `pearson_results.csv` | Gene-by-contaminant Pearson correlation coefficients, p-values, and significance labels. |
| `mlr_results.csv` | Multiple linear regression coefficients and model summaries for selected Sb species. |
