# Data Schema

This module reads shared toy data from:

```text
../data/toy_shared/
```

The standard run mode is:

```bash
cd 04_faprotax_functional_profile
Rscript scripts/run_demo.R
```

## `sample_metadata.csv`

Required columns:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier. Must match `sample_id` in the functional annotation table. |
| `group` | Sample group used for group-level summaries and plotting. |

## `functional_annotation_table.csv`

Required columns:

| Column | Description |
|---|---|
| `gene_id` | Simulated gene identifier. |
| `feature_id` | Linked microbial feature identifier. |
| `sample_id` | Sample identifier. |
| `annotation_source` | Functional annotation source, such as `FAPROTAX`, `SulfurCycle`, `NitrogenCycle`, or `CarbonCycle`. |
| `function_category` | Broad functional category. |
| `pathway` | Pathway or function name used as the default function label. |
| `ko_id` | KO-like identifier. |
| `count` | Simulated functional pseudo-count. |
| `abundance` | Optional precomputed abundance. This module recalculates relative abundance from `count`. |

## Output Tables

| Output | Description |
|---|---|
| `functional_relative_abundance_by_sample.csv` | Function-by-sample relative abundance matrix. |
| `functional_group_summary.csv` | Mean, SD, SE, and sample count by group and function. |
| `functional_bubble_plotting_table.csv` | Plot-ready sample/function table for the bubble profile. |
| `functional_barplot_table.csv` | Plot-ready group/function table for the grouped barplot. |

## Notes

The shared functional annotation table is simulated. It is designed to support reproducible visualization workflows and does not represent real FAPROTAX database output.
