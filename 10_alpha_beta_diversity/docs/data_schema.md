# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are used only to demonstrate a reproducible workflow.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match sample columns in `abundance_table.csv`. |
| `group` | Experimental or site group used for boxplots, PERMANOVA, PCoA/NMDS colors, and beta-distance comparisons. |

Optional columns from the shared toy dataset, such as `site_type`, `contamination_level`, `replicate`, and `batch`, are not required by this module.

## `abundance_table.csv`

Required format:

| Column | Description |
| --- | --- |
| `feature_id` | ASV/OTU/feature identifier. |
| sample columns | One numeric pseudo-count column per `sample_id`. |

The table is expected to be feature by sample. The script transposes it to sample by feature for `vegan`.

## Sanity Checks

The script checks that:

- input files exist;
- required columns are present;
- all sample IDs in `sample_metadata.csv` and `abundance_table.csv` match;
- abundance values are numeric and non-negative;
- configured groups are present in the metadata.

## Supported Outputs From These Inputs

From these two raw-like tables, the module computes:

- rarefied count matrix for alpha diversity, if `rarefy_counts = TRUE`;
- alpha diversity metrics;
- pairwise alpha diversity tests;
- relative abundance matrix for Bray-Curtis beta diversity;
- PCoA coordinates;
- PERMANOVA results;
- NMDS coordinates;
- pairwise and focal beta-distance tables.
