# Data Schema

This module reads shared simulated/desensitized toy data from `../data/toy_shared/`.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match abundance, functional, and environmental tables. |
| `group` | Site group. Used for sample alignment and retained for user extension. |

## `environmental_variables.csv`

Required columns are controlled by the settings in `scripts/run_demo.R`.

Default VPA groups:

| Group | Default columns | Purpose |
| --- | --- | --- |
| Sb | `Sb_total` | Focal antimony predictor group. |
| Cu & As | `Cu`, `As` | Co-metal predictor group. |
| Nutrients | `TOC`, `SO4`, `pH` | Nutrient/background chemistry predictor group. |

Default Mantel heatmap variables:

`pH`, `EC`, `TOC`, `TN`, `SO4`, `NO3`, `Sb_total`, `Sb_III`, `Sb_V`, `As`, `Cu`

## `abundance_table.csv`

Required format:

| Column | Description |
| --- | --- |
| `feature_id` | ASV/OTU/feature identifier. |
| sample columns | Numeric pseudo-count abundance for each sample. |

The module uses the top configured features to build the species/community response matrix.

## `functional_annotation_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Sample identifier. |
| `ko_id` | Default functional feature used to build the functional response matrix. |
| `abundance` | Default numeric functional abundance used for VPA/Mantel. |

The feature and value columns can be changed with `functional_feature_column` and `functional_value_column`.

## Sanity Checks

The script checks that:

- all input files exist;
- required columns are present;
- sample IDs can be aligned across shared tables;
- abundance values are numeric and non-negative;
- selected environmental variables are numeric.
