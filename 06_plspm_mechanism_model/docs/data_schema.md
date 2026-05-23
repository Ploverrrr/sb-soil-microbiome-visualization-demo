# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are designed to mimic a metal-contaminated soil microbiome study.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match the other shared tables. |
| `group` | Site or treatment group. Preserved in output tables for interpretation. |

## `environmental_variables.csv`

Required columns depend on the user-editable `latent_variable_blocks`. The default model uses:

| Column | Default latent variable | Description |
| --- | --- | --- |
| `sample_id` | NA | Unique sample identifier. |
| `pH` | `Environment` | Simulated soil pH. |
| `EC` | `Environment` | Simulated electrical conductivity. |
| `TOC` | `Nutrient` | Simulated total organic carbon. |
| `TN` | `Nutrient` | Simulated total nitrogen. |
| `TP` | `Nutrient` | Simulated total phosphorus. |
| `Sb_total` | `Sb_As` | Simulated total antimony. |
| `As` | `Sb_As` | Simulated arsenic. |

Other environmental variables in the shared table can be used by editing `latent_variable_blocks` and `path_matrix`.

## `abundance_table.csv`

Required structure:

| Column | Description |
| --- | --- |
| `feature_id` | Unique microbial feature or ASV/OTU identifier. |
| sample columns | One column per `sample_id`, containing raw counts or pseudo-count abundance. |

The script derives alpha-diversity indicators internally:

| Derived column | Description |
| --- | --- |
| `Chao1` | Chao1-style richness estimator from pseudo-rarefied feature counts. |
| `Shannon` | Shannon diversity. |
| `Simpson` | Gini-Simpson diversity, `1 - sum(p^2)`. |
| `Pielou` | Shannon evenness divided by log richness. |

## `functional_annotation_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Sample identifier. |
| `function_category` | Functional category used for derived indicators. |
| `count` | Simulated functional pseudo-count. |

The default script derives:

| Derived column | Source `function_category` |
| --- | --- |
| `MetalRes` | `Metal resistance` |
| `Sulfur` | `Sulfur cycling` |
| `Nitrogen` | `Nitrogen cycling` |
| `Carbon` | `Carbon cycling` |

## Default PLS-PM Blocks

| Latent variable | Manifest variables |
| --- | --- |
| `Environment` | `pH`, `EC` |
| `Nutrient` | `TOC`, `TN`, `TP` |
| `Sb_As` | `Sb_total`, `As` |
| `Function` | `MetalRes`, `Sulfur`, `Nitrogen`, `Carbon` |
| `Alpha` | `Chao1`, `Shannon`, `Simpson`, `Pielou` |

## Output Tables

| Output | Description |
| --- | --- |
| `plspm_model_input_table.csv` | Unscaled derived model input table with sample metadata. |
| `plspm_scaled_input_table.csv` | Scaled manifest variables used by `plspm()`. |
| `plspm_path_matrix.csv` | 0/1 latent-variable path matrix. |
| `plspm_latent_variable_blocks.csv` | Long table linking latent variables to manifest variables. |
| `plspm_path_coefficients.csv` | Path coefficient matrix from `plspm_model$path_coefs`. |
| `plspm_inner_model.csv` | Flattened inner model statistics. |
| `plspm_inner_summary.csv` | Inner model summary from `plspm`. |
| `plspm_outer_model.csv` | Outer model loadings/weights from `plspm`. |
| `plspm_effects.csv` | Direct, indirect, and total effects. |
| `plspm_latent_scores.csv` | Sample-level latent variable scores. |
| `plspm_model_metrics.csv` | Goodness-of-fit and basic model dimensions. |
