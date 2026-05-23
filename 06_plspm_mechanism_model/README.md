# 06_plspm_mechanism_model

This module demonstrates a reproducible PLS-PM mechanism-model workflow for a simulated metal-contaminated soil microbiome study.

The workflow follows the original reference scripts: it builds a manifest-variable table, scales it, defines latent-variable blocks and a 0/1 path matrix, runs `plspm::plspm()`, and draws the native `plspm` inner and outer model plots. The demo does not use a precomputed model table or manually edited path diagram.

## What This Module Shows

- Construction of a mechanism-oriented PLS-PM input table from shared raw-like toy data.
- Environmental, nutrient, contaminant, functional, and alpha-diversity latent variables.
- Native `plspm::innerplot()` path diagram.
- Native `plspm::outerplot(what = "loadings")` and `outerplot(what = "weights")` measurement-model diagrams.
- Exported path coefficients, inner model statistics, outer model loadings/weights, effects, latent scores, and model metrics.

## Shared Toy Inputs

The module reads:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/environmental_variables.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/functional_annotation_table.csv
```

The toy data are simulated/desensitized and are only intended to demonstrate the workflow. They do not represent real study results.

## How To Run

From this module directory:

```bash
cd 06_plspm_mechanism_model
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## Main Outputs

Results:

- `results/plspm_model_input_table.csv`
- `results/plspm_scaled_input_table.csv`
- `results/plspm_path_matrix.csv`
- `results/plspm_latent_variable_blocks.csv`
- `results/plspm_path_coefficients.csv`
- `results/plspm_inner_model.csv`
- `results/plspm_inner_summary.csv`
- `results/plspm_outer_model.csv`
- `results/plspm_effects.csv`
- `results/plspm_latent_scores.csv`
- `results/plspm_model_metrics.csv`

Figures:

- `figures/plspm_inner_path_model.pdf`
- `figures/plspm_inner_path_model.png`
- `figures/plspm_outer_loadings.pdf`
- `figures/plspm_outer_loadings.png`
- `figures/plspm_outer_weights.pdf`
- `figures/plspm_outer_weights.png`
- `figures/plspm_total_effects.pdf`
- `figures/plspm_total_effects.png`

## User-Editable Settings

Edit the settings block at the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `latent_variable_blocks`
- `path_matrix`
- `plspm_modes`
- positive/negative path colors
- output file names
- figure width and height

## Replacing With Your Own Data

Use the same shared input structure:

- `sample_metadata.csv` must contain `sample_id` and `group`.
- `environmental_variables.csv` must contain `sample_id` and the environmental indicators used in `latent_variable_blocks`.
- `abundance_table.csv` must contain `feature_id` plus one column per sample. The script calculates Chao1, Shannon, Simpson, and Pielou internally.
- `functional_annotation_table.csv` must contain `sample_id`, `function_category`, and `count`. The script derives functional indicators for metal resistance, sulfur cycling, nitrogen cycling, and carbon cycling.

If you rename indicators or latent variables, update both `latent_variable_blocks` and `path_matrix` so their names and order remain identical.

## Original-Script Features Preserved

The original PLS-PM scripts used `plspm`, scaled the model input, defined manifest-variable blocks, created a 0/1 latent-variable path matrix, and drew `innerplot()` plus `outerplot()` figures with red positive paths, blue negative paths, grey links, displayed path values, and no inner box border. This module preserves that plotting workflow while replacing private inputs with simulated shared toy data.
