# 01 RF Correlation Heatmap

This module demonstrates a publication-style environmental association figure combining:

- correlation heatmap between top microbial features and soil environmental variables;
- random forest feature importance for predicting one target environmental variable.

The demo reads the repository-level shared toy dataset. It does not use private data, real sample values, or precomputed `Importance(%).CSV`-style tables.

## What This Module Shows

The workflow starts from raw-like toy tables and recomputes all analysis results inside the module:

1. calculate sample-wise relative abundance;
2. aggregate features to a selected taxonomic level, such as `Genus`;
3. select top features by mean relative abundance;
4. calculate feature-environment correlations;
5. train a random forest model to predict a selected environmental variable;
6. visualize correlation strength and RF importance.

By default, the module uses top genera and predicts `Sb_total`.

## Shared Toy Inputs

The default input files are:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/environmental_variables.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/taxonomy_table.csv
../data/toy_shared/functional_annotation_table.csv
```

The module can use either:

- `feature_source = "taxonomy"`: aggregate microbial features from `abundance_table.csv` and `taxonomy_table.csv`;
- `feature_source = "function"`: aggregate functional rows from `functional_annotation_table.csv`.

The `data/toy/` folder inside this module is reserved for optional module-specific inputs and should not duplicate the shared toy dataset.

## How To Run

From the repository root:

```bash
cd 01_rf_correlation_heatmap
Rscript scripts/run_demo.R
```

Outputs are written only inside this module:

```text
results/
figures/
```

The primary portfolio figure is:

```text
figures/rf_correlation_combined.pdf
figures/rf_correlation_combined.png
```

The separate heatmap and RF importance plots are retained as auxiliary outputs for inspection and reuse:

```text
figures/correlation_heatmap.pdf
figures/correlation_heatmap.png
figures/rf_importance_plot.pdf
figures/rf_importance_plot.png
```

## User-Editable Parameters

Open `scripts/run_demo.R` and edit the settings block near the top. Common settings include:

- `feature_source`
- `target_taxonomic_level`
- `environmental_variables`
- `target_environmental_variable_for_rf`
- `top_n_features`
- `correlation_method`
- `p_adjust_method`
- output file names
- figure size
- heatmap colors
- RF bubble size range

## Interpreting The Outputs

Correlation values describe monotonic or linear associations between feature abundance and environmental variables, depending on `correlation_method`.

Random forest importance describes how much each selected feature contributes to predicting the chosen target environmental variable in a nonlinear regression model. The reported model metric is the out-of-bag pseudo R2 from the fitted random forest model; it is calculated from the model and is not manually hard-coded.

The combined figure uses the same selected features and the same feature order for both panels. The left panel shows feature-environment correlations; the right panel shows RF importance recalculated from the model.

## Why Not Use A Precomputed Importance Table?

A hand-edited `Importance(%).CSV` table is useful during figure polishing, but it hides the modeling step. This demo recomputes feature abundance, correlations, p-value adjustment, random forest importance, and model performance from raw-like toy inputs so the workflow remains reproducible.

## Using Your Own Data

To use your own data, preserve the shared schema or update the file paths and column settings in `scripts/run_demo.R`.

Required relationships:

- `sample_id` values must match across metadata, environmental variables, and abundance sample columns.
- `feature_id` values must match between abundance and taxonomy tables when using `feature_source = "taxonomy"`.
- functional rows must include `sample_id`, `feature_id`, `pathway`, `count`, and related annotation columns when using `feature_source = "function"`.

The shared toy data are simulated/desensitized and do not represent real research results.
