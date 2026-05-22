# 03 Ternary Taxa Distribution

This module demonstrates a publication-style ternary plot for visualizing how microbial taxa are distributed across three site groups in a simulated metal-contaminated soil microbiome study.

The demo reads the shared toy dataset from the repository-level `data/toy_shared/` folder. It does not use private data, real sample values, or a precomputed final plotting table.

## What This Module Shows

The workflow starts from raw-like toy tables:

- sample metadata
- feature-by-sample pseudo-count abundance
- feature taxonomy

It then calculates relative abundance, aggregates features to a user-selected taxonomic level, computes group mean abundance, selects the top taxa, converts the three selected groups into ternary proportions, and draws a ternary plot.

In the final plot:

- each point is a taxon;
- the three axes represent `Control`, `Mining`, and `Smelting` by default;
- point size represents overall mean relative abundance;
- point color identifies the taxon.

## Input Data

By default, the module reads:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/taxonomy_table.csv
```

The `data/toy/` folder inside this module is reserved for optional module-specific input and is not used to duplicate the shared toy dataset.

## How To Run

From the repository root:

```bash
cd 03_ternary_taxa_distribution
Rscript scripts/run_demo.R
```

Outputs are written only inside this module:

```text
results/
figures/
```

## User-Editable Parameters

Open `scripts/run_demo.R` and edit the parameter block near the top. Common settings include:

- `shared_data_dir`
- `target_taxonomic_level`
- `group_for_axis_1`
- `group_for_axis_2`
- `group_for_axis_3`
- `top_n`
- output file names
- figure size
- color palette

## Using Your Own Data

To use your own data, either replace the files in `../data/toy_shared/` with files using the same schema, or update the file paths in the user-editable settings section of `scripts/run_demo.R`.

Your abundance table should be feature-by-sample, with `feature_id` in the first column and one column per sample. The sample IDs must match `sample_metadata.csv`, and taxonomy `feature_id` values must match the abundance table.

## Data Privacy

The shared toy dataset is simulated/desensitized and is provided only to demonstrate the workflow. It does not represent real research results.
