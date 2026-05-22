# 02 Microbe Env Network

This module demonstrates a reproducible microbe-environment association network for a simulated metal-contaminated soil microbiome study.

The workflow reads the repository-level shared toy dataset, calculates microbial relative abundance, aggregates microbial features to a selected taxonomic level, computes microbe-environment associations, and builds a bipartite network linking microbial taxa to environmental variables.

## What This Module Shows

The demo starts from raw-like toy tables:

- sample metadata;
- environmental variables;
- feature-by-sample abundance;
- taxonomy annotations.

It does not use precomputed edge lists, Gephi files, or private original data.

The output includes:

- a network edge list;
- a network node table;
- a correlation matrix;
- a publication-style network plot.

## Shared Toy Inputs

Default inputs:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/environmental_variables.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/taxonomy_table.csv
```

The module does not duplicate shared toy data in `data/toy/`.

## How To Run

From the repository root:

```bash
cd 02_microbe_env_network
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## User-Editable Parameters

Edit the settings block near the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `target_taxonomic_level`
- `top_n_taxa`
- `environmental_variables`
- `correlation_method`
- `p_adjust_method`
- `correlation_threshold`
- `adjusted_p_threshold`
- output file names
- figure width/height
- node and edge color settings

## Interpreting The Network

Edges represent statistically filtered associations between mean-normalized microbial taxon abundance and environmental variables. Edge color indicates positive or negative correlation. Edge width scales with correlation strength.

This network is an exploratory visualization of associations. It should not be interpreted as causal evidence.

## Using Your Own Data

To use your own data, keep the same shared table schema or update the file paths and settings in `scripts/run_demo.R`.

Required relationships:

- sample IDs must match across metadata, environmental variables, and abundance columns;
- feature IDs must match between abundance and taxonomy tables;
- selected environmental variables must be numeric.

The shared toy dataset is simulated/desensitized and does not represent real research results.
