# 04 FAPROTAX Functional Profile

This module demonstrates a reproducible FAPROTAX-style functional profile workflow for a simulated metal-contaminated soil microbiome study.

The demo uses the repository-level shared toy dataset. It does not use real functional annotation results, private data, or manually prepared plotting tables.

## What This Module Shows

The workflow starts from a long-format functional annotation table and sample metadata. It:

1. filters selected functional annotation sources;
2. aggregates functional counts by sample and pathway;
3. calculates both count-weighted relative abundance and unweighted feature-assignment percentage, matching the original FAPROTAX `abundance_weighted = FALSE` idea;
4. summarizes functions by group;
5. selects top functions;
6. draws an unweighted functional bubble profile, a z-score bubble profile, and a grouped barplot.

The visual style is inspired by FAPROTAX functional profile figures, but the data are simulated/desensitized and not real FAPROTAX outputs.

## Shared Toy Inputs

Default inputs:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/functional_annotation_table.csv
```

The module also accepts the shared abundance and taxonomy tables in the project, but this demo only needs the functional annotation long table and metadata.

The module does not duplicate shared toy data in `data/toy/`.

## How To Run

From the repository root:

```bash
cd 04_faprotax_functional_profile
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## User-Editable Parameters

Edit the settings block in `scripts/run_demo.R` to change:

- `shared_data_dir`
- `annotation_sources_to_include`
- `function_group_column`
- `top_n_functions`
- group order and colors
- output file names
- figure sizes

## Outputs

This module writes:

- sample-level functional relative abundance table;
- sample-level unweighted feature-assignment percentage table;
- group-level mean/sd/standard error table;
- plot-ready bubble, z-score bubble, and barplot tables;
- functional bubble profile PDF/PNG;
- scaled count-weighted FAPROTAX bubble profile PDF/PNG;
- grouped functional barplot PDF/PNG.

The scaled bubble plot is included because the original FAPROTAX reference script scaled function profiles before plotting point size and color. In this toy workflow, the scaled bubble plot uses the count-weighted functional profile so the simulated gradient remains visible.

## Using Your Own Data

To use your own data, preserve the same long-format functional annotation schema or update the path and column settings in `scripts/run_demo.R`.

At minimum, your functional annotation table should include:

- `sample_id`
- `annotation_source`
- `function_category`
- `pathway`
- `count`

The shared toy data are simulated/desensitized and do not represent real research results.
