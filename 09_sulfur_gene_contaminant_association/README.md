# 09_sulfur_gene_contaminant_association

This module demonstrates a reproducible sulfur gene / contaminant association workflow for a simulated metal-contaminated soil microbiome study.

The workflow follows the original sulfur gene script: target sulfur KO extraction, Pearson correlation analysis, `ggcorrplot` heatmap, significant scatter regression plots, and multiple linear regression summary figures. It does not use a precomputed correlation table or manually edited plotting table as input.

## What This Module Shows

- Extraction of sulfur cycling KO abundances from long-format functional annotation data.
- Construction of dsr/sox-style combined indicators from available target KOs.
- Pearson correlation testing between sulfur gene indicators and contaminants.
- Original-style Pearson heatmap using `ggcorrplot()`.
- Original-style scatter regression panels for the strongest sulfur gene / contaminant associations.
- MLR fitted-vs-measured and standardized coefficient plots for Sb species.

## Shared Toy Inputs

The module reads:

```text
../data/toy_shared/environmental_variables.csv
../data/toy_shared/functional_annotation_table.csv
```

The toy data are simulated/desensitized and are only intended to demonstrate the workflow. They do not represent real study results.

## How To Run

From this module directory:

```bash
cd 09_sulfur_gene_contaminant_association
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## Main Outputs

Results:

- `results/sulfur_gene_abundance_by_sample.csv`
- `results/sulfur_gene_contaminant_analysis_table.csv`
- `results/pearson_results.csv`
- `results/mlr_results.csv`

Figures:

- `figures/fig1_pearson_heatmap.pdf`
- `figures/fig1_pearson_heatmap.png`
- `figures/fig2_pearson_scatter.pdf`
- `figures/fig2_pearson_scatter.png`
- `figures/fig3_mlr_results.pdf`
- `figures/fig3_mlr_results.png`

## User-Editable Settings

Edit the settings block at the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `target_ko_names`
- `functional_value_column`
- `contaminant_vars`
- `primary_response_vars`
- `preferred_gene_order`
- `scatter_gene_candidates`
- `max_scatter_pairs`
- output file names
- figure width and height
- heatmap, scatter, and regression colors

## Replacing With Your Own Data

Use the same input structure:

- `environmental_variables.csv` must contain `sample_id` and the selected contaminant variables.
- `functional_annotation_table.csv` must contain `sample_id`, `ko_id`, and the selected abundance/count value column.

For real sulfur gene work, update `target_ko_names` to match the KO IDs present in your annotation table. The original script targeted `K11180`, `K11181`, `K17218`, `K17222`, and `K17223` for dsr/sox genes.

## Toy Data Limitation

The shared toy dataset includes simulated sulfur-related KOs such as `K11180`, `K17218`, and `K00958`, but it does not include a full dsrAB/soxBCD gene set. The script keeps the original combined-indicator logic and builds `dsrAB` or `soxBCD` from whichever target genes are available. This is appropriate for demonstrating the workflow, but users should provide complete gene annotations for real analysis.

## Original-Script Features Preserved

The original sulfur script produced `Pearson_results.csv`, a Pearson heatmap, scatter regression plots for selected significant pairs, and MLR result plots. This module preserves that analysis and plotting structure while replacing private inputs and absolute paths with shared simulated toy data and relative paths.
