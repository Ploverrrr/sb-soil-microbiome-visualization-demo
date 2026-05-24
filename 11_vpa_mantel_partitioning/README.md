# 11 VPA Mantel Partitioning

This independent module demonstrates variation partitioning analysis (VPA) and Mantel-based environmental association visualization for a simulated metal-contaminated soil microbiome study.

The workflow follows the original project scripts:

- build species/community and functional-gene response matrices;
- Hellinger-transform response matrices with `vegan::decostand()`;
- standardize environmental variables;
- run `vegan::varpart()` for `Sb`, `Cu & As`, and `Nutrients`;
- run partial RDA tests for the pure Sb fraction;
- draw the original-style two-panel native `vegan` VPA plot;
- run Mantel tests using Bray-Curtis response distances and Euclidean environmental distances;
- draw a `ggcor::quickcor()` environmental Spearman heatmap with Mantel links.

All data are simulated/desensitized toy data. They demonstrate workflow structure and figure style, not real research results.

## Input Data

Run this module from its own folder. By default, the script reads:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/environmental_variables.csv
../data/toy_shared/abundance_table.csv
../data/toy_shared/functional_annotation_table.csv
```

The module does not create a private copy of shared toy data.

## How To Run

```bash
cd 11_vpa_mantel_partitioning
Rscript scripts/run_demo.R
```

## Outputs

Results:

- `results/species_response_matrix.csv`
- `results/function_response_matrix.csv`
- `results/vpa_fraction_table.csv`
- `results/partial_rda_tests.csv`
- `results/mantel_results.csv`
- `results/environment_spearman_correlation.csv`

Figures:

- `figures/fig1_vpa_combined.pdf/png`
- `figures/fig2_mantel_env_correlation.pdf/png`

## User-Editable Settings

Open `scripts/run_demo.R` and edit the settings block near the top. Common settings include:

- `shared_data_dir`
- `top_n_taxa_for_vpa`
- `top_n_functions_for_vpa`
- `functional_feature_column`
- `functional_value_column`
- `sb_variables`
- `co_metal_variables`
- `nutrient_variables`
- `mantel_environmental_variables`
- `mantel_links_to_plot_per_response`
- plot colors, output file names, and figure sizes

The full Mantel result table is always exported. The figure can show a smaller number of representative links per response block because the simulated toy gradient makes nearly every Mantel test significant; this keeps the portfolio figure readable without hiding the complete computed results.

## Replacing With Your Own Data

To use your own data, provide the same raw-like table structure:

- feature-by-sample abundance table;
- long-format functional annotation abundance table;
- sample metadata table;
- environmental variables table.

The script recalculates response matrices, VPA fractions, partial RDA tests, Mantel tests, and plotting tables from these raw-like inputs. Do not use final VPA or Mantel result tables as the main input if you want the workflow to remain reproducible.
