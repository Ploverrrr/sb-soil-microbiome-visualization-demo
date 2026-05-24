# 10 Alpha Beta Diversity

This independent module demonstrates alpha diversity, beta diversity, and ordination analysis for a simulated metal-contaminated soil microbiome study.

The workflow follows the original project scripts for alpha/beta diversity figures:

- optional even-depth rarefaction before alpha diversity calculation;
- `vegan` alpha diversity indices including Chao1 and Shannon;
- grouped `ggpubr` boxplots with pairwise significance labels;
- Bray-Curtis community dissimilarity;
- PCoA with confidence ellipses and PERMANOVA caption;
- NMDS with stress value;
- beta-distance violin/box/jitter plot using the original half-violin style.

All data used here are simulated/desensitized toy data. They are designed to show the workflow and plotting style, not to represent real research results.

## Input Data

Run the module from this folder. By default, the script reads shared toy data from:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/abundance_table.csv
```

This module does not create or store its own copy of the shared microbiome tables.

## How To Run

```bash
cd 10_alpha_beta_diversity
Rscript scripts/run_demo.R
```

The script writes intermediate tables to `results/` and figures to `figures/`.

## Main Outputs

Results:

- `results/alpha_diversity_indices.csv`
- `results/alpha_diversity_pairwise_tests.csv`
- `results/pcoa_coordinates.csv`
- `results/nmds_coordinates.csv`
- `results/permanova_results.csv`
- `results/beta_pairwise_distances.csv`
- `results/beta_distance_by_sample_group.csv`

Figures:

- `figures/fig1_alpha_diversity_boxplots.pdf/png`
- `figures/fig2_pcoa_bray_curtis.pdf/png`
- `figures/fig3_nmds_bray_curtis.pdf/png`
- `figures/fig4_beta_distance_boxplot.pdf/png`

## User-Editable Settings

Open `scripts/run_demo.R` and edit the settings block near the top. Common settings include:

- `shared_data_dir`
- `group_column`
- `group_order`
- `group_palette`
- `rarefy_counts`
- `alpha_min_count_threshold`
- `alpha_metrics_to_plot`
- `beta_distance_method`
- `use_relative_abundance_for_beta`
- output file names and figure sizes

The default palette keeps the teal-blue-purple style used in the original alpha diversity figure while extending it to the four simulated groups.

Because the shared toy abundance table uses pseudo-counts, `alpha_min_count_threshold` treats very low counts as absent before alpha diversity calculation. This keeps richness-style metrics interpretable without changing the shared input files.

## Replacing With Your Own Data

To use your own data, replace the shared input files or point `shared_data_dir` to another folder with the same schema:

- `sample_metadata.csv` must include `sample_id` and the configured grouping column.
- `abundance_table.csv` must contain one `feature_id` column followed by numeric sample columns.

The script recalculates alpha indices, Bray-Curtis distances, ordination coordinates, PERMANOVA, NMDS, and plotting tables from these raw-like inputs. Do not use a final plotting table as the main input if you want the workflow to remain reproducible.
