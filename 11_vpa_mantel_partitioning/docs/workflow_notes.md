# Workflow Notes

This module adapts the original VPA and Mantel scripts into a reproducible public demo that starts from raw-like toy tables.

## Original VPA Logic Preserved

The original VPA script:

1. read species, functional gene, and environmental matrices;
2. transposed and aligned samples;
3. applied Hellinger transformation to species and gene matrices;
4. standardized environmental predictors;
5. grouped predictors into `Sb`, `Cu & As`, and `Nutrients`;
6. ran `vegan::varpart()` separately for species community and functional genes;
7. plotted two native `vegan` VPA diagrams side by side;
8. ran partial RDA to test the pure Sb fraction.

This module keeps that same analysis structure and plotting style, while replacing absolute paths and private input tables with shared simulated toy data.

## Original Mantel Logic Preserved

The original Mantel script:

1. read a response matrix and an environmental matrix;
2. used `ggcor::mantel_test()` with Bray-Curtis response distance and Euclidean environmental distance;
3. binned Mantel `r` and `p` values into legend categories;
4. drew an environmental Spearman correlation heatmap using `ggcor::quickcor()`;
5. overlaid Mantel links with `anno_link()`;
6. used a green Spearman heatmap palette and link colors based on Mantel significance.

This module preserves that visual grammar and computes two response blocks using the reference workflow's `spec.select` structure: `callvulg = 1:2` and `B = 3:44`. The plot keeps a bare `geom_square()` call in the plotting block; the script defines a tiny compatibility wrapper beforehand so that call works with the installed `ggcor`/`ggplot2` versions.

## Computational Workflow

1. Read shared toy metadata, environmental variables, taxonomic abundance, and functional annotation data.
2. Build a sample-by-feature species/community matrix from top taxonomic features.
3. Build a sample-by-function matrix from long-format functional annotations.
4. Save both response matrices to `results/`.
5. Apply Hellinger transformation to both response matrices.
6. Scale environmental variables.
7. Run `vegan::varpart()` for species community and functional genes.
8. Save cleaned and raw adjusted R2 fractions.
9. Run partial RDA tests for the pure Sb fraction while conditioning on co-metals and nutrients.
10. Draw the native two-panel VPA figure.
11. Run Mantel tests for the reference `callvulg` and `B` response blocks against selected environmental variables.
12. Draw the Mantel-link environmental correlation heatmap from the complete Mantel table.

## Interpretation Notes

The shared toy data contain a simplified contamination gradient. VPA fractions, Mantel statistics, and significance values are generated only to demonstrate a reproducible workflow and figure style. They should not be interpreted as real ecological evidence.
