# Workflow Notes

This module adapts the original alpha/beta diversity scripts into a reproducible demo that starts from raw-like feature counts rather than precomputed diversity tables.

## Original Script Logic Preserved

The original alpha diversity script used:

- `vegan` for diversity calculation;
- even-depth rarefaction before alpha index calculation;
- Chao1 and Shannon boxplots;
- `ggpubr::ggboxplot`;
- pairwise significance labels;
- teal/blue/purple group colors.

The original beta diversity script used:

- Bray-Curtis distance;
- `cmdscale()` PCoA;
- `adonis2()` PERMANOVA;
- `metaMDS()` NMDS;
- confidence ellipses;
- half-violin, boxplot, and jitter layers for beta-distance visualization.

This module keeps those methods and visual conventions while replacing private inputs with shared simulated toy data.

## Computational Workflow

1. Read `sample_metadata.csv` and `abundance_table.csv`.
2. Check required columns and sample ID consistency.
3. Convert the abundance table from feature by sample to sample by feature.
4. Optionally treat very low pseudo-counts as absent for alpha diversity, then rarefy counts to the minimum sample depth with a fixed random seed.
5. Calculate alpha diversity indices:
   - observed species;
   - richness;
   - Chao1;
   - ACE;
   - Shannon;
   - Simpson;
   - Pielou evenness;
   - Good's coverage.
6. Run pairwise group tests for configured alpha metrics.
7. Convert community counts to relative abundance for Bray-Curtis beta diversity.
8. Calculate Bray-Curtis distances with `vegan::vegdist()`.
9. Run PCoA with `cmdscale()` and PERMANOVA with `vegan::adonis2()`.
10. Run NMDS with `vegan::metaMDS()`.
11. Build beta-distance tables for pairwise sample distances and sample-group distributions.
12. Save result tables, PDF figures, and PNG figures.

## Figure Design Notes

The alpha diversity figure keeps the original grouped boxplot style and pairwise significance marks. The PCoA and NMDS figures retain the original point-plus-ellipse ordination layout. The beta-distance figure keeps the original half-violin plus boxplot plus jitter composition.

The simulated toy data have a clear but simplified contamination gradient. The statistical results should be interpreted only as workflow examples.
