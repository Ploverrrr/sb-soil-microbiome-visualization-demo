# Workflow Notes

This module follows the original sulfur gene / Sb script while replacing private inputs with shared simulated toy data.

## Analysis Steps

1. Read shared environmental variables and functional annotation tables.
2. Check that required files and columns exist.
3. Check that sample IDs match between the two input tables.
4. Extract target sulfur KO IDs from `functional_annotation_table.csv`.
5. Aggregate the selected value column, default `abundance`, by KO and sample.
6. Rename target KOs to gene labels such as `dsrA`, `soxB`, and `sat`.
7. Create combined indicators following the original logic:

```text
dsrAB  = mean(dsrA, dsrB), using available columns
soxBCD = mean(soxB, soxC, soxD), using available columns
```

8. Merge sulfur gene indicators with selected contaminant variables.
9. Run pairwise Pearson correlations between sulfur indicators and contaminants.
10. Export the Pearson result table.
11. Draw the original-style Pearson heatmap with:

```r
ggcorrplot(
  cor_mat,
  method = "square",
  type = "full",
  lab = TRUE,
  p.mat = cor_pmat,
  sig.level = 0.05,
  insig = "blank",
  colors = c("#A9D1E8", "white", "#E3A8A8")
)
```

12. Select significant gene-contaminant pairs for scatter plots. The default scatter panels prioritize combined sulfur indicators (`dsrAB`, `soxBCD`, `sat`) and select diverse gene indicators before filling remaining panels, avoiding repetitive scatter plots when one toy variable dominates the correlations.
13. Draw scatter plots with `geom_point()`, `geom_smooth(method = "lm")`, and Pearson `r` labels.
14. Fit multiple linear regression models for selected Sb species using sulfur indicators as predictors.
15. Export MLR coefficient and model summary tables.
16. Draw fitted-vs-measured plots and standardized coefficient plots.
17. Save PDF and PNG figures to `figures/`.

## Interpretation

The Pearson heatmap shows pairwise associations among sulfur gene indicators and contaminant variables. Blank cells indicate non-significant correlations at the default p-value threshold.

The scatter plots show the strongest individual sulfur gene / contaminant relationships with linear regression ribbons.

The MLR plots summarize whether sulfur indicators jointly explain selected Sb species in the simulated toy dataset. Coefficients are standardized for comparability.

## Scope

The shared toy dataset is simulated/desensitized and does not contain a complete real sulfur gene catalog. This module demonstrates a reproducible workflow and original figure style. For real analysis, users should supply a complete KO abundance table and carefully review gene selection, normalization, and model assumptions.
