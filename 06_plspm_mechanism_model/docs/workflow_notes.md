# Workflow Notes

This module follows the original PLS-PM scripts closely: build a scaled manifest-variable table, define latent-variable blocks, define a 0/1 path matrix, run `plspm()`, inspect model outputs, and draw native `innerplot()` and `outerplot()` figures.

## Analysis Steps

1. Read shared sample metadata, environmental variables, abundance table, and functional annotation table.
2. Check that required files and columns exist.
3. Check that sample IDs match across all input tables.
4. Derive alpha-diversity indicators from `abundance_table.csv`. Because the shared toy abundance table uses positive pseudo-counts for all features, the Chao1-style richness indicator is calculated after deterministic pseudo-rarefaction to avoid a constant richness value:
   - Chao1
   - Shannon
   - Simpson
   - Pielou
5. Derive functional indicators from `functional_annotation_table.csv`:
   - metal resistance
   - sulfur cycling
   - nitrogen cycling
   - carbon cycling
6. Combine environmental indicators, derived functional indicators, and alpha-diversity indicators into one model input table.
7. Scale manifest variables with `scale()`, matching the original script pattern.
8. Define default latent-variable blocks:

```r
Environment = c("pH", "EC")
Nutrient    = c("TOC", "TN", "TP")
Sb_As       = c("Sb_total", "As")
Function    = c("MetalRes", "Sulfur", "Nitrogen", "Carbon")
Alpha       = c("Chao1", "Shannon", "Simpson", "Pielou")
```

9. Define the default path matrix:

```text
Environment -> Nutrient
Environment -> Sb_As
Nutrient    -> Sb_As
Environment -> Function
Nutrient    -> Function
Sb_As       -> Function
Environment -> Alpha
Nutrient    -> Alpha
Sb_As       -> Alpha
Function    -> Alpha
```

10. Run `plspm::plspm()` with reflective mode `"A"` for all blocks.
11. Export path coefficients, inner model statistics, outer model values, effects, latent scores, and model metrics.
12. Draw the native inner path diagram with:

```r
innerplot(dat_pls, colpos = "red", colneg = "blue",
          show.values = TRUE, lcol = "gray", box.lwd = 0)
```

13. Draw native measurement-model diagrams with:

```r
outerplot(dat_pls, what = "loadings", arr.width = 0.1,
          colpos = "red", colneg = "blue",
          show.values = TRUE, lcol = "gray")

outerplot(dat_pls, what = "weights", arr.width = 0.1,
          colpos = "red", colneg = "blue",
          show.values = TRUE, lcol = "gray")
```

14. Save PDF and PNG figures to `figures/`.

## Interpretation

The inner path model summarizes hypothesized mechanism links among environmental conditions, nutrient status, Sb/As contamination, microbial functional potential, and alpha diversity.

Positive paths are colored red and negative paths are colored blue, preserving the original `plspm` plotting style. Displayed path values are model-estimated coefficients from the simulated toy dataset.

The outer loading and weight plots show how each manifest variable contributes to its latent variable. These are measurement-model diagnostics rather than final ecological conclusions.

The total-effects plot is an auxiliary summary derived from `dat_pls$effects`. The main preserved original figures are the native `innerplot()` and `outerplot()` outputs.

## Scope

This is a model-structure and visualization demo. The toy data are simulated/desensitized, and the resulting coefficients should not be interpreted as real causal evidence. Users replacing the toy data should evaluate model assumptions, block definitions, sample size, indicator scaling, and path structure for their own study.
