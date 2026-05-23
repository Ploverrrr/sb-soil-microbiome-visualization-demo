# Workflow Notes

This module follows the original differential abundance scripts closely while replacing private input tables with shared simulated toy data.

## Analysis Steps

1. Read shared sample metadata, abundance table, and taxonomy table.
2. Check that required files and columns exist.
3. Check that sample IDs match between metadata and abundance table.
4. Check that feature IDs match between abundance and taxonomy tables.
5. Aggregate feature-level pseudo-counts to the selected taxonomic level. The default is `Genus`.
6. For each contaminated group, subset samples into:

```r
condition = factor(c("control", "Contaminant"),
                   levels = c("control", "Contaminant"))
```

7. Run the original DESeq2-style workflow:

```r
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = col_data,
                              design = ~ condition)
dds <- DESeq(dds, fitType = "parametric")
res <- results(dds, contrast = c("condition", "Contaminant", "control"))
```

8. Remove rows with missing `pvalue`.
9. Classify taxa as `Up`, `Down`, or `Stable` using:

```text
pvalue < 0.05 and abs(log2FoldChange) >= 1
```

10. Select the top up-regulated and down-regulated taxa for volcano labels.
11. Preserve the full `-log10(pvalue)` in the result table, and cap only the plot-specific y value when toy-data p values are extremely small.
12. Draw per-contrast volcano plots with the original visual grammar:
    - `geom_point(aes(color = change), alpha = 1, size = 1.5)`
    - cyan/grey/pink color mapping
    - dashed vertical log2FC thresholds
    - dashed horizontal p-value threshold
    - `theme_bw()`
    - `geom_text_repel()` top labels
13. Combine the three volcano panels with `patchwork`.
14. Select heatmap taxa from the focus comparison, default `Smelting` vs `Control`.
15. Build a row-scaled log2 abundance matrix:

```text
z-score(log2(count + 1)) by taxon
```

16. Draw a rectangular heatmap with `ComplexHeatmap::Heatmap()`.
17. Draw a circular heatmap with `circlize::circos.heatmap()`, preserving the original circular heatmap style.
18. Save PDF and PNG figures to `figures/`.
19. Save all intermediate result tables to `results/`.

## Interpretation

The volcano plots show differential abundance effect size and statistical evidence for each treatment-vs-control contrast. Taxa to the right are enriched in the contaminated group; taxa to the left are enriched in the control group.

The heatmaps show abundance patterns for selected differential taxa across all sample groups. Values are row-scaled log2 counts, so colors represent within-taxon relative enrichment across samples rather than absolute abundance.

## Scope

This is a portfolio visualization demo based on simulated/desensitized toy data. The DESeq2 output and heatmap patterns are designed to demonstrate the workflow and should not be interpreted as real biological results.

For real microbiome count data, users should consider library-size distributions, filtering, compositional effects, design formula choices, and appropriate multiple-testing interpretation.
