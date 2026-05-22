# Workflow Notes

This module implements a simplified LEfSe-style biomarker workflow from raw-like toy data.

## Analysis Steps

1. Read shared sample metadata, abundance table, and taxonomy table.
2. Check that required files and columns exist.
3. Check that sample IDs match between metadata and abundance table.
4. Check that feature IDs match between abundance and taxonomy tables.
5. Convert feature counts to within-sample relative abundance.
6. Aggregate feature-level abundance to all configured taxonomic ranks by default, matching the original LEfSe `taxa_level = "all"` workflow.
7. Calculate prevalence and overall mean abundance for each taxon.
8. For each taxon, run a Kruskal-Wallis test across groups.
9. Adjust Kruskal-Wallis p values with the selected FDR method.
10. Assign each taxon to the group with the highest mean abundance.
11. Run a one-vs-rest Wilcoxon test for that enriched group.
12. Calculate a LEfSe-style log-ratio effect score:

```text
log10((mean abundance in enriched group + pseudocount) /
      (mean abundance in other groups + pseudocount))
```

13. Filter biomarkers by prevalence, mean abundance, FDR, and effect score.
14. If strict filters return too few biomarkers for a useful demo figure, keep the strongest ranked taxa and mark them as relaxed demo selections.
15. Write intermediate result tables to `results/`.
16. Build a LEfSe-style concentric-ring cladogram from the taxonomy hierarchy.
17. Draw a group-colored biomarker barplot, cladogram, KW abundance plot, combined LDA+cladogram figure, combined LDA+KW figure, and auxiliary group mean heatmap.
18. Save PDF and PNG figures to `figures/`.

## Interpretation

The barplot ranks taxa by a transparent effect-size score. The color indicates the group where the taxon has the highest mean relative abundance.

The heatmap shows the same selected biomarkers across all groups, making it easier to see whether the detected signal follows the simulated contamination gradient.

The cladogram follows the original figure intent: nodes are arranged by taxonomic hierarchy, rings represent taxonomic depth from Phylum to Genus, and colored nodes mark taxa enriched in different groups. It is a ggplot-based implementation rather than a direct `microeco::plot_diff_cladogram()` call, but its visual grammar is intentionally aligned with the original LEfSe cladogram.

The KW abundance plot uses the same selected biomarkers as the LDA-style barplot and displays group mean relative abundance with standard error, paralleling the original `plot_diff_abund()` figure.

## Scope

This module is intentionally reproducible and dependency-light. It does not call the original LEfSe software and does not require a precomputed LDA table. For formal biological analysis, users should run the original LEfSe or an equivalent validated biomarker method and adapt this visualization workflow carefully.
