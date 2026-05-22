# Workflow Notes

This module implements a simplified LEfSe-style biomarker workflow from raw-like toy data.

## Analysis Steps

1. Read shared sample metadata, abundance table, and taxonomy table.
2. Check that required files and columns exist.
3. Check that sample IDs match between metadata and abundance table.
4. Check that feature IDs match between abundance and taxonomy tables.
5. Convert feature counts to within-sample relative abundance.
6. Aggregate feature-level abundance to the selected taxonomy rank.
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
16. Draw a group-colored biomarker barplot and a group mean abundance heatmap.
17. Save PDF and PNG figures to `figures/`.

## Interpretation

The barplot ranks taxa by a transparent effect-size score. The color indicates the group where the taxon has the highest mean relative abundance.

The heatmap shows the same selected biomarkers across all groups, making it easier to see whether the detected signal follows the simulated contamination gradient.

## Scope

This module is intentionally reproducible and dependency-light. It does not call the original LEfSe software and does not require a precomputed LDA table. For formal biological analysis, users should run the original LEfSe or an equivalent validated biomarker method and adapt this visualization workflow carefully.
