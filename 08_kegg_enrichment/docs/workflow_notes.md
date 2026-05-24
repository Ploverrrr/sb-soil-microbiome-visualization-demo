# Workflow Notes

This module follows the original KEGG enrichment plotting workflow while making the upstream KO list reproducible from shared toy data.

## Analysis Steps

1. Read shared sample metadata and functional annotation table.
2. Check that required files and columns exist.
3. Check that sample IDs match between metadata and functional annotation.
4. Filter KO-like IDs with `ko_id_pattern`, default `^K`.
5. Aggregate long-format functional counts into a KO-by-sample count matrix.
6. Run DESeq2 for `treatment_group` vs `control_group` to generate a foreground KO list from raw-like toy input.
7. Select foreground KOs using log2 fold-change and p-value thresholds. If the strict threshold returns too few KOs for a small toy enrichment universe, keep the top ranked KOs and mark them as `top_ranked_for_toy_demo`.
8. Run enrichment with one of two backends:

```text
toy_offline:
  clusterProfiler::enricher() with simulated KEGG-like TERM2GENE tables

clusterprofiler_kegg:
  clusterProfiler::enrichKEGG()
  clusterProfiler::enrichMKEGG()
```

9. Export pathway and module enrichment result tables.
10. Draw the original-style combined plots:

```r
ko1 <- barplot(result1)
ko2 <- dotplot(result1)
ko1 + ko2

mo1 <- barplot(result2)
mo2 <- dotplot(result2)
mo1 + mo2
```

11. Draw an auxiliary ggplot bubble plot using `GeneRatio`, `Count`, and `p.adjust`.
12. Save PDF and PNG figures to `figures/`.
13. Save all intermediate result tables to `results/`.

## Interpretation

The pathway and module figures show which KEGG-like terms are over-represented in the foreground KO list compared with the KO universe.

The default offline result is a demonstration of enrichment workflow mechanics. It is not a real KEGG biological result because both the toy KO abundances and the TERM2GENE mapping are simulated.

## Why Use Offline Toy Enrichment By Default?

The original code used `enrichKEGG()` and `enrichMKEGG()` directly. Those functions can require online KEGG access, and the shared toy dataset intentionally contains only a small simulated KO universe. The offline backend keeps the demo reproducible on GitHub while preserving the `clusterProfiler` result object and the original `barplot()` / `dotplot()` plotting workflow.

## Scope

This is a portfolio visualization demo. Users applying the workflow to real data should use valid KO IDs, a complete universe definition, appropriate foreground selection, and a real KEGG annotation source.
