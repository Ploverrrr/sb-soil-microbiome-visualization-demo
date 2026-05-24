# 08_kegg_enrichment

This module demonstrates a reproducible KEGG-style enrichment workflow for a simulated metal-contaminated soil microbiome study.

The demo follows the original enrichment scripts: a KO foreground list is generated, `clusterProfiler` enrichment is run, `barplot()` and `dotplot()` are drawn, and the two plots are combined with `patchwork`. It does not use a precomputed enrichment result table or manually edited plotting table as input.

## What This Module Shows

- KO count matrix construction from shared functional annotation toy data.
- Differential KO screening for a treatment-vs-control comparison.
- KEGG pathway-style enrichment using `clusterProfiler::enricher()` with an offline toy TERM2GENE mapping.
- Optional original-style `enrichKEGG()` / `enrichMKEGG()` backend for users with KEGG access.
- Original-style combined barplot + dotplot outputs for pathway and module-like enrichment.
- Auxiliary ggplot bubble plot using GeneRatio, Count, and adjusted p-value.

## Shared Toy Inputs

The module reads:

```text
../data/toy_shared/sample_metadata.csv
../data/toy_shared/functional_annotation_table.csv
```

The toy data are simulated/desensitized and are only intended to demonstrate the workflow. They do not represent real study results.

## How To Run

From this module directory:

```bash
cd 08_kegg_enrichment
Rscript scripts/run_demo.R
```

Outputs are written only to:

```text
results/
figures/
```

## Main Outputs

Results:

- `results/ko_count_matrix.csv`
- `results/differential_ko_statistics.csv`
- `results/foreground_ko_list.csv`
- `results/kegg_pathway_enrichment_result.csv`
- `results/kegg_module_enrichment_result.csv`

Figures:

- `figures/kegg_pathway_enrichment_combined.pdf`
- `figures/kegg_pathway_enrichment_combined.png`
- `figures/kegg_module_enrichment_combined.pdf`
- `figures/kegg_module_enrichment_combined.png`
- `figures/kegg_pathway_bubble_plot.pdf`
- `figures/kegg_pathway_bubble_plot.png`

## User-Editable Settings

Edit the settings block at the top of `scripts/run_demo.R` to change:

- `shared_data_dir`
- `control_group`
- `treatment_group`
- `ko_id_pattern`
- `enrichment_backend`
- `foreground_log2fc_threshold`
- `foreground_pvalue_threshold`
- `minimum_foreground_kos`
- `enrichment_pvalue_cutoff`
- `p_adjust_method`
- `qvalue_cutoff`
- `show_category`
- output file names
- figure width and height
- bubble plot colors and sizes

## Backend Choice

The default `enrichment_backend = "toy_offline"` is intentionally reproducible. It uses the simulated KO IDs in the shared toy data and a small KEGG-like TERM2GENE mapping generated inside the script.

`enrichment_backend = "clusterprofiler_kegg"` calls `clusterProfiler::enrichKEGG()` and `clusterProfiler::enrichMKEGG()`, matching the original scripts more directly. This mode may require KEGG online access and can be unstable in offline environments.

## Replacing With Your Own Data

Use the same input structure:

- `sample_metadata.csv` must contain `sample_id` and `group`.
- `functional_annotation_table.csv` must contain `sample_id`, `ko_id`, `pathway`, and `count`.

For real KEGG analysis, use valid KEGG Orthology IDs such as `K02111`. You can keep the offline backend with a project-specific TERM2GENE mapping, or switch to `clusterprofiler_kegg` if online KEGG annotation is available.

## Original-Script Features Preserved

The original scripts used `clusterProfiler`, `enrichKEGG()`, `enrichMKEGG()`, `barplot()`, `dotplot()`, `patchwork`, and `ggsave()` to produce combined KO and module enrichment figures. This module preserves that plotting workflow and adds reproducible upstream KO selection from shared toy data.
