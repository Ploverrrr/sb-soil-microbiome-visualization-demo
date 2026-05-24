# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are designed to mimic a metal-contaminated soil microbiome study.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match `sample_id` values in the functional annotation table. |
| `group` | Group used for treatment-vs-control KO screening. Defaults are `Control`, `Tailing`, `Mining`, and `Smelting`. |

## `functional_annotation_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Sample identifier. |
| `ko_id` | KO-like identifier. Default filtering keeps IDs matching `^K`. |
| `pathway` | Pathway or function label. Preserved for interpretation. |
| `count` | Simulated functional pseudo-count. |

Optional columns such as `gene_id`, `feature_id`, `annotation_source`, and `function_category` can be present and are preserved in the shared dataset.

## Derived Tables

| Output | Description |
| --- | --- |
| `ko_count_matrix.csv` | KO-by-sample count matrix generated from the long functional annotation table. |
| `differential_ko_statistics.csv` | DESeq2 statistics for `treatment_group` vs `control_group`. |
| `foreground_ko_list.csv` | KO IDs selected for enrichment, including the selection reason. |
| `kegg_pathway_enrichment_result.csv` | KEGG pathway-style enrichment table from `clusterProfiler`. |
| `kegg_module_enrichment_result.csv` | KEGG module-like enrichment table from `clusterProfiler`. |

The enrichment result tables include the standard `clusterProfiler` columns such as `ID`, `Description`, `GeneRatio`, `BgRatio`, `pvalue`, `p.adjust`, `qvalue`, `geneID`, and `Count`.

## Toy TERM2GENE Mapping

The default offline demo generates a small KEGG-like TERM2GENE mapping inside `scripts/utils.R`. This is necessary because the shared toy dataset contains only a small set of simulated KO IDs. The mapping is designed only to demonstrate the enrichment workflow and should not be interpreted as a real KEGG database.

Users can switch `enrichment_backend` to `clusterprofiler_kegg` for real KEGG calls, or replace the offline TERM2GENE mapping with their own curated project mapping.
