# Data Schema

This module reads shared toy data from `../data/toy_shared/`. The data are simulated/desensitized and are designed to mimic a metal-contaminated soil microbiome study.

## `sample_metadata.csv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Unique sample identifier. Must match sample columns in `abundance_table.csv`. |
| `group` | Group used for DESeq2 contrasts. Defaults are `Control`, `Tailing`, `Mining`, and `Smelting`. |

Optional columns such as `replicate`, `site_type`, and `batch` can be present. The default script uses `replicate` only for sample ordering when it is available.

## `abundance_table.csv`

Required structure:

| Column | Description |
| --- | --- |
| `feature_id` | Unique microbial feature or ASV/OTU identifier. |
| sample columns | One column per `sample_id`, containing raw counts or pseudo-count abundance. |

The script aggregates feature counts to the selected `target_taxonomic_level` before DESeq2 analysis. Counts are rounded to integer values for `DESeqDataSetFromMatrix()`.

## `taxonomy_table.csv`

Required columns:

| Column | Description |
| --- | --- |
| `feature_id` | Feature identifier matching `abundance_table.csv`. |
| `Genus` | Default taxonomic level used for differential abundance testing. |

Other taxonomy levels such as `Phylum`, `Class`, `Order`, or `Family` can be used by changing `target_taxonomic_level` in `scripts/run_demo.R`.

## Output Tables

| Output | Description |
| --- | --- |
| `taxon_count_matrix.csv` | Taxon-by-sample count matrix used as DESeq2 input. |
| `deseq2_all_contrasts.csv` | Combined DESeq2 output for all treatment-vs-control comparisons. |
| `deseq2_tailing_vs_control.csv` | DESeq2 result table for `Tailing` vs `Control`. |
| `deseq2_mining_vs_control.csv` | DESeq2 result table for `Mining` vs `Control`. |
| `deseq2_smelting_vs_control.csv` | DESeq2 result table for `Smelting` vs `Control`. |
| `heatmap_zscore_matrix.csv` | Row-scaled log2 abundance matrix used for heatmap plotting. |
| `heatmap_selected_taxa.csv` | Differential taxa selected for heatmap visualization. |

The DESeq2 result tables include `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`, `taxon`, `change`, `comparison`, `gene`, `negative_log10_pvalue`, and `negative_log10_pvalue_for_plot`. The full p-value transform is preserved; the plot-specific column is capped only to keep toy-data volcano labels readable.
