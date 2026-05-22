# Demo Selection

The portfolio uses numbered independent modules. Each module represents one publication-style figure type or closely related figure family.

## Confirmed Modules

| Module | Figure Capability |
|---|---|
| `01_rf_correlation_heatmap` | Random forest importance plus environmental correlation heatmap |
| `02_microbe_env_network` | Microbe-environment association network |
| `03_ternary_taxa_distribution` | Ternary plot for taxa distribution across three groups |
| `04_faprotax_functional_profile` | FAPROTAX-style functional profile and group summary |
| `05_lefse_biomarker` | LEfSe biomarker visualization |
| `06_plspm_mechanism_model` | PLS-PM mechanism model |
| `07_differential_volcano_heatmap` | Differential abundance, volcano plot, and heatmap |
| `08_kegg_enrichment` | KEGG enrichment bubble and bar plots |
| `09_sulfur_gene_contaminant_association` | Sulfur gene and contaminant association visualization |
| `10_alpha_beta_diversity` | Alpha and beta diversity with ordination |
| `11_vpa_mantel_partitioning` | VPA and Mantel environmental partitioning |

## Priority Logic

The first build phase should prioritize the required paper-style figure types:

1. `01_rf_correlation_heatmap`
2. `02_microbe_env_network`
3. `03_ternary_taxa_distribution`
4. `04_faprotax_functional_profile`
5. `05_lefse_biomarker`
6. `06_plspm_mechanism_model`

The second build phase should complete the complementary modules:

1. `07_differential_volcano_heatmap`
2. `08_kegg_enrichment`
3. `09_sulfur_gene_contaminant_association`
4. `10_alpha_beta_diversity`
5. `11_vpa_mantel_partitioning`
