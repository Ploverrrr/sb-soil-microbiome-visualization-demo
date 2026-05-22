# Soil Microbiome Visualization Portfolio

This repository is a GitHub portfolio project for environmental microbiome and metagenomic data visualization. It is organized around publication-style figure types inspired by an original research workflow on soil microbiomes, environmental factors, contaminant gradients, taxonomic profiles, functional annotations, and mechanism-oriented statistical graphics.

Each numbered folder is an independent, reproducible demo module.

All public data used in this repository will be simulated, toy, or desensitized data. Public demo modules must not depend on private raw data, manually edited final plotting tables, or intermediate result files copied from the original research project.

The `_private_original/` directory is a git-ignored private reference folder. It is not part of the public project and must not be copied, tracked, or used as an input source for public demo modules.

Each module starts from raw-like toy input data, performs the required calculations inside the module, and writes reproducible outputs to its own `results/` and `figures/` folders.

## Module Structure

Each numbered module follows the same layout:

```text
XX_module_name/
├── README.md
├── docs/
│   ├── data_schema.md
│   └── workflow_notes.md
├── data/
│   └── toy/
├── scripts/
├── results/
└── figures/
```

## Confirmed Demo Modules

1. `01_rf_correlation_heatmap`
2. `02_microbe_env_network`
3. `03_ternary_taxa_distribution`
4. `04_faprotax_functional_profile`
5. `05_lefse_biomarker`
6. `06_plspm_mechanism_model`
7. `07_differential_volcano_heatmap`
8. `08_kegg_enrichment`
9. `09_sulfur_gene_contaminant_association`
10. `10_alpha_beta_diversity`
11. `11_vpa_mantel_partitioning`
