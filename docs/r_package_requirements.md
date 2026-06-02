# R Package Requirements

This project is an R-based public portfolio demo for environmental microbiome and contaminant-gradient visualization. The modules use simulated/desensitized toy data and write reproducible results and figures from their own `scripts/run_demo.R` files.

The repository does not automatically install packages. Install requirements in your local R environment before running the relevant modules.

## Core Package Groups

### CRAN Packages

The module scripts currently require or use the following CRAN packages:

- `ggplot2`
- `dplyr`
- `tidyr`
- `tibble`
- `patchwork`
- `ggrepel`
- `randomForest`
- `igraph`
- `ggtern`
- `ggpubr`
- `gghalves`
- `ggsignif`
- `vegan`
- `ggcorrplot`
- `plspm`
- `microeco`
- `aplot`
- `gridExtra`
- `circlize`

`ggcor` is also required by `11_vpa_mantel_partitioning`. Its installation source may vary by local R setup and should be verified before publishing installation instructions for a specific environment.

### Bioconductor Packages

The differential and enrichment modules require Bioconductor packages:

- `DESeq2`
- `ComplexHeatmap`
- `clusterProfiler`
- `enrichplot`

## Module-Specific Dependencies

| Module | Packages required by `scripts/run_demo.R` |
|---|---|
| `01_rf_correlation_heatmap` | `ggplot2`, `randomForest`, `patchwork` |
| `02_microbe_env_network` | `igraph` |
| `03_ternary_taxa_distribution` | `ggtern` |
| `04_faprotax_functional_profile` | `ggplot2` |
| `05_lefse_biomarker` | `ggplot2`, `microeco`, `aplot`, `gridExtra` |
| `06_plspm_mechanism_model` | `plspm` |
| `07_differential_volcano_heatmap` | `DESeq2`, `ggplot2`, `ggrepel`, `patchwork`, `ComplexHeatmap`, `circlize` |
| `08_kegg_enrichment` | `DESeq2`, `clusterProfiler`, `enrichplot`, `ggplot2`, `patchwork` |
| `09_sulfur_gene_contaminant_association` | `dplyr`, `tidyr`, `tibble`, `ggplot2`, `ggcorrplot`, `ggrepel`, `patchwork` |
| `10_alpha_beta_diversity` | `vegan`, `ggplot2`, `ggpubr`, `patchwork`, `gghalves`, `ggsignif` |
| `11_vpa_mantel_partitioning` | `vegan`, `ggplot2`, `ggcor`, `dplyr` |

## Suggested Installation Pattern

Install CRAN packages from an R session:

```r
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "patchwork",
  "ggrepel",
  "randomForest",
  "igraph",
  "ggtern",
  "ggpubr",
  "gghalves",
  "ggsignif",
  "vegan",
  "ggcorrplot",
  "plspm",
  "microeco",
  "aplot",
  "gridExtra",
  "circlize"
))
```

Install Bioconductor packages from an R session:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "DESeq2",
  "ComplexHeatmap",
  "clusterProfiler",
  "enrichplot"
))
```

Before running all modules, verify any locally uncommon packages such as `ggcor` and confirm package versions if strict reproducibility is required.

## Running After Installation

Generate shared toy data from the repository root:

```bash
Rscript scripts/create_shared_toy_data.R
```

Run a module from its own folder:

```bash
cd 01_rf_correlation_heatmap
Rscript scripts/run_demo.R
```

Or run the optional root-level batch helper:

```bash
Rscript scripts/run_all_demos.R
```
