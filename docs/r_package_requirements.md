# R Package Requirements

This project is an R-based public portfolio demo for environmental microbiome and contaminant-gradient visualization. The modules use simulated/desensitized toy data and write reproducible results and figures from their own `scripts/run_demo.R` files.

The module scripts do not install packages automatically. They check for required packages and stop with a clear message if something is missing.

For a fresh public-style R environment, use the repository-level helper from the project root:

```bash
Rscript scripts/install_r_dependencies.R
```

Then regenerate toy data and run the full demo set:

```bash
Rscript scripts/create_shared_toy_data.R
Rscript scripts/run_all_demos.R
```

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

## Automated Installation Helper

The recommended installation path is:

```bash
Rscript scripts/install_r_dependencies.R
```

The helper installs the CRAN package group, the Bioconductor package group, and `ggcor`. It installs only required dependency types (`Depends`, `Imports`, and `LinkingTo`) rather than optional `Suggests`, which keeps public CI focused on packages actually used by the demo scripts.

The helper first tries to install `gghalves`, `aplot`, and `ggcor` from the configured CRAN repository. If `gghalves` or `aplot` are unavailable from that repository, it can try GitHub fallbacks. If `ggcor` is unavailable, it can try `remotes::install_github("houyunhuang/ggcor")`.

Useful environment-variable options:

```bash
CRAN_REPO=https://packagemanager.posit.co/cran/latest Rscript scripts/install_r_dependencies.R
NCPUS=4 Rscript scripts/install_r_dependencies.R
INSTALL_BIOC_PACKAGES=false Rscript scripts/install_r_dependencies.R
INSTALL_CRAN_GITHUB_FALLBACKS=false Rscript scripts/install_r_dependencies.R
INSTALL_GGCOR_FROM_GITHUB=false Rscript scripts/install_r_dependencies.R
```

Use `INSTALL_BIOC_PACKAGES=false` only for partial testing, because modules `07_differential_volcano_heatmap` and `08_kegg_enrichment` require Bioconductor packages.

## Manual Installation Pattern

If you prefer to install manually, install CRAN packages from an R session:


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

Install `ggcor` if it is available from your configured repository:

```r
install.packages("ggcor")
```

If `ggcor` is unavailable from your CRAN mirror, install it from its upstream source according to the package maintainer's current instructions. The automated helper uses `remotes::install_github("houyunhuang/ggcor")` as its fallback because the Mantel module follows the original `ggcor` plotting workflow.

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
