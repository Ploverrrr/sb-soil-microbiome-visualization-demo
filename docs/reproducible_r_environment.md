# Reproducible R Environment

This repository is designed so a public user can recreate all demo outputs from simulated toy inputs. The expected full run is:

```bash
Rscript scripts/install_r_dependencies.R
Rscript scripts/create_shared_toy_data.R
Rscript scripts/run_all_demos.R
```

The module scripts do not install packages automatically. They check for required packages and stop with a clear message if something is missing. Package installation is handled by `scripts/install_r_dependencies.R` so environment setup remains explicit and reproducible.

## Fresh Environment Test

For a local smoke test without mixing packages into your usual R library, you can use a temporary user library:

```bash
mkdir -p /tmp/sb_demo_r_lib
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/install_r_dependencies.R
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/create_shared_toy_data.R
R_LIBS_USER=/tmp/sb_demo_r_lib Rscript scripts/run_all_demos.R
```

On macOS or Linux, this is a practical way to check whether the repository can run in a clean public-style R environment.

## Lightweight Smoke Test

For a faster check that avoids the Bioconductor-heavy modules, run:

```bash
Rscript scripts/run_smoke_test.R
```

The smoke test regenerates the shared toy data, runs modules `01` through `04`, checks key CSV/PDF/PNG outputs, runs the project integrity checker, and refreshes `docs/output_manifest.csv`.

## Dependency Groups

Most modules use CRAN packages such as `ggplot2`, `patchwork`, `vegan`, `randomForest`, `igraph`, `ggtern`, `microeco`, and `plspm`.

Two modules use Bioconductor packages:

- `07_differential_volcano_heatmap`: `DESeq2`, `ComplexHeatmap`
- `08_kegg_enrichment`: `DESeq2`, `clusterProfiler`, `enrichplot`

The install helper installs only required dependency types (`Depends`, `Imports`, and `LinkingTo`) instead of optional `Suggests`. This avoids failing a public CI run because a package recommends optional helper packages that are not needed by these demo scripts.

The `gghalves` package is installed from the configured CRAN repository first. If it is unavailable from that repository, the helper can try a GitHub fallback. The LEfSe module can use `aplot` for its native combined cladogram layout when `aplot` is available, but falls back to `patchwork` when it is not.

The Mantel module uses `ggcor`. Depending on the R version and repository mirror, `ggcor` may not be available from CRAN. The install helper first tries the configured CRAN repository and then, by default, tries `remotes::install_github("houyunhuang/ggcor")`.

## Install Script Options

The dependency helper can be configured with environment variables:

```bash
CRAN_REPO=https://packagemanager.posit.co/cran/latest Rscript scripts/install_r_dependencies.R
NCPUS=4 Rscript scripts/install_r_dependencies.R
INSTALL_BIOC_PACKAGES=false Rscript scripts/install_r_dependencies.R
INSTALL_CRAN_GITHUB_FALLBACKS=false Rscript scripts/install_r_dependencies.R
INSTALL_GGCOR_FROM_GITHUB=false Rscript scripts/install_r_dependencies.R
```

Use `INSTALL_BIOC_PACKAGES=false` only for partial testing, because modules `07` and `08` need Bioconductor packages to run.

## GitHub Actions

The workflow in `.github/workflows/run-demos.yml` runs the same public reproducibility path:

1. Check out the repository.
2. Set up R.
3. Install Linux system libraries commonly needed by R graphics and bioinformatics packages.
4. Run `Rscript scripts/install_r_dependencies.R`.
5. Regenerate shared toy data.
6. Run all 11 demo modules with `Rscript scripts/run_all_demos.R`.
7. Upload generated toy CSV files, result CSV files, and PDF/PNG figures as a short-lived workflow artifact.

This CI job is intentionally a full check rather than a minimal smoke test. It may take several minutes because Bioconductor packages are relatively large.

The workflow in `.github/workflows/smoke-test.yml` runs the lightweight smoke test. It is intended as a quick structural check and does not replace the full demo workflow.

## Expected Outputs

After a successful full run:

- `data/toy_shared/` contains regenerated simulated toy CSV files.
- Each numbered module writes CSV outputs to its own `results/` folder.
- Each numbered module writes PDF/PNG figures to its own `figures/` folder.

The generated data and figures are demonstration outputs from simulated/desensitized data. They do not represent real study results.
