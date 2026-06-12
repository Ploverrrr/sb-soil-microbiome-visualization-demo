# Purpose:
#   Install the R packages needed to run all public demo modules.
#
# Usage:
#   Rscript scripts/install_r_dependencies.R
#
# Notes:
#   This script is intentionally separate from the module scripts. Individual
#   demos report missing packages but do not install anything automatically.
#   Run this helper explicitly when preparing a fresh R environment.

default_cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(default_cran_repo) || is.na(default_cran_repo) || default_cran_repo %in% c("", "@CRAN@")) {
  default_cran_repo <- "https://cloud.r-project.org"
}
cran_repo <- Sys.getenv("CRAN_REPO", unset = default_cran_repo)
options(repos = c(CRAN = cran_repo))

truthy <- function(value) {
  tolower(value) %in% c("1", "true", "yes", "y")
}

install_bioc_packages <- truthy(Sys.getenv("INSTALL_BIOC_PACKAGES", unset = "true"))
install_ggcor_from_github <- truthy(Sys.getenv("INSTALL_GGCOR_FROM_GITHUB", unset = "true"))

ncpus <- suppressWarnings(as.integer(Sys.getenv("NCPUS", unset = "")))
if (is.na(ncpus) || ncpus < 1) {
  ncpus <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
}

cran_packages <- c(
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
)

bioc_packages <- c(
  "DESeq2",
  "ComplexHeatmap",
  "clusterProfiler",
  "enrichplot"
)

install_missing_from_cran <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    message("All requested CRAN packages are already installed.")
    return(invisible(TRUE))
  }

  message("Installing CRAN package(s): ", paste(missing, collapse = ", "))
  install.packages(missing, dependencies = TRUE, Ncpus = ncpus)

  still_missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "The following CRAN package(s) are still missing after installation: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

install_missing_from_bioc <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) {
    message("All requested Bioconductor packages are already installed.")
    return(invisible(TRUE))
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager from CRAN.")
    install.packages("BiocManager", dependencies = TRUE, Ncpus = ncpus)
  }

  message("Installing Bioconductor package(s): ", paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE, Ncpus = ncpus)

  still_missing <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still_missing) > 0) {
    stop(
      "The following Bioconductor package(s) are still missing after installation: ",
      paste(still_missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

install_ggcor <- function() {
  if (requireNamespace("ggcor", quietly = TRUE)) {
    message("ggcor is already installed.")
    return(invisible(TRUE))
  }

  message("Trying to install ggcor from the configured CRAN repository.")
  cran_result <- try(
    install.packages("ggcor", dependencies = TRUE, Ncpus = ncpus),
    silent = TRUE
  )

  if (requireNamespace("ggcor", quietly = TRUE)) {
    message("ggcor installed successfully from CRAN.")
    return(invisible(TRUE))
  }

  if (inherits(cran_result, "try-error")) {
    message("CRAN installation attempt for ggcor failed.")
  } else {
    message("ggcor was not available from the configured CRAN repository.")
  }

  if (!install_ggcor_from_github) {
    stop(
      "ggcor is required by 11_vpa_mantel_partitioning but is not installed. ",
      "Set INSTALL_GGCOR_FROM_GITHUB=true to let this script try remotes::install_github('houyunhuang/ggcor'), ",
      "or install ggcor manually following the package maintainer's current instructions.",
      call. = FALSE
    )
  }

  if (!requireNamespace("remotes", quietly = TRUE)) {
    message("Installing remotes from CRAN so ggcor can be installed from GitHub.")
    install.packages("remotes", dependencies = TRUE, Ncpus = ncpus)
  }

  message("Installing ggcor from GitHub: houyunhuang/ggcor")
  remotes::install_github("houyunhuang/ggcor", dependencies = TRUE, upgrade = "never")

  if (!requireNamespace("ggcor", quietly = TRUE)) {
    stop("ggcor is still missing after the GitHub installation attempt.", call. = FALSE)
  }

  invisible(TRUE)
}

message("Using CRAN repository: ", cran_repo)
message("Using Ncpus: ", ncpus)

install_missing_from_cran(cran_packages)

if (install_bioc_packages) {
  install_missing_from_bioc(bioc_packages)
} else {
  message("Skipping Bioconductor package installation because INSTALL_BIOC_PACKAGES is false.")
}

install_ggcor()

all_required <- c(cran_packages, if (install_bioc_packages) bioc_packages else character(), "ggcor")
still_missing <- all_required[!vapply(all_required, requireNamespace, logical(1), quietly = TRUE)]

if (length(still_missing) > 0) {
  stop(
    "Dependency installation did not complete. Missing package(s): ",
    paste(still_missing, collapse = ", "),
    call. = FALSE
  )
}

message("R dependency installation check completed successfully.")
