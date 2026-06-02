# Purpose:
#   Run all public demo modules from the repository root.
#
# Notes:
#   Each module's run_demo.R is executed from inside that module directory.
#   Failures are recorded and the runner continues with the next module.

modules <- c(
  "01_rf_correlation_heatmap",
  "02_microbe_env_network",
  "03_ternary_taxa_distribution",
  "04_faprotax_functional_profile",
  "05_lefse_biomarker",
  "06_plspm_mechanism_model",
  "07_differential_volcano_heatmap",
  "08_kegg_enrichment",
  "09_sulfur_gene_contaminant_association",
  "10_alpha_beta_diversity",
  "11_vpa_mantel_partitioning"
)

root_dir <- getwd()
results <- data.frame(
  module = modules,
  status = "pending",
  exit_code = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_along(modules)) {
  module <- modules[i]
  module_dir <- file.path(root_dir, module)
  demo_script <- file.path("scripts", "run_demo.R")

  cat("\n============================================================\n")
  cat("Running module:", module, "\n")
  cat("============================================================\n")

  if (!dir.exists(module_dir)) {
    cat("FAILED: module directory does not exist.\n")
    results$status[i] <- "failed"
    results$exit_code[i] <- NA_integer_
    next
  }

  if (!file.exists(file.path(module_dir, demo_script))) {
    cat("FAILED: scripts/run_demo.R does not exist.\n")
    results$status[i] <- "failed"
    results$exit_code[i] <- NA_integer_
    next
  }

  old_dir <- getwd()
  setwd(module_dir)
  exit_code <- tryCatch(
    system2("Rscript", demo_script),
    error = function(error) {
      cat("FAILED:", conditionMessage(error), "\n")
      1
    }
  )
  setwd(old_dir)

  if (identical(exit_code, 0L)) {
    results$status[i] <- "passed"
    results$exit_code[i] <- 0L
    cat("PASSED:", module, "\n")
  } else {
    results$status[i] <- "failed"
    results$exit_code[i] <- as.integer(exit_code)
    cat("FAILED:", module, "(exit code:", exit_code, ")\n")
  }
}

cat("\n============================================================\n")
cat("Demo run summary\n")
cat("============================================================\n")

passed <- results$module[results$status == "passed"]
failed <- results$module[results$status == "failed"]

cat("Passed modules:", length(passed), "\n")
if (length(passed) > 0) {
  cat(paste0("  - ", passed, collapse = "\n"), "\n")
}

cat("Failed modules:", length(failed), "\n")
if (length(failed) > 0) {
  cat(paste0("  - ", failed, collapse = "\n"), "\n")
}

if (length(failed) > 0) {
  quit(status = 1)
}
