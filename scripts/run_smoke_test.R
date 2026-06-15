# Purpose:
#   Run a lightweight public smoke test without Bioconductor-heavy modules.
#
# Usage:
#   Rscript scripts/run_smoke_test.R

run_command <- function(command, args, working_dir = getwd()) {
  old_dir <- getwd()
  on.exit(setwd(old_dir), add = TRUE)
  setwd(working_dir)
  status <- system2(command, args)
  if (!identical(status, 0L)) {
    stop(
      "Command failed in ",
      working_dir,
      ": ",
      paste(c(command, args), collapse = " "),
      " (exit code ",
      status,
      ")",
      call. = FALSE
    )
  }
}

check_outputs <- function(paths) {
  missing <- paths[!file.exists(paths)]
  empty <- paths[file.exists(paths) & file.info(paths)$size <= 0]

  if (length(missing) > 0 || length(empty) > 0) {
    problems <- character()
    if (length(missing) > 0) {
      problems <- c(problems, paste("Missing outputs:", paste(missing, collapse = ", ")))
    }
    if (length(empty) > 0) {
      problems <- c(problems, paste("Empty outputs:", paste(empty, collapse = ", ")))
    }
    stop(paste(problems, collapse = "\n"), call. = FALSE)
  }
}

smoke_modules <- c(
  "01_rf_correlation_heatmap",
  "02_microbe_env_network",
  "03_ternary_taxa_distribution",
  "04_faprotax_functional_profile"
)

expected_outputs <- list(
  "01_rf_correlation_heatmap" = c(
    "results/correlation_results.csv",
    "results/rf_importance.csv",
    "figures/rf_correlation_combined.pdf",
    "figures/rf_correlation_combined.png"
  ),
  "02_microbe_env_network" = c(
    "results/network_edge_list.csv",
    "results/network_node_list.csv",
    "figures/microbe_env_network.pdf",
    "figures/microbe_env_network.png"
  ),
  "03_ternary_taxa_distribution" = c(
    "results/group_mean_abundance_by_Genus.csv",
    "results/ternary_plotting_table_Genus.csv",
    "figures/ternary_taxa_distribution.pdf",
    "figures/ternary_taxa_distribution.png"
  ),
  "04_faprotax_functional_profile" = c(
    "results/functional_group_summary.csv",
    "results/functional_bubble_plotting_table.csv",
    "figures/faprotax_function_bubble_profile.pdf",
    "figures/faprotax_function_bubble_profile.png"
  )
)

cat("Regenerating shared toy data...\n")
run_command("Rscript", file.path("scripts", "create_shared_toy_data.R"))

for (module in smoke_modules) {
  cat("\nRunning smoke module:", module, "\n")
  run_command("Rscript", file.path("scripts", "run_demo.R"), working_dir = module)
  check_outputs(file.path(module, expected_outputs[[module]]))
  cat("Smoke outputs verified for:", module, "\n")
}

cat("\nRunning project integrity checker...\n")
run_command("Rscript", file.path("scripts", "check_project_integrity.R"))

cat("\nRefreshing output manifest...\n")
run_command("Rscript", file.path("scripts", "write_output_manifest.R"))

cat("\nSmoke test completed successfully.\n")
