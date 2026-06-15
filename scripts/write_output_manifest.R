# Purpose:
#   Write a machine-readable manifest of declared module outputs.
#
# Usage:
#   Rscript scripts/write_output_manifest.R
#
# Output:
#   docs/output_manifest.csv

registry_path <- file.path("docs", "module_registry.csv")
output_path <- file.path("docs", "output_manifest.csv")

if (!file.exists(registry_path)) {
  stop("Missing module registry: ", registry_path, call. = FALSE)
}

split_list <- function(value) {
  value <- ifelse(is.na(value), "", value)
  out <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
  out[nzchar(out)]
}

format_mtime <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  format(file.info(path)$mtime, "%Y-%m-%d %H:%M:%S %z")
}

file_type_from_path <- function(path) {
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c("csv", "tsv")) return("result_table")
  if (extension %in% c("png", "pdf")) return("figure")
  "other"
}

registry <- read.csv(registry_path, stringsAsFactors = FALSE, check.names = FALSE)
rows <- list()
row_index <- 1L

for (i in seq_len(nrow(registry))) {
  module <- registry$module_folder[i]
  declared_files <- c(
    split_list(registry$main_result_files[i]),
    split_list(registry$main_figure_files[i])
  )

  for (relative_file in declared_files) {
    path <- file.path(module, relative_file)
    exists <- file.exists(path)
    rows[[row_index]] <- data.frame(
      module = module,
      file_path = path,
      file_type = file_type_from_path(path),
      generated_by = file.path(module, "scripts", "run_demo.R"),
      exists = exists,
      file_size = if (exists) file.info(path)$size else NA_real_,
      last_modified = format_mtime(path),
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}

manifest <- do.call(rbind, rows)
write.csv(manifest, output_path, row.names = FALSE)

missing_outputs <- manifest$file_path[!manifest$exists]
message("Output manifest written to: ", output_path)
message("Declared outputs: ", nrow(manifest))
message("Missing declared outputs: ", length(missing_outputs))
if (length(missing_outputs) > 0) {
  message(paste0("  - ", missing_outputs, collapse = "\n"))
}
