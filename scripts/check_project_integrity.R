# Purpose:
#   Run lightweight integrity checks for the public-safe toy-data demo project.
#
# Usage:
#   Rscript scripts/check_project_integrity.R

message_log <- data.frame(
  level = character(),
  check = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

add_message <- function(level, check, detail) {
  message_log <<- rbind(
    message_log,
    data.frame(level = level, check = check, detail = detail, stringsAsFactors = FALSE)
  )
}

pass <- function(check, detail) add_message("PASS", check, detail)
warn <- function(check, detail) add_message("WARNING", check, detail)
fail <- function(check, detail) add_message("FAIL", check, detail)

split_list <- function(value) {
  value <- ifelse(is.na(value), "", value)
  out <- trimws(unlist(strsplit(value, ";", fixed = TRUE)))
  out[nzchar(out)]
}

read_csv_required <- function(path, check_name) {
  if (!file.exists(path)) {
    fail(check_name, paste("Missing file:", path))
    return(NULL)
  }
  tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(error) {
      fail(check_name, paste("Could not read", path, "-", conditionMessage(error)))
      NULL
    }
  )
}

check_set_equal <- function(name, left, right, left_name, right_name) {
  missing_from_right <- setdiff(left, right)
  missing_from_left <- setdiff(right, left)

  if (length(missing_from_right) == 0 && length(missing_from_left) == 0) {
    pass(name, paste(left_name, "and", right_name, "match."))
    return(invisible(TRUE))
  }

  details <- character()
  if (length(missing_from_right) > 0) {
    details <- c(
      details,
      paste0(left_name, " values missing from ", right_name, ": ", paste(head(missing_from_right, 10), collapse = ", "))
    )
  }
  if (length(missing_from_left) > 0) {
    details <- c(
      details,
      paste0(right_name, " values missing from ", left_name, ": ", paste(head(missing_from_left, 10), collapse = ", "))
    )
  }
  fail(name, paste(details, collapse = " | "))
  invisible(FALSE)
}

required_shared_files <- file.path(
  "data", "toy_shared",
  c(
    "sample_metadata.csv",
    "environmental_variables.csv",
    "taxonomy_table.csv",
    "abundance_table.csv",
    "functional_annotation_table.csv"
  )
)

missing_shared <- required_shared_files[!file.exists(required_shared_files)]
if (length(missing_shared) == 0) {
  pass("shared toy data files", "All expected shared toy-data CSV files exist.")
} else {
  fail("shared toy data files", paste("Missing:", paste(missing_shared, collapse = ", ")))
}

sample_metadata <- read_csv_required(file.path("data", "toy_shared", "sample_metadata.csv"), "read sample metadata")
environmental_variables <- read_csv_required(file.path("data", "toy_shared", "environmental_variables.csv"), "read environmental variables")
taxonomy_table <- read_csv_required(file.path("data", "toy_shared", "taxonomy_table.csv"), "read taxonomy table")
abundance_table <- read_csv_required(file.path("data", "toy_shared", "abundance_table.csv"), "read abundance table")
functional_annotation <- read_csv_required(file.path("data", "toy_shared", "functional_annotation_table.csv"), "read functional annotation")

if (!is.null(sample_metadata) && "sample_id" %in% colnames(sample_metadata)) {
  if (anyDuplicated(sample_metadata$sample_id) == 0) {
    pass("sample metadata uniqueness", "sample_metadata.csv has unique sample_id values.")
  } else {
    fail("sample metadata uniqueness", "sample_metadata.csv contains duplicated sample_id values.")
  }
} else {
  fail("sample metadata schema", "sample_metadata.csv must contain sample_id.")
}

if (!is.null(abundance_table) && "feature_id" %in% colnames(abundance_table)) {
  if (anyDuplicated(abundance_table$feature_id) == 0) {
    pass("abundance feature uniqueness", "abundance_table.csv has unique feature_id values.")
  } else {
    fail("abundance feature uniqueness", "abundance_table.csv contains duplicated feature_id values.")
  }
} else {
  fail("abundance schema", "abundance_table.csv must contain feature_id.")
}

if (!is.null(taxonomy_table) && "feature_id" %in% colnames(taxonomy_table)) {
  if (anyDuplicated(taxonomy_table$feature_id) == 0) {
    pass("taxonomy feature uniqueness", "taxonomy_table.csv has unique feature_id values.")
  } else {
    fail("taxonomy feature uniqueness", "taxonomy_table.csv contains duplicated feature_id values.")
  }
} else {
  fail("taxonomy schema", "taxonomy_table.csv must contain feature_id.")
}

if (!is.null(sample_metadata) && !is.null(environmental_variables) &&
    all(c("sample_id") %in% colnames(sample_metadata)) &&
    all(c("sample_id") %in% colnames(environmental_variables))) {
  check_set_equal(
    "metadata/environment sample_id relationship",
    sample_metadata$sample_id,
    environmental_variables$sample_id,
    "sample_metadata.csv",
    "environmental_variables.csv"
  )
}

if (!is.null(sample_metadata) && !is.null(abundance_table) &&
    "sample_id" %in% colnames(sample_metadata) && "feature_id" %in% colnames(abundance_table)) {
  abundance_samples <- setdiff(colnames(abundance_table), "feature_id")
  check_set_equal(
    "metadata/abundance sample_id relationship",
    sample_metadata$sample_id,
    abundance_samples,
    "sample_metadata.csv",
    "abundance_table.csv columns"
  )
}

if (!is.null(sample_metadata) && !is.null(functional_annotation) &&
    "sample_id" %in% colnames(sample_metadata) && "sample_id" %in% colnames(functional_annotation)) {
  check_set_equal(
    "metadata/function sample_id relationship",
    sample_metadata$sample_id,
    unique(functional_annotation$sample_id),
    "sample_metadata.csv",
    "functional_annotation_table.csv"
  )
}

if (!is.null(abundance_table) && !is.null(taxonomy_table) &&
    "feature_id" %in% colnames(abundance_table) && "feature_id" %in% colnames(taxonomy_table)) {
  check_set_equal(
    "abundance/taxonomy feature_id relationship",
    abundance_table$feature_id,
    taxonomy_table$feature_id,
    "abundance_table.csv",
    "taxonomy_table.csv"
  )
}

if (!is.null(abundance_table) && !is.null(functional_annotation) &&
    "feature_id" %in% colnames(abundance_table) && "feature_id" %in% colnames(functional_annotation)) {
  missing_function_features <- setdiff(unique(functional_annotation$feature_id), abundance_table$feature_id)
  if (length(missing_function_features) == 0) {
    pass("functional feature_id relationship", "All functional feature_id values exist in abundance_table.csv.")
  } else {
    fail(
      "functional feature_id relationship",
      paste("Functional feature_id values missing from abundance_table.csv:", paste(head(missing_function_features, 10), collapse = ", "))
    )
  }
}

registry_path <- file.path("docs", "module_registry.csv")
module_registry <- read_csv_required(registry_path, "read module registry")
if (!is.null(module_registry)) {
  expected_columns <- c(
    "module_id",
    "module_folder",
    "main_capability",
    "main_input_tables",
    "main_result_files",
    "main_figure_files",
    "required_packages",
    "uses_bioconductor",
    "scientific_boundary"
  )
  missing_columns <- setdiff(expected_columns, colnames(module_registry))
  if (length(missing_columns) == 0) {
    pass("module registry schema", "docs/module_registry.csv has the expected columns.")
  } else {
    fail("module registry schema", paste("Missing columns:", paste(missing_columns, collapse = ", ")))
  }
} else {
  module_registry <- data.frame(module_folder = Sys.glob("[0-9][0-9]_*"), stringsAsFactors = FALSE)
}

module_folders <- module_registry$module_folder
numbered_folders <- Sys.glob("[0-9][0-9]_*")
check_set_equal(
  "registered numbered modules",
  sort(module_folders),
  sort(numbered_folders),
  "docs/module_registry.csv",
  "repository numbered folders"
)

for (module in module_folders) {
  readme_path <- file.path(module, "README.md")
  run_demo_path <- file.path(module, "scripts", "run_demo.R")

  if (file.exists(readme_path)) {
    pass("module README", paste(module, "has README.md."))
  } else {
    fail("module README", paste(module, "is missing README.md."))
  }

  if (file.exists(run_demo_path)) {
    pass("module run_demo", paste(module, "has scripts/run_demo.R."))
  } else {
    fail("module run_demo", paste(module, "is missing scripts/run_demo.R."))
  }

  for (folder_name in c("results", "figures")) {
    folder_path <- file.path(module, folder_name)
    if (dir.exists(folder_path)) {
      pass("output folder", paste(folder_path, "exists."))
    } else {
      created <- tryCatch(dir.create(folder_path, recursive = TRUE), warning = function(w) FALSE, error = function(e) FALSE)
      if (created && dir.exists(folder_path)) {
        warn("output folder", paste(folder_path, "was missing and was created."))
      } else {
        fail("output folder", paste(folder_path, "does not exist and could not be created."))
      }
    }
  }
}

script_files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
script_files <- script_files[!grepl("^\\./?\\.git/", script_files)]
script_files <- script_files[basename(script_files) != "check_project_integrity.R"]
script_patterns <- c(
  "/Users/",
  "/home/",
  "C:/",
  "~/",
  "_private_original",
  "raw/",
  "private/",
  "original_data",
  "raw_data",
  "real_data",
  "private_data",
  "data/raw",
  "data/private"
)

path_hits <- character()
for (script_file in script_files) {
  lines <- readLines(script_file, warn = FALSE)
  hit_index <- unique(unlist(lapply(script_patterns, function(pattern) grep(pattern, lines, fixed = TRUE))))
  if (length(hit_index) > 0) {
    path_hits <- c(path_hits, paste0(script_file, ":", hit_index, ": ", trimws(lines[hit_index])))
  }
}

if (length(path_hits) == 0) {
  pass("R script path hygiene", "No absolute local paths or raw/private/original-data path references were found in R scripts.")
} else {
  fail("R script path hygiene", paste(path_hits, collapse = "\n"))
}

if (!is.null(module_registry) && all(c("module_folder", "main_result_files", "main_figure_files") %in% colnames(module_registry))) {
  for (i in seq_len(nrow(module_registry))) {
    module <- module_registry$module_folder[i]
    declared_files <- c(split_list(module_registry$main_result_files[i]), split_list(module_registry$main_figure_files[i]))
    declared_files <- file.path(module, declared_files)
    missing_declared <- declared_files[!file.exists(declared_files)]
    if (length(missing_declared) == 0) {
      pass("declared output presence", paste(module, "has all declared output files currently present."))
    } else {
      warn(
        "declared output presence",
        paste(module, "is missing declared output file(s):", paste(missing_declared, collapse = ", "))
      )
    }
  }
}

cat("\nProject integrity check summary\n")
cat("================================\n")
for (level in c("FAIL", "WARNING", "PASS")) {
  subset_log <- message_log[message_log$level == level, , drop = FALSE]
  cat(level, ":", nrow(subset_log), "\n", sep = "")
  if (nrow(subset_log) > 0) {
    for (i in seq_len(nrow(subset_log))) {
      cat("  - [", subset_log$check[i], "] ", subset_log$detail[i], "\n", sep = "")
    }
  }
}

fail_count <- sum(message_log$level == "FAIL")
warning_count <- sum(message_log$level == "WARNING")

cat("\nFinal status: ")
if (fail_count > 0) {
  cat("FAIL\n")
  quit(status = 1)
} else if (warning_count > 0) {
  cat("WARNING\n")
  quit(status = 0)
} else {
  cat("PASS\n")
  quit(status = 0)
}
