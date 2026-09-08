#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "tests/run_all.R"
# Rscript encodes spaces as ~+~ in --file on some launch paths (including the
# desktop runner), so decode it before resolving the project root.
script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

test_files <- c(
  "test_api_self_tests.R",
  "test_ddm_inspector.R"
)
for (test_file in test_files) {
  status <- system2("Rscript", shQuote(file.path(root, "tests", test_file)))
  if (status != 0) quit(status = status)
}
cat("All R test suites passed.\n")
