#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "tests/test_api_self_tests.R"
script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(file.path(root, "api"))

source("simulation.R.dtd-backup")
source("ddm_inspector.R")
source("helpers.R")
source("templates.R")
source("self_tests.R")

failures <- character()
checks <- 0L
expect_true <- function(condition, message) {
  checks <<- checks + 1L
  if (!isTRUE(condition)) failures <<- c(failures, message)
}

cat("Built-in API self-test suite\n")
catalog <- DDM.self_test_catalog()
expect_true(length(catalog) == 6L, "The installed app must expose six default tests")
expect_true(identical(formals(Simulate.DBP)$disc, 0.0015), "Engine default discrepancy is 0.0015")
expect_true(identical(formals(DDM.simulate_inspected)$disc, 0.0015), "Inspector default matches the engine")
expect_true(identical(names(catalog), c("engine", "schedule", "activation", "signals", "reproducibility", "inspector")),
            "Only DDM tests must be exposed")
expect_true(all(vapply(catalog, function(x) x$model == "DDM", logical(1))), "All tests describe the DDM")

set.seed(810)
rng_before <- .Random.seed
result <- DDM.run_self_tests()
expect_true(identical(rng_before, .Random.seed), "Running installation checks must preserve the caller RNG")
if (!isTRUE(result$success)) print(result$tests)
expect_true(isTRUE(result$success), "All built-in application tests must pass")
expect_true(result$summary$total == 6L, "Self-test summary must report six tests")
expect_true(result$summary$passed == 6L && result$summary$failed == 0L,
            "Self-test summary must report 6/6 passed")
expect_true(length(result$tests) == 6L, "Self-test response must include every result")
expect_true(all(vapply(result$tests, function(x) isTRUE(x$pass), logical(1))),
            "Each individual installed-app test must pass")
expect_true(all(vapply(result$tests, function(x) nzchar(x$observed), logical(1))),
            "Each installed-app test must report its observed value")
expect_true(all(vapply(result$tests, function(x) nzchar(x$expected), logical(1))),
            "Each installed-app test must report its criterion")
expect_true(all(vapply(result$tests, function(x) is.numeric(x$durationSeconds) && x$durationSeconds >= 0,
                       logical(1))),
            "Each installed-app test must report a nonnegative duration")

# Rejected model requests must not be silently reinterpreted or consume RNG.
accepts <- function(body, npes = NULL) tryCatch({ DDM.assert_request(body, npes); TRUE }, error = function(e) FALSE)
expect_true(accepts(list()), "Legacy untagged DDM requests remain supported")
expect_true(accepts(list(model = "DDM")), "DDM name is accepted")
expect_true(!accepts(list(model = "foreign")), "Foreign model is rejected")
expect_true(!accepts(list(modelKind = "foreign")), "Foreign modelKind is rejected")
expect_true(!accepts(list(model = c("dtd", "foreign"))), "Ambiguous model list is rejected")
expect_true(!accepts(list(), data.frame(Layer = "Outcome")), "Foreign layers are rejected")
expect_true(identical(rng_before, .Random.seed), "Request validation preserves random state")
expect_true(all(vapply(get_all_templates(), function(x) !grepl("^carta", x$id), logical(1))),
            "The template catalog contains no retired-model entries")

cat(sprintf("Completed %d checks.\n", checks))
if (length(failures)) {
  cat("FAILURES:\n", paste0("- ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("ALL BUILT-IN API SELF-TESTS PASSED\n")
