#!/usr/bin/env Rscript
# Technical equivalence tests, not new scientific simulations or validation data.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "tests/test_ddm_inspector.R"
script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(root)
source("api/simulation.R.dtd-backup")
source("api/ddm_inspector.R")

checks <- 0L
failures <- character()
expect <- function(condition, message) {
  checks <<- checks + 1L
  if (!isTRUE(condition)) failures <<- c(failures, message)
}
expect_error <- function(expression, pattern, message) {
  error <- tryCatch({ force(expression); "" }, error = conditionMessage)
  expect(nzchar(error) && grepl(pattern, error, fixed = TRUE), message)
}
close <- function(actual, expected) isTRUE(all.equal(actual, expected, tolerance = 1e-14))

npes <- data.frame(
  NPE = c("US", "CS", "S2", "H", "M", "I", "R", "D"),
  Type = c(rep("Excitatory", 5), "Inhibitory", "Excitatory", "Excitatory"),
  Layer = c("US", "PrimarySensory", "AssociativeSensory", "Hippocampal",
            "AssociativeMotor", "AssociativeMotor", "PrimaryMotor", "Dopaminergic"),
  Activation = 0, Temporal.Summation = 0.1, Activation.Decay = 0.1,
  mu = 0.2, sigma = 0.15, logisSigma = 0.1)
connections <- data.frame(
  PreSinapticNPE = c("US", "CS", "CS", "S2", "H", "S2", "M", "I", "M", "R", "US", "US"),
  PostSinapticNPE = c("US", "CS", "S2", "H", "M", "M", "I", "R", "R", "D", "D", "R"),
  Weight = c(1, 1, .6, .3, .2, .4, .5, .3, .6, .2, 1, 1),
  alpha = .5, beta = .1, alpha_prime = .5, beta_prime = .1)
# The external-unit self-connections are only fixture plumbing: externally
# assigned activations ignore them. They avoid the historical engine's 1:0
# zero-input iteration issue in fixed-order, learning-off bookkeeping; that
# unrelated behavior is deliberately not repaired by this UI feature.
trials <- list(
  Acquisition = c("CS,1,US,0,TRUE", "CS,1,US,0,TRUE", "CS,1,US,.8,TRUE", "CS,0,US,0,TRUE"),
  Extinction = c("CS,1,US,0,TRUE", "CS,1,US,0,TRUE", "CS,1,US,0,TRUE", "CS,0,US,0,TRUE"),
  ITI = "CS,0,US,0,TRUE")
branches <- character()
weight_branches <- character()
cases <- 0L

for (threshold in c("gaussian", "beta")) {
  for (pupdate in c("async_random", "sync_random", "async_fixed", "sync_fixed")) {
    for (iti in c(FALSE, TRUE)) {
      for (learning in c(FALSE, TRUE)) {
        cases <- cases + 1L
        suffix <- if (iti) ", True, 2, 3, ITI" else ", False"
        phases <- paste0(c("acquisition, Random, Acquisition, 3", "extinction, Random, Extinction, 3"), suffix)
        set.seed(700 + cases)
        time <- Create.Phases(phases, trials)
        time[, ncol(time)] <- learning
        args <- list(NPEs = npes, Connections = connections, TimeSteps = time,
                     HasITI = rep(iti, 2), threshold = threshold, pupdate = pupdate)
        set.seed(1000 + cases)
        baseline <- do.call(Simulate.DBP, args)
        rng_baseline <- .Random.seed
        set.seed(1000 + cases)
        inspected <- do.call(DDM.simulate_inspected, args)
        label <- paste(threshold, pupdate, iti, learning)
        expect(identical(baseline, inspected$result), paste(label, "result identity"))
        expect(identical(rng_baseline, .Random.seed), paste(label, "RNG identity"))
        trace <- inspected$inspector
        expect(trace$recordedTimesteps == nrow(time) && !trace$truncated, paste(label, "complete trace"))

        for (step in trace$steps) {
          row <- step$rowIndex + 1L
          expect(identical(as.character(step$activationOrder), vapply(step$units, `[[`, "", "name")), "Actual activation order retained")
          expect(length(step$units) == nrow(npes), "Every unit recorded once")
          expect(step$resetApplied == (!iti && step$timestep == 1), "Only original reset marked")
          expect(step$learningEnabled == learning, "Learning flag retained")
          if (!learning) {
            expect(is.null(step$dD) && is.null(step$dH) && length(step$learningOrder) == 0,
                   "Disabled learning signals are not invented zeros")
          } else {
            expect(identical(step$dD, baseline[row, "dVTA"]) && identical(step$dH, baseline[row, "dH"]),
                   "Actual computed signals retained")
          }
          for (unit in step$units) {
            branches <- union(branches, unit$branch)
            expect(identical(unit$activation, baseline[row, unit$name]), "Unit activation equals result")
            if (unit$branch %in% c("external", "unconditional")) {
              expect(is.null(unit$threshold) && is.null(unit$excInput) && is.null(unit$inhInput),
                     "Unused threshold and inputs explicitly null")
            } else {
              logistic <- function(value) 1 / (1 + exp((.5 - value) / unit$logisSigma))
              expect(close(unit$logisticExc, logistic(unit$excInput)), "Excitatory logistic operands match")
              expect(close(unit$logisticInh, logistic(unit$inhInput)), "Inhibitory logistic operands match")
              previous <- logistic(unit$previousExcitatoryInput)
              expected <- switch(unit$branch,
                suprathreshold = unit$logisticExc + unit$temporalSummation * previous * (1 - unit$logisticExc) - unit$logisticInh,
                subthreshold = previous - unit$activationDecay * previous,
                inhibited = 0)
              expect(close(unit$activation, expected), "Recorded activation equation evaluates to actual output")
              expected_branch <- if (unit$logisticExc <= unit$logisticInh) "inhibited" else
                if (unit$logisticExc >= unit$threshold) "suprathreshold" else "subthreshold"
              expect(identical(unit$branch, expected_branch), "Actual random threshold determines recorded branch")
            }
          }
          for (connection in step$connections) {
            weight_branches <- union(weight_branches, connection$branch)
            expect(identical(connection$weightAfter, baseline[row, connection$name]), "Weight equals result")
            expected <- switch(connection$branch,
              potentiation = connection$weightBefore + connection$alpha * connection$capacity *
                connection$postActivation * connection$signal * connection$proportion,
              decrement = connection$weightBefore - connection$beta * connection$weightBefore *
                connection$preActivation * connection$postActivation,
              fixedUS = connection$weightBefore,
              learningOff = connection$weightBefore)
            expect(close(connection$weightUnclipped, expected), "Recorded learning equation evaluates to pre-clip weight")
            expect(close(connection$weightAfter, min(max(expected, 0), 1)), "Recorded clipped weight correct")
            expect(close(connection$deltaWeight, connection$weightAfter - connection$weightBefore), "Actual weight difference retained")
          }
        }

        set.seed(1000 + cases)
        limited <- do.call(DDM.simulate_inspected, c(args, list(inspector = list(maxTimesteps = 2))))
        expect(identical(baseline, limited$result) && identical(rng_baseline, .Random.seed),
               paste(label, "truncation does not truncate simulation or consume RNG"))
        expect(limited$inspector$recordedTimesteps == 2 && limited$inspector$truncated &&
                 identical(limited$inspector$steps, trace$steps[1:2]), "Explicit prefix trace is truthful")
      }
    }
  }
}

expect(all(c("external", "unconditional", "suprathreshold", "subthreshold", "inhibited") %in% branches),
       "All five activation branches covered")
expect(all(c("potentiation", "decrement", "fixedUS", "learningOff") %in% weight_branches),
       "All four weight branches covered")
expect(DDM.inspector_options()$maxTimesteps == 2000L, "Default limit 2000")
expect(DDM.inspector_options(list(maxTimesteps = 10000))$maxTimesteps == 10000L, "Hard limit allowed")
for (invalid in list(0, -1, 10001, Inf, NA_real_, 1.5, "2000", c(1, 2))) {
  expect_error(DDM.inspector_options(list(maxTimesteps = invalid)), "must be an integer", "Invalid trace limit rejected")
}
modified <- Simulate.DBP
body(modified) <- quote(stop("Do not run this modified engine"))
set.seed(19)
rng_before_guard <- .Random.seed
expect_error(DDM.inspector_instrument(modified), "loaded function differs", "Modified loaded engine fails closed")
expect_error(DDM.inspector_instrument(engine_path = "missing-engine-for-inspector-test"), "engine source does not match", "Missing engine fails closed")
expect(identical(rng_before_guard, .Random.seed), "Guards do not consume RNG")
expect(identical(unname(tools::md5sum("api/simulation.R.dtd-backup")), DDM.inspector_engine_md5), "Published engine untouched")

# Exercise the actual /simulate-one handler without opening a server or touching
# the application currently running. Only the original DDM engine is loaded.
api_environment <- new.env(parent = globalenv())
old_directory <- getwd()
setwd(file.path(root, "api"))
sys.source("plumber.R", envir = api_environment)
expressions <- parse("plumber.R")
handler_expression <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("function")) &&
  grepl("DDM.simulate_inspected", paste(deparse(x), collapse = " "), fixed = TRUE), as.list(expressions))
expect(length(handler_expression) == 1L, "One inspected API handler")
handler <- eval(handler_expression[[1]], api_environment)
request_body <- list(npes = npes, connections = connections, trials = trials,
  contingencies = c("acquisition, Random, Acquisition, 2, False"), hasITI = FALSE,
  model = "dtd", inspector = list(enabled = TRUE, maxTimesteps = 3))
set.seed(84)
api_inspected <- handler(list(body = request_body))
api_rng <- .Random.seed
request_body$inspector <- NULL
set.seed(84)
api_regular <- handler(list(body = request_body))
expect(isTRUE(api_inspected$success) && isTRUE(api_regular$success), "Both API modes succeed")
expect(identical(api_inspected$result, api_regular$result) && identical(api_rng, .Random.seed),
       "API inspection preserves schedule, result, and RNG")
expect(is.null(api_regular$inspector) && api_inspected$inspector$recordedTimesteps == 3L,
       "No trace on normal API requests")
serialized <- jsonlite::toJSON(api_inspected$inspector, auto_unbox = TRUE, null = "null")
roundtrip <- jsonlite::fromJSON(serialized, simplifyVector = FALSE)
expect(length(roundtrip$steps) == 3L && is.null(roundtrip$steps[[1]]$units[[1]]$threshold) ==
         is.null(api_inspected$inspector$steps[[1]]$units[[1]]$threshold), "JSON nulls and arrays retained")
expect(grepl('@serializer json list(auto_unbox=TRUE, null="null", digits=NA)',
             paste(readLines("plumber.R"), collapse = "\n"), fixed = TRUE),
       "HTTP route explicitly preserves scalar types and numeric precision")
expect(identical(jsonlite::toJSON(I("single-unit"), auto_unbox = TRUE),
                 structure('["single-unit"]', class = "json")), "Singleton update order remains an array")
http_base <- Sys.getenv("DDM_INSPECTOR_TEST_API", "")
if (nzchar(http_base)) {
  request_body$inspector <- list(enabled = TRUE, maxTimesteps = 3)
  handle <- curl::new_handle(post = TRUE, postfields = jsonlite::toJSON(request_body, auto_unbox = TRUE, dataframe = "columns"),
                             httpheader = "Content-Type: application/json")
  response <- curl::curl_fetch_memory(paste0(http_base, "/api/simulate-one"), handle)
  wire <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  if (!isTRUE(wire$success)) cat("HTTP fixture response:", rawToChar(response$content), "\n")
  expect(response$status_code == 200L && isTRUE(wire$success), "Real HTTP request succeeds")
  expect(is.numeric(wire$inspector$recordedTimesteps) && identical(wire$inspector$truncated, TRUE),
         "HTTP trace scalars are not one-element arrays")
  expect(identical(wire$inspector$steps[[1]]$rowIndex, 0L), "HTTP playback row is a scalar number")
  for (s in wire$inspector$steps) for (u in s$units) {
    expect(is.character(u$name) && is.numeric(u$activation), "HTTP operands are scalars")
    if (!is.null(u$excInput)) expect(close(u$logisticExc, 1 / (1 + exp((.5 - u$excInput) / u$logisSigma))),
                                   "HTTP serializer retains equation precision")
  }
}
setwd(old_directory)

cat(sprintf("DDM inspector: %d factorial cases; %d checks.\n", cases, checks))
if (length(failures)) {
  cat(paste0("FAIL: ", unique(failures), collapse = "\n"), "\n")
  quit(status = 1)
}
cat("ALL DDM INSPECTOR TESTS PASSED: exact results, exact RNG, truthful equations.\n")
