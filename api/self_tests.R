# Small, deterministic installation checks for the original DDM engine.
# Technical fixtures only: no new article data or claims of behavioral validation.
DDM.self_test_catalog <- function() {
  definitions <- list(
    engine = c("DDM engine integrity", "The scientific source matches the audited DDM and accepts observation hooks."),
    schedule = c("Temporal protocol", "Training and extinction keep their trial and timestep coordinates."),
    activation = c("Activations and fixed US connections", "Outputs are finite and bounded; the original fixed US-to-D weight is retained."),
    signals = c("Diffuse learning signals", "The original dVTA/dH signals and learned weights are present and finite."),
    reproducibility = c("Stochastic reproducibility", "Repeating a technical seed reproduces the result and final random state."),
    inspector = c("Inspector transparency", "Recording, including a limited trace, preserves the complete result and random state.")
  )
  lapply(names(definitions), function(id) list(id = id, model = "DDM",
    name = definitions[[id]][1], description = definitions[[id]][2])) |>
    setNames(names(definitions))
}

DDM.self_test_case <- function(meta, expression) {
  started <- proc.time()[["elapsed"]]
  outcome <- tryCatch(force(expression), error = function(e)
    list(pass = FALSE, observed = "Test could not complete", expected = "No runtime error", error = conditionMessage(e)))
  c(meta, list(pass = isTRUE(outcome$pass), observed = as.character(outcome$observed),
    expected = as.character(outcome$expected),
    durationSeconds = round(proc.time()[["elapsed"]] - started, 3),
    error = outcome$error %||% NULL))
}

DDM.run_self_tests <- function() {
  # Tests must not disturb a subsequent user's random simulation sequence.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) saved_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", saved_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  started <- proc.time()[["elapsed"]]
  catalog <- DDM.self_test_catalog()
  cases <- list()
  fixture <- new.env(parent = emptyenv())
  prepare <- function() {
    if (exists("args", envir = fixture, inherits = FALSE)) return(fixture$args)
    template <- get_template_data("extinction")
    # Names come from the shipped DDM template, not a replacement architecture.
    trial_names <- names(template$trials)
    phases <- c(paste0("training, Random, ", trial_names[1], ", 3, False"),
                paste0("extinction, Random, ", trial_names[2], ", 3, False"))
    set.seed(1993)
    time <- Create.Phases(phases, template$trials)
    fixture$args <- list(NPEs = template$npes, Connections = template$connections,
      TimeSteps = time, HasITI = c(FALSE, FALSE), threshold = "gaussian",
      disc = 0.0015, pupdate = "async_random")
    fixture$args
  }
  baseline <- function() {
    args <- prepare()
    if (!exists("result", envir = fixture, inherits = FALSE)) {
      set.seed(1993)
      fixture$result <- do.call(Simulate.DBP, args)
      fixture$rng <- .Random.seed
    }
    fixture$result
  }
  cases$engine <- DDM.self_test_case(catalog$engine, {
    observed_hash <- unname(tools::md5sum(.ddm_inspector_engine_path))
    checked_function <- DDM.inspector_instrument()
    list(pass = identical(observed_hash, DDM.inspector_engine_md5) && is.function(checked_function),
      observed = paste("MD5", observed_hash), expected = paste("MD5", DDM.inspector_engine_md5))
  })
  cases$schedule <- DDM.self_test_case(catalog$schedule, {
    time <- prepare()$TimeSteps
    coordinates <- paste(time$Phase, time$Trial, time$TimeStep)
    trials <- unique(paste(time$Phase, time$Trial))
    list(pass = length(trials) == 6L && !anyDuplicated(coordinates) &&
      identical(unique(as.character(time$Phase)), c("training", "extinction")),
      observed = paste(nrow(time), "timesteps;", length(trials), "trials"),
      expected = "6 distinct trials, two ordered phases, unique temporal coordinates")
  })
  cases$activation <- DDM.self_test_case(catalog$activation, {
    result <- baseline()
    units <- as.matrix(result[, prepare()$NPEs$NPE, drop = FALSE])
    fixed <- result[["US-D"]]
    list(pass = all(is.finite(units)) && all(units >= 0 & units <= 1) &&
      length(fixed) == nrow(result) && all(fixed == 1),
      observed = paste("Activation range", paste(signif(range(units), 6), collapse = " to "),
                       "; US-D =", paste(unique(fixed), collapse = ",")),
      expected = "Finite activations in [0,1]; fixed US-D = 1")
  })
  cases$signals <- DDM.self_test_case(catalog$signals, {
    result <- baseline()
    signals <- as.matrix(result[, c("dVTA", "dH"), drop = FALSE])
    weights <- as.matrix(result[, grepl("-", names(result)), drop = FALSE])
    list(pass = all(is.finite(signals)) && ncol(weights) > 0 &&
      all(is.finite(weights)) && all(weights >= 0 & weights <= 1),
      observed = paste(nrow(signals), "rows of dVTA/dH;", ncol(weights), "weight columns"),
      expected = "Both diffuse signals finite; all stored weights in [0,1]")
  })
  cases$reproducibility <- DDM.self_test_case(catalog$reproducibility, {
    reference <- baseline()
    set.seed(1993)
    repeated <- do.call(Simulate.DBP, prepare())
    exact <- identical(reference, repeated) && identical(fixture$rng, .Random.seed)
    list(pass = exact, observed = paste("Exact result and RNG identity:", exact),
      expected = "Identical result and final random state")
  })
  cases$inspector <- DDM.self_test_case(catalog$inspector, {
    reference <- baseline()
    set.seed(1993)
    observed <- do.call(DDM.simulate_inspected, c(prepare(), list(inspector = list(maxTimesteps = 2))))
    exact <- identical(reference, observed$result) && identical(fixture$rng, .Random.seed)
    list(pass = exact && observed$inspector$recordedTimesteps == 2L &&
      observed$inspector$truncated && observed$inspector$totalTimesteps == nrow(reference),
      observed = paste("Exact result/RNG:", exact, "; recorded 2 of", nrow(reference), "timesteps"),
      expected = "Complete original result and RNG retained despite a two-timestep recording limit")
  })
  passed <- sum(vapply(cases, function(x) isTRUE(x$pass), logical(1)))
  list(success = passed == length(cases),
    summary = list(total = length(cases), passed = passed, failed = length(cases) - passed,
      durationSeconds = round(proc.time()[["elapsed"]] - started, 3)),
    modelVersion = "DDM · DDM-UI 3.2", kind = "DDM installation and regression checks",
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE), tests = unname(cases))
}
