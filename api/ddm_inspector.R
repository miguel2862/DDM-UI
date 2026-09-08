# Read-only pedagogical instrumentation of the published DDM engine.
# The scientific source is not edited. Hooks are inserted into a private AST
# copy, and every expected insertion is checked before that copy can execute.

.ddm_inspector_source <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.ddm_inspector_directory <- if (is.null(.ddm_inspector_source)) getwd() else
  dirname(normalizePath(.ddm_inspector_source, mustWork = TRUE))
.ddm_inspector_engine_path <- file.path(.ddm_inspector_directory, "simulation.R.dtd-backup")
rm(.ddm_inspector_source, .ddm_inspector_directory)

DDM.inspector_engine_hash <- "fea83a5b6458c421cfe55756f791cc052a9ce2855b54b843c121c2191d7d2871"
DDM.inspector_engine_md5 <- "24c54494063d7d5445aac05f054e701a"

DDM.inspector_options <- function(options = list()) {
  if (!is.list(options)) stop("inspector must be an object")
  limit <- options$maxTimesteps
  if (is.null(limit)) limit <- 2000L
  if (length(limit) != 1L || !is.numeric(limit) || !is.finite(limit) ||
      limit != floor(limit) || limit < 1 || limit > 10000) {
    stop("inspector.maxTimesteps must be an integer between 1 and 10000")
  }
  list(maxTimesteps = as.integer(limit))
}

DDM.inspector_instrument <- function(engine = Simulate.DBP,
                                      engine_path = .ddm_inspector_engine_path) {
  if (!file.exists(engine_path) ||
      !identical(unname(tools::md5sum(engine_path)), DDM.inspector_engine_md5)) {
    stop("DDM inspector disabled: the engine source does not match the audited version")
  }
  reference <- new.env(parent = environment(engine))
  sys.source(engine_path, envir = reference)
  if (!identical(body(engine), body(reference$Simulate.DBP)) ||
      !identical(formals(engine), formals(reference$Simulate.DBP))) {
    stop("DDM inspector disabled: the loaded function differs from the audited engine")
  }

  counts <- c(step_loop = 0L, activation_order = 0L, fixed_order = 0L,
              unit_start = 0L, activation_branch = 0L, signals = 0L,
              connection_loop = 0L, clip = 0L, frozen_loop = 0L)
  append_block <- function(expr, before = NULL, after = NULL) {
    statements <- if (is.call(expr) && identical(expr[[1]], as.name("{")))
      as.list(expr)[-1] else list(expr)
    as.call(c(list(as.name("{")), if (!is.null(before)) list(before),
              statements, if (!is.null(after)) list(after)))
  }
  visit <- function(expr) {
    if (!is.call(expr)) return(expr)
    original <- expr
    for (index in seq_along(expr)[-1]) {
      # Missing formal arguments in arbitrary expressions must remain missing.
      if (!identical(expr[[index]], quote(expr = ))) expr[index] <- list(visit(expr[[index]]))
    }
    if (identical(original[[1]], as.name("for")) &&
        identical(original[[2]], as.name("ts")) &&
        identical(original[[3]], quote(1:nrow(TimeSteps)))) {
      counts[["step_loop"]] <<- counts[["step_loop"]] + 1L
      expr[[4]] <- append_block(expr[[4]],
        quote(.ddm_recorder$begin(environment())),
        quote(.ddm_recorder$finish(environment())))
    }
    if (identical(original, quote(scrambledNPEs <- sample(1:length(network), length(network), replace = F)))) {
      counts[["activation_order"]] <<- counts[["activation_order"]] + 1L
      expr <- append_block(expr, after = quote(.ddm_recorder$order(environment())))
    }
    if (identical(original, quote(scrambledNPEs <- 1:length(network)))) {
      counts[["fixed_order"]] <<- counts[["fixed_order"]] + 1L
      expr <- append_block(expr, after = quote(.ddm_recorder$order(environment())))
    }
    if (identical(original, quote(network[[i]]@PreviousExcitatoryInput <- network[[i]]@ExcitatoryInput))) {
      counts[["unit_start"]] <<- counts[["unit_start"]] + 1L
      expr <- append_block(expr, after = quote(.ddm_recorder$unitStart(environment())))
    }
    if (identical(original[[1]], as.name("if")) &&
        identical(original[[2]], as.name("npe.is.unconditionally.activated"))) {
      counts[["activation_branch"]] <<- counts[["activation_branch"]] + 1L
      expr <- append_block(expr, after = quote(.ddm_recorder$unitFinish(environment())))
    }
    if (identical(original, quote(PreviousdCA1 <- dH))) {
      counts[["signals"]] <<- counts[["signals"]] + 1L
      expr <- append_block(expr, after = quote(.ddm_recorder$signals(environment())))
    }
    if (identical(original[[1]], as.name("for")) &&
        identical(original[[2]], as.name("j")) &&
        identical(original[[3]], as.name("scrambledConnections"))) {
      counts[["connection_loop"]] <<- counts[["connection_loop"]] + 1L
      expr[[4]] <- append_block(expr[[4]],
        quote(.ddm_recorder$connectionStart(environment())))
    }
    if (identical(original, quote(network[[i]]@InputConnections[[j]]@weight <- min(max(network[[i]]@InputConnections[[j]]@weight, 0), 1)))) {
      counts[["clip"]] <<- counts[["clip"]] + 1L
      expr <- append_block(expr,
        quote(.ddm_recorder$beforeClip(environment())),
        quote(.ddm_recorder$connectionFinish(environment())))
    }
    if (identical(original[[1]], as.name("for")) &&
        identical(original[[2]], as.name("j")) &&
        identical(original[[3]], as.name("scrambleConnections"))) {
      counts[["frozen_loop"]] <<- counts[["frozen_loop"]] + 1L
      expr[[4]] <- append_block(expr[[4]],
        quote(.ddm_recorder$frozenConnection(environment())))
    }
    expr
  }
  instrumented <- engine
  body(instrumented) <- visit(body(engine))
  expected <- c(step_loop = 1L, activation_order = 2L, fixed_order = 1L,
                unit_start = 1L, activation_branch = 1L, signals = 1L,
                connection_loop = 1L, clip = 1L, frozen_loop = 1L)
  if (!identical(counts, expected)) {
    stop("DDM inspector disabled: instrumentation anchors changed (",
         paste(names(counts), counts, collapse = ", "), ")")
  }
  formals(instrumented) <- c(formals(instrumented), alist(.ddm_recorder = ))
  instrumented
}

DDM.inspector_recorder <- function(max_timesteps, total_timesteps) {
  state <- new.env(parent = emptyenv())
  state$steps <- vector("list", min(max_timesteps, total_timesteps))
  state$active <- FALSE
  state$stage <- "activation"
  state$step <- NULL
  state$connection <- NULL

  begin <- function(frame) {
    state$active <- frame$ts <= max_timesteps
    if (!state$active) return(invisible(NULL))
    state$stage <- "activation"
    phase_index <- frame$current.phase
    if (frame$ts > 1 && frame$TimeSteps[frame$ts - 1, 1] != frame$TimeSteps[frame$ts, 1]) {
      phase_index <- phase_index + 1L
    }
    state$step <- list(
      rowIndex = as.integer(frame$ts - 1L),
      phase = as.character(frame$TimeSteps[frame$ts, "Phase"]),
      trial = as.integer(frame$TimeSteps[frame$ts, "Trial"]),
      timestep = as.integer(frame$TimeSteps[frame$ts, "TimeStep"]),
      resetApplied = isTRUE(frame$TimeSteps[frame$ts, 3] == 1 && !frame$HasITI[phase_index]),
      learningEnabled = as.logical(frame$TimeSteps[frame$ts, ncol(frame$TimeSteps)]),
      activationOrder = character(), learningOrder = character(),
      dD = NULL, dH = NULL, previousDH = frame$PreviousdCA1,
      units = list(), connections = list())
  }
  order <- function(frame) {
    if (!state$active) return(invisible(NULL))
    field <- if (state$stage == "activation") "activationOrder" else "learningOrder"
    state$step[[field]] <- names(frame$network)[frame$scrambledNPEs]
  }
  unit_record <- function(frame, external = FALSE) {
    npe <- frame$network[[frame$i]]
    forced <- !external && frame$npe.is.unconditionally.activated
    used <- !external && !forced
    branch <- if (external) "external" else if (forced) "unconditional" else
      if (frame$p_epsp <= frame$p_ipsp) "inhibited" else
      if (frame$p_epsp >= npe@Threshold) "suprathreshold" else "subthreshold"
    list(name = npe@Name, layer = npe@Layer,
         order = match(npe@Name, state$step$activationOrder), branch = branch,
         previousActivation = npe@PreviousActivation, activation = npe@Activation,
         previousExcitatoryInput = npe@PreviousExcitatoryInput,
         excInput = if (used) npe@ExcitatoryInput else NULL,
         inhInput = if (used) npe@InhibitoryInput else NULL,
         logisticExc = if (used) frame$p_epsp else NULL,
         logisticInh = if (used) frame$p_ipsp else NULL,
         threshold = if (used) npe@Threshold else NULL,
         mu = npe@mu, sigma = npe@sigma, logisSigma = npe@logisSigma,
         temporalSummation = npe@TemporalSummation, activationDecay = npe@ActivationDecay)
  }
  unitStart <- function(frame) {
    if (!state$active) return(invisible(NULL))
    if (frame$network[[frame$i]]@Layer %in% c("US", "PrimarySensory")) {
      state$step$units[[length(state$step$units) + 1L]] <- unit_record(frame, TRUE)
    }
  }
  unitFinish <- function(frame) {
    if (!state$active) return(invisible(NULL))
    state$step$units[[length(state$step$units) + 1L]] <- unit_record(frame)
  }
  signals <- function(frame) {
    if (!state$active) return(invisible(NULL))
    state$stage <- "learning"
    state$step$dD <- frame$dD
    state$step$dH <- frame$dH
  }
  connectionStart <- function(frame, frozen = FALSE) {
    if (!state$active) return(invisible(NULL))
    post <- frame$network[[frame$i]]
    connection <- post@InputConnections[[frame$j]]
    pre <- frame$network[[connection@preSinapticNPE]]
    fixed <- pre@Layer == "US"
    increasing <- !frozen && !fixed && frame$d >= frame$disc
    excitatory <- pre@Type == "Excitatory"
    denominator <- if (excitatory) post@ExcitatoryInput else post@InhibitoryInput
    state$connection <- list(
      name = connection@Name, pre = pre@Name, post = post@Name,
      order = length(state$step$connections) + 1L, preType = pre@Type,
      preActivation = pre@Activation, postActivation = post@Activation,
      weightBefore = connection@weight, weightUnclipped = connection@weight,
      weightAfter = connection@weight, deltaWeight = 0,
      branch = if (frozen) "learningOff" else if (fixed) "fixedUS" else
        if (increasing) "potentiation" else "decrement",
      signalKind = if (frozen || fixed) NULL else
        if (post@Layer %in% c("AssociativeSensory", "Hippocampal")) "dH" else "dD",
      signal = if (frozen || fixed) NULL else frame$d,
      disc = frame$disc,
      capacity = if (increasing) post@r[if (excitatory) 1 else 2] else NULL,
      proportion = NULL,
      alpha = if (increasing) if (excitatory) connection@alpha else connection@alpha_prime else NULL,
      beta = if (!frozen && !fixed && !increasing)
        if (excitatory) connection@beta else connection@beta_prime else NULL,
      inputDenominator = if (increasing) denominator else NULL)
  }
  beforeClip <- function(frame) {
    if (!state$active) return(invisible(NULL))
    connection <- frame$network[[frame$i]]@InputConnections[[frame$j]]
    state$connection$weightUnclipped <- connection@weight
    if (state$connection$branch == "potentiation") state$connection$proportion <- connection@p
  }
  connectionFinish <- function(frame) {
    if (!state$active) return(invisible(NULL))
    state$connection$weightAfter <- frame$network[[frame$i]]@InputConnections[[frame$j]]@weight
    state$connection$deltaWeight <- state$connection$weightAfter - state$connection$weightBefore
    state$step$connections[[length(state$step$connections) + 1L]] <- state$connection
  }
  frozenConnection <- function(frame) {
    if (!state$active) return(invisible(NULL))
    connectionStart(frame, TRUE)
    connectionFinish(frame)
  }
  finish <- function(frame) {
    if (!state$active) return(invisible(NULL))
    state$steps[[frame$ts]] <- state$step
  }
  final <- function() {
    # Preserve array shape even for a one-unit order under auto_unbox=TRUE.
    steps <- lapply(state$steps, function(step) {
      step$activationOrder <- I(step$activationOrder)
      step$learningOrder <- I(step$learningOrder)
      step
    })
    list(schemaVersion = 1L, engine = "DDM", engineHash = DDM.inspector_engine_hash,
         recordedTimesteps = length(state$steps), totalTimesteps = total_timesteps,
         truncated = total_timesteps > max_timesteps, steps = steps)
  }
  list(begin = begin, order = order, unitStart = unitStart, unitFinish = unitFinish,
       signals = signals, connectionStart = connectionStart, beforeClip = beforeClip,
       connectionFinish = connectionFinish, frozenConnection = frozenConnection,
       finish = finish, final = final)
}

DDM.simulate_inspected <- function(NPEs, Connections, TimeSteps, HasITI,
                                  threshold = "gaussian", saveData = list(),
                                  disc = 0.0015, pupdate = "async_random",
                                  inspector = list()) {
  options <- DDM.inspector_options(inspector)
  instrumented <- DDM.inspector_instrument()
  recorder <- DDM.inspector_recorder(options$maxTimesteps, nrow(TimeSteps))
  result <- instrumented(NPEs = NPEs, Connections = Connections, TimeSteps = TimeSteps,
                         HasITI = HasITI, threshold = threshold, saveData = saveData,
                         disc = disc, pupdate = pupdate, .ddm_recorder = recorder)
  list(result = result, inspector = recorder$final())
}
