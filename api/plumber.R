library(plumber)
library(jsonlite)

source("simulation.R")
source("helpers.R")
source("templates.R")

#* @apiTitle DiffDiscM Simulation API
#* @apiDescription REST API for the Diffuse Discrepancy Model (DiffDiscM) simulator

#* Enable CORS
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin", "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Accept")

  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  }

  plumber::forward()
}


#* Health check
#* @get /api/health
function() {
  list(status = "ok", timestamp = Sys.time())
}


#* Get all phenomenon templates (metadata only)
#* @get /api/templates
#* @serializer json
function() {
  get_all_templates()
}


#* Get full template data by ID
#* @get /api/templates/<id>
#* @param id The template identifier
#* @serializer json
function(id) {
  get_template_data(id)
}


#* Validate network configuration
#* @post /api/validate
#* @serializer json
function(req) {
  body <- req$body

  tryCatch({
    npes <- as.data.frame(body$npes, stringsAsFactors = FALSE)
    connections <- as.data.frame(body$connections, stringsAsFactors = FALSE)

    # Basic validation
    if (nrow(npes) < 2) {
      return(list(valid = FALSE, error = "At least 2 NPEs are required (including US and D)"))
    }

    if (!"US" %in% npes$NPE) {
      return(list(valid = FALSE, error = "A US (Unconditioned Stimulus) unit is required"))
    }

    if (!any(npes$Layer == "Dopaminergic")) {
      return(list(valid = FALSE, error = "A Dopaminergic (D) unit is required"))
    }

    # Check US->D connection exists with weight 1
    us_d <- connections[connections$PreSinapticNPE == "US" &
                        connections$PostSinapticNPE %in% npes$NPE[npes$Layer == "Dopaminergic"], ]
    if (nrow(us_d) == 0) {
      return(list(valid = FALSE, error = "A connection from US to D with weight 1.0 is required"))
    }

    return(list(valid = TRUE, message = "Network configuration is valid"))
  }, error = function(e) {
    list(valid = FALSE, error = e$message)
  })
}


#* Create timesteps from phases and trials
#* @post /api/create-phases
#* @serializer json
function(req) {
  body <- req$body

  tryCatch({
    phases <- body$contingencies
    trials <- body$trials

    timesteps <- Create.Phases(phases, trials)
    list(
      success = TRUE,
      timesteps = timesteps,
      totalRows = nrow(timesteps),
      phases = unique(timesteps$Phase)
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}


#* Run simulation (legacy — all networks at once, no real progress)
#* @post /api/simulate
#* @serializer json
function(req) {
  body <- req$body

  tryCatch({
    npes <- as.data.frame(body$npes, stringsAsFactors = FALSE)
    connections <- as.data.frame(body$connections, stringsAsFactors = FALSE)
    trials <- body$trials
    contingencies <- body$contingencies
    hasITI <- body$hasITI
    numNetworks <- ifelse(is.null(body$numNetworks), 1, as.integer(body$numNetworks))
    threshold <- ifelse(is.null(body$threshold), "gaussian", body$threshold)
    pupdate <- ifelse(is.null(body$pupdate), "async_random", body$pupdate)
    saveData <- if (!is.null(body$saveData)) body$saveData else list()
    disc <- ifelse(is.null(body$disc), 0.001, as.numeric(body$disc))
    # Ensure correct column names
    if (is.null(colnames(npes)) || !("NPE" %in% colnames(npes))) {
      colnames(npes) <- c("NPE", "Type", "Layer", "Activation", "Temporal.Summation",
                          "Activation.Decay", "mu", "sigma", "logisSigma")
    }
    if (is.null(colnames(connections)) || !("PreSinapticNPE" %in% colnames(connections))) {
      colnames(connections) <- c("PreSinapticNPE", "PostSinapticNPE", "Weight",
                                 "alpha", "beta", "alpha_prime", "beta_prime")
    }

    # Ensure numeric columns
    npes[, 4:9] <- lapply(npes[, 4:9, drop = FALSE], as.numeric)
    connections[, 3:7] <- lapply(connections[, 3:7, drop = FALSE], as.numeric)

    # Create timesteps
    timesteps <- Create.Phases(contingencies, trials)

    if (is.null(hasITI)) {
      hasITI <- rep(FALSE, length(contingencies))
    }

    # Run simulations
    results <- vector("list", numNetworks)
    start_time <- Sys.time()

    for (i in 1:numNetworks) {
      results[[i]] <- Simulate.DBP(
        NPEs = npes,
        Connections = connections,
        TimeSteps = timesteps,
        HasITI = hasITI,
        threshold = threshold,
        disc = disc,
        pupdate = pupdate,
        saveData = saveData
      )
    }

    end_time <- Sys.time()
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Identify unit columns vs connection columns vs signal columns
    all_cols <- colnames(results[[1]])
    meta_cols <- c("Phase", "Trial", "TimeStep")
    signal_cols <- c("dVTA", "dH")
    data_cols <- setdiff(all_cols, c(meta_cols, signal_cols))
    unit_cols <- data_cols[!grepl("-", data_cols)]
    connection_cols <- data_cols[grepl("-", data_cols)]

    list(
      success = TRUE,
      results = results,
      metadata = list(
        numNetworks = numNetworks,
        phases = unique(results[[1]]$Phase),
        units = unit_cols,
        connections = connection_cols,
        signals = signal_cols,
        totalTrials = max(as.numeric(results[[1]]$Trial), na.rm = TRUE),
        totalTimesteps = nrow(results[[1]]),
        duration = round(duration, 2),
        disc = disc
      )
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}


#* Run a single network simulation (used for real-time progress)
#* @post /api/simulate-one
#* @serializer json
function(req) {
  body <- req$body

  tryCatch({
    npes <- as.data.frame(body$npes, stringsAsFactors = FALSE)
    connections <- as.data.frame(body$connections, stringsAsFactors = FALSE)
    trials <- body$trials
    contingencies <- body$contingencies
    hasITI <- body$hasITI
    threshold <- ifelse(is.null(body$threshold), "gaussian", body$threshold)
    pupdate <- ifelse(is.null(body$pupdate), "async_random", body$pupdate)
    saveData <- if (!is.null(body$saveData)) body$saveData else list()
    disc <- ifelse(is.null(body$disc), 0.001, as.numeric(body$disc))

    # Ensure correct column names
    if (is.null(colnames(npes)) || !("NPE" %in% colnames(npes))) {
      colnames(npes) <- c("NPE", "Type", "Layer", "Activation", "Temporal.Summation",
                          "Activation.Decay", "mu", "sigma", "logisSigma")
    }
    if (is.null(colnames(connections)) || !("PreSinapticNPE" %in% colnames(connections))) {
      colnames(connections) <- c("PreSinapticNPE", "PostSinapticNPE", "Weight",
                                 "alpha", "beta", "alpha_prime", "beta_prime")
    }

    # Ensure numeric columns
    npes[, 4:9] <- lapply(npes[, 4:9, drop = FALSE], as.numeric)
    connections[, 3:7] <- lapply(connections[, 3:7, drop = FALSE], as.numeric)

    # Create timesteps
    timesteps <- Create.Phases(contingencies, trials)

    if (is.null(hasITI)) {
      hasITI <- rep(FALSE, length(contingencies))
    }

    # Run ONE simulation
    result <- Simulate.DBP(
      NPEs = npes,
      Connections = connections,
      TimeSteps = timesteps,
      HasITI = hasITI,
      threshold = threshold,
      disc = disc,
      pupdate = pupdate,
      saveData = saveData
    )

    list(success = TRUE, result = result)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}
