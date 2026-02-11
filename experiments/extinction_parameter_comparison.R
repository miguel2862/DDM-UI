# =============================================================================
# Systematic Parameter Comparison for Extinction
# Compares: DDM-UI current, Donahoe et al 1993, Burgos 2007, and hybrids
# =============================================================================

# Source the simulation engine
source("/Users/miguel/Documents/DDM-UI/api/simulation.R")
source("/Users/miguel/Documents/DDM-UI/api/helpers.R")

NUM_NETWORKS <- 20
TRAINING_TRIALS <- 100
EXTINCTION_TRIALS <- 100

# =============================================================================
# Helper: Build extinction experiment with specified parameters
# =============================================================================
build_extinction_experiment <- function(
  beta_val = 0.12,
  decay_val = 0.1,
  mu_thresh = 0.2,
  sigma_thresh = 0.15,
  threshold_type = "gaussian",  # "gaussian" or "beta"
  init_weight = 0.1,
  disc_val = 0.001
) {
  npes <- data.frame(
    NPE = c("US", "D", "S1", "S..1", "H1", "M..1", "M.1"),
    Type = rep("Excitatory", 7),
    Layer = c("US", "Dopaminergic", "PrimarySensory", "AssociativeSensory",
              "Hippocampal", "AssociativeMotor", "PrimaryMotor"),
    Activation = rep(0, 7),
    Temporal.Summation = rep(0.1, 7),
    Activation.Decay = rep(decay_val, 7),
    mu = rep(mu_thresh, 7),
    sigma = rep(sigma_thresh, 7),
    logisSigma = rep(0.1, 7),
    stringsAsFactors = FALSE
  )

  connections <- data.frame(
    PreSinapticNPE = c("S1", "S..1", "S..1", "M..1", "M..1", "US"),
    PostSinapticNPE = c("S..1", "H1", "M..1", "D", "M.1", "D"),
    Weight = c(init_weight, init_weight, init_weight, init_weight, init_weight, 1.0),
    alpha = rep(0.5, 6),
    beta = rep(beta_val, 6),
    alpha_prime = rep(0.5, 6),
    beta_prime = rep(beta_val, 6),
    stringsAsFactors = FALSE
  )

  trials <- list(
    Training = c(
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,1.00,S1,1.00,True"
    ),
    Extinction = c(
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True",
      "US,0.00,S1,1.00,True"
    )
  )

  contingencies <- c(
    paste0("training, Random, Training, ", TRAINING_TRIALS, ", False"),
    paste0("extinction, Random, Extinction, ", EXTINCTION_TRIALS, ", False")
  )

  hasITI <- c(FALSE, FALSE)

  list(
    npes = npes,
    connections = connections,
    trials = trials,
    contingencies = contingencies,
    hasITI = hasITI,
    threshold_type = threshold_type,
    disc_val = disc_val
  )
}

# =============================================================================
# Helper: Run N networks and collect M.1 activation at last extinction trial
# =============================================================================
run_condition <- function(experiment, n_networks = NUM_NETWORKS) {
  phases <- Create.Phases(experiment$contingencies, experiment$trials)

  all_units <- experiment$npes$NPE
  all_conns <- paste(experiment$connections$PreSinapticNPE,
                     experiment$connections$PostSinapticNPE, sep = "-")
  all_conns <- all_conns[experiment$connections$Weight < 1.0]  # exclude US-D

  saveData <- list(
    Elements = c(all_units, all_conns),
    TimeSteps = c(5)  # Only save last timestep of each trial
  )

  results <- list()
  for (net in 1:n_networks) {
    tryCatch({
      res <- Simulate.DBP(
        NPEs = experiment$npes,
        Connections = experiment$connections,
        TimeSteps = phases,
        HasITI = experiment$hasITI,
        threshold = experiment$threshold_type,
        saveData = saveData,
        disc = experiment$disc_val,
        pupdate = "async_random"
      )
      results[[net]] <- res
    }, error = function(e) {
      cat("  Network", net, "failed:", conditionMessage(e), "\n")
      results[[net]] <<- NULL
    })
  }

  results <- Filter(Negate(is.null), results)
  return(results)
}

# =============================================================================
# Helper: Analyze results — extinction success rate + summary stats
# =============================================================================
analyze_results <- function(results, label) {
  n <- length(results)
  if (n == 0) {
    cat(sprintf("\n=== %s === NO RESULTS\n", label))
    return(NULL)
  }

  # Get M.1 activation at last training trial and last extinction trial
  m1_end_training <- numeric(n)
  m1_end_extinction <- numeric(n)
  # Also get key weight at end of extinction
  sm_weight_end <- numeric(n)

  for (i in 1:n) {
    df <- results[[i]]
    training_rows <- df[tolower(df$Phase) == "training", ]
    extinction_rows <- df[tolower(df$Phase) == "extinction", ]

    if (nrow(training_rows) > 0) {
      m1_end_training[i] <- training_rows[nrow(training_rows), "M.1"]
    }
    if (nrow(extinction_rows) > 0) {
      m1_end_extinction[i] <- extinction_rows[nrow(extinction_rows), "M.1"]

      # Try to get S..1-M..1 weight
      if ("S..1-M..1" %in% colnames(extinction_rows)) {
        sm_weight_end[i] <- extinction_rows[nrow(extinction_rows), "S..1-M..1"]
      }
    }
  }

  # Define extinction success: M.1 < 0.1 at end of extinction
  extinguished <- m1_end_extinction < 0.1
  success_rate <- mean(extinguished) * 100

  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  Networks:           %d\n", n))
  cat(sprintf("  Extinction success: %.0f%% (%d/%d)\n", success_rate, sum(extinguished), n))
  cat(sprintf("  M.1 end training:   mean=%.3f, sd=%.3f\n", mean(m1_end_training), sd(m1_end_training)))
  cat(sprintf("  M.1 end extinction: mean=%.3f, sd=%.3f\n", mean(m1_end_extinction), sd(m1_end_extinction)))
  cat(sprintf("  M.1 ext. range:     [%.3f, %.3f]\n", min(m1_end_extinction), max(m1_end_extinction)))
  cat(sprintf("  S..1-M..1 end ext:  mean=%.3f, sd=%.3f\n", mean(sm_weight_end), sd(sm_weight_end)))
  cat(sprintf("  Extinguished M.1:   %s\n",
              paste(round(m1_end_extinction[extinguished], 3), collapse=", ")))
  cat(sprintf("  NOT extinct. M.1:   %s\n",
              paste(round(m1_end_extinction[!extinguished], 3), collapse=", ")))

  return(list(
    label = label,
    n = n,
    success_rate = success_rate,
    m1_training = m1_end_training,
    m1_extinction = m1_end_extinction,
    sm_weight = sm_weight_end,
    extinguished = extinguished
  ))
}

# =============================================================================
# Define experimental conditions
# =============================================================================

conditions <- list(

  # --- Condition 1: DDM-UI Current (baseline) ---
  list(
    label = "1. DDM-UI Current (β=0.12, κ=0.1, θ~Beta(0.2,0.15), disc=0.001)",
    params = list(beta_val=0.12, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 2: Donahoe et al 1993 (original) ---
  # β=0.035, κ=0.05, θ~Gaussian(0,1), w0=0.01
  list(
    label = "2. Donahoe 1993 (β=0.035, κ=0.05, θ~Gauss(0,1), w0=0.01)",
    params = list(beta_val=0.035, decay_val=0.05, mu_thresh=0.0, sigma_thresh=1.0,
                  threshold_type="gaussian", init_weight=0.01, disc_val=0.001)
  ),

  # --- Condition 3: Burgos 2007 ---
  # β=0.1, w0=0.01
  list(
    label = "3. Burgos 2007 (β=0.1, θ~Gauss(0.2,0.15), w0=0.01)",
    params = list(beta_val=0.1, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.01, disc_val=0.001)
  ),

  # --- Condition 4: Only change β to Donahoe 1993 value ---
  list(
    label = "4. DDM-UI + β=0.035 only",
    params = list(beta_val=0.035, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 5: Only change threshold to Donahoe 1993 ---
  list(
    label = "5. DDM-UI + θ~Gauss(0,1) only",
    params = list(beta_val=0.12, decay_val=0.1, mu_thresh=0.0, sigma_thresh=1.0,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 6: Only change decay to Donahoe 1993 ---
  list(
    label = "6. DDM-UI + κ=0.05 only",
    params = list(beta_val=0.12, decay_val=0.05, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 7: Only change initial weight to Donahoe 1993 ---
  list(
    label = "7. DDM-UI + w0=0.01 only",
    params = list(beta_val=0.12, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.01, disc_val=0.001)
  ),

  # --- Condition 8: DDM-UI + disc=0.005 ---
  list(
    label = "8. DDM-UI + disc=0.005",
    params = list(beta_val=0.12, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.005)
  ),

  # --- Condition 9: β=0.035 + θ~Gauss(0,1) (Donahoe combo) ---
  list(
    label = "9. β=0.035 + θ~Gauss(0,1)",
    params = list(beta_val=0.035, decay_val=0.1, mu_thresh=0.0, sigma_thresh=1.0,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 10: Donahoe 1993 full but with w0=0.1 (DDM-UI init weight) ---
  list(
    label = "10. Donahoe 1993 + w0=0.1",
    params = list(beta_val=0.035, decay_val=0.05, mu_thresh=0.0, sigma_thresh=1.0,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 11: Burgos 2007 + w0=0.1 ---
  list(
    label = "11. Burgos 2007 + w0=0.1",
    params = list(beta_val=0.1, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.001)
  ),

  # --- Condition 12: DDM-UI + disc=0.05 (Burgos 2015 value) ---
  list(
    label = "12. DDM-UI + disc=0.05",
    params = list(beta_val=0.12, decay_val=0.1, mu_thresh=0.2, sigma_thresh=0.15,
                  threshold_type="gaussian", init_weight=0.1, disc_val=0.05)
  )
)

# =============================================================================
# RUN ALL CONDITIONS
# =============================================================================

cat("================================================================\n")
cat("SYSTEMATIC EXTINCTION PARAMETER COMPARISON\n")
cat(sprintf("Networks per condition: %d\n", NUM_NETWORKS))
cat(sprintf("Training trials: %d, Extinction trials: %d\n", TRAINING_TRIALS, EXTINCTION_TRIALS))
cat("================================================================\n")

all_results <- list()

for (cond in conditions) {
  cat(sprintf("\nRunning: %s ...\n", cond$label))
  start_time <- Sys.time()

  exp <- do.call(build_extinction_experiment, cond$params)
  raw_results <- run_condition(exp)
  analysis <- analyze_results(raw_results, cond$label)
  all_results[[cond$label]] <- analysis

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf("  (completed in %.1f seconds)\n", as.numeric(elapsed)))
}

# =============================================================================
# SUMMARY TABLE
# =============================================================================
cat("\n\n================================================================\n")
cat("SUMMARY TABLE\n")
cat("================================================================\n")
cat(sprintf("%-55s %8s %10s %10s\n", "Condition", "Success%", "M.1 mean", "M.1 sd"))
cat(paste(rep("-", 87), collapse=""), "\n")

for (res in all_results) {
  if (!is.null(res)) {
    cat(sprintf("%-55s %7.0f%% %10.3f %10.3f\n",
                res$label, res$success_rate,
                mean(res$m1_extinction), sd(res$m1_extinction)))
  }
}
cat("\n")
