# Experiment B: Reset OFF (HasITI=TRUE) vs Reset OFF + ITI vs Default
# Tests the effect of the activation reset mechanism

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

set.seed(42)

template <- get_template_data("extinction")
n_networks <- 5

cat("=== EXPERIMENT B: Reset Mechanism Effects ===\n\n")

results_B <- data.frame(
  Condition = character(),
  Network = integer(),
  M1D_end_training = numeric(),
  M1D_end_extinction = numeric(),
  M1D_change = numeric(),
  Pct_dVTA_above_disc = numeric(),
  stringsAsFactors = FALSE
)

# --- Condition 1: HasITI=TRUE for extinction phase, but no actual ITI trial ---
# Just disables the reset. No ITI timesteps added.
cat("--- Condition 1: HasITI=TRUE (no reset), no ITI timesteps ---\n")
for (net in 1:n_networks) {
  trials <- template$trials
  # Contingencies: training has HasITI=FALSE (default), extinction has HasITI=TRUE (no reset)
  # When HasITI=TRUE but no ITI params provided, we need the contingency to have True
  # but Create.Phases will need ITI params. Looking at the code:
  # if (has.iti) { min.ITI <- ...; max.ITI <- ...; ITITimestep <- ... }
  # So we need to provide them but set ITI length to 0... but that will crash sample(0:0)
  # Actually let's look more carefully - when has.iti=TRUE, it always tries to add ITI timesteps
  # before each trial. If min.ITI=max.ITI=0 => current.ITI=0 => for (ts in 1:0) runs 0 times...
  # Actually in R, 1:0 = c(1, 0) which would run twice. We need to handle this differently.
  #
  # Alternative: We'll create a custom HasITI vector and call Create.Phases manually for each phase.

  # Phase 1: training with HasITI=FALSE (normal)
  training_contingency <- "training, Random, Training, 100, False"
  # Phase 2: extinction with HasITI=TRUE, min=0, max=0, ITI trial needed but won't really run
  # We need an ITI trial type with 1 timestep to satisfy the format, and set ITI count = 1
  # Actually, let's use HasITI=TRUE with min=1, max=1 but make the ITI trial have zero stimuli
  # That way there's a minimal 1-timestep "ITI" that does nothing

  # Better approach: Directly set HasITI=TRUE but use the original contingencies.
  # The issue is Create.Phases parses the contingency string for ITI params.
  # Let's just create TimeSteps with HasITI=FALSE for both, then override HasITI when calling Simulate.

  TimeSteps <- Create.Phases(
    phases = template$contingencies,
    trials = trials
  )

  # Override HasITI: training=FALSE, extinction=TRUE (disables reset)
  hasITI_cond1 <- c(FALSE, TRUE)

  result <- Simulate.DBP(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = TimeSteps,
    HasITI = hasITI_cond1,
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  extinction_rows <- which(result$Phase == "extinction")

  m1d_col <- "M..1-D"
  m1d_end_training <- result[max(training_rows), m1d_col]
  m1d_end_extinction <- result[max(extinction_rows), m1d_col]
  m1d_change <- m1d_end_extinction - m1d_end_training

  ext_dvta <- result[extinction_rows, "dVTA"]
  pct_above <- sum(ext_dvta > 0.001, na.rm = TRUE) / length(ext_dvta) * 100

  results_B <- rbind(results_B, data.frame(
    Condition = "NoReset_NoITI",
    Network = net,
    M1D_end_training = round(m1d_end_training, 6),
    M1D_end_extinction = round(m1d_end_extinction, 6),
    M1D_change = round(m1d_change, 6),
    Pct_dVTA_above_disc = round(pct_above, 2),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Net %d: M''1->D training=%.4f, extinction=%.4f, change=%.4f, %%dVTA>disc=%.1f%%\n",
              net, m1d_end_training, m1d_end_extinction, m1d_change, pct_above))
}

# --- Condition 2: HasITI=TRUE + actual ITI of 30 timesteps ---
cat("\n--- Condition 2: HasITI=TRUE (no reset) + 30-timestep ITI ---\n")
for (net in 1:n_networks) {
  trials <- template$trials
  # Add an ITI trial type: 1 timestep with all zeros, learning active
  trials$ITI <- c("US,0.00,S1,0.00,True")

  # Modify extinction contingency to include ITI
  # Format: "phase_name, order, trial_types, trial_counts, hasITI, minITI, maxITI, ITItrialtype"
  contingencies_cond2 <- c(
    "training, Random, Training, 100, False",
    "extinction, Random, Extinction, 100, True, 30, 30, ITI"
  )

  TimeSteps <- Create.Phases(
    phases = contingencies_cond2,
    trials = trials
  )

  hasITI_cond2 <- c(FALSE, TRUE)

  result <- Simulate.DBP(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = TimeSteps,
    HasITI = hasITI_cond2,
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  extinction_rows <- which(result$Phase == "extinction")

  m1d_col <- "M..1-D"
  m1d_end_training <- result[max(training_rows), m1d_col]
  m1d_end_extinction <- result[max(extinction_rows), m1d_col]
  m1d_change <- m1d_end_extinction - m1d_end_training

  ext_dvta <- result[extinction_rows, "dVTA"]
  pct_above <- sum(ext_dvta > 0.001, na.rm = TRUE) / length(ext_dvta) * 100

  results_B <- rbind(results_B, data.frame(
    Condition = "NoReset_ITI30",
    Network = net,
    M1D_end_training = round(m1d_end_training, 6),
    M1D_end_extinction = round(m1d_end_extinction, 6),
    M1D_change = round(m1d_change, 6),
    Pct_dVTA_above_disc = round(pct_above, 2),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Net %d: M''1->D training=%.4f, extinction=%.4f, change=%.4f, %%dVTA>disc=%.1f%%\n",
              net, m1d_end_training, m1d_end_extinction, m1d_change, pct_above))
}

# --- Condition 3: Default (HasITI=FALSE, reset active) ---
cat("\n--- Condition 3: Default (HasITI=FALSE, reset active) ---\n")
for (net in 1:n_networks) {
  trials <- template$trials

  TimeSteps <- Create.Phases(
    phases = template$contingencies,
    trials = trials
  )

  result <- Simulate.DBP(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = TimeSteps,
    HasITI = template$hasITI,
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  extinction_rows <- which(result$Phase == "extinction")

  m1d_col <- "M..1-D"
  m1d_end_training <- result[max(training_rows), m1d_col]
  m1d_end_extinction <- result[max(extinction_rows), m1d_col]
  m1d_change <- m1d_end_extinction - m1d_end_training

  ext_dvta <- result[extinction_rows, "dVTA"]
  pct_above <- sum(ext_dvta > 0.001, na.rm = TRUE) / length(ext_dvta) * 100

  results_B <- rbind(results_B, data.frame(
    Condition = "Default_Reset",
    Network = net,
    M1D_end_training = round(m1d_end_training, 6),
    M1D_end_extinction = round(m1d_end_extinction, 6),
    M1D_change = round(m1d_change, 6),
    Pct_dVTA_above_disc = round(pct_above, 2),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("  Net %d: M''1->D training=%.4f, extinction=%.4f, change=%.4f, %%dVTA>disc=%.1f%%\n",
              net, m1d_end_training, m1d_end_extinction, m1d_change, pct_above))
}

cat("\n=== EXPERIMENT B SUMMARY ===\n")
for (cond in c("NoReset_NoITI", "NoReset_ITI30", "Default_Reset")) {
  subset <- results_B[results_B$Condition == cond, ]
  label <- switch(cond,
    "NoReset_NoITI" = "HasITI=TRUE, no ITI timesteps (no reset)",
    "NoReset_ITI30" = "HasITI=TRUE + 30-ts ITI (no reset)",
    "Default_Reset" = "Default (HasITI=FALSE, reset active)"
  )
  cat(sprintf("\n%s:\n", label))
  cat(sprintf("  Mean M''1->D end training:    %.4f (SD=%.4f)\n", mean(subset$M1D_end_training), sd(subset$M1D_end_training)))
  cat(sprintf("  Mean M''1->D end extinction:  %.4f (SD=%.4f)\n", mean(subset$M1D_end_extinction), sd(subset$M1D_end_extinction)))
  cat(sprintf("  Mean change:                  %.4f (SD=%.4f)\n", mean(subset$M1D_change), sd(subset$M1D_change)))
  cat(sprintf("  Mean %%dVTA > disc:            %.1f%% (SD=%.1f%%)\n", mean(subset$Pct_dVTA_above_disc), sd(subset$Pct_dVTA_above_disc)))
  cat(sprintf("  Pct change from training:     %.1f%%\n", mean(subset$M1D_change) / mean(subset$M1D_end_training) * 100))
}
