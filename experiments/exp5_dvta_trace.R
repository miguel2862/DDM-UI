#!/usr/bin/env Rscript
# Experiment 5: Track dVTA at every timestep during first 20 extinction trials
# Also track D activation and M''1 activation to understand the dynamics

cat("=== EXPERIMENT 5: dVTA PER-TIMESTEP TRACE ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

template <- get_template_data("extinction")

set.seed(42)

# Run 1 network
timesteps <- Create.Phases(template$contingencies, template$trials)

result <- Simulate.DBP(
  NPEs = template$npes,
  Connections = template$connections,
  TimeSteps = timesteps,
  HasITI = template$hasITI,
  threshold = "gaussian",
  disc = 0.001
)

md_col <- "M..1-D"

# --- Part A: Training phase, last 5 trials ---
cat("--- Part A: Last 5 training trials (dVTA trace) ---\n")
training_data <- result[result$Phase == "training", ]
train_trials <- unique(training_data$Trial)
last_5_train <- tail(train_trials, 5)

for (tr in last_5_train) {
  cat(sprintf("\nTraining Trial %s:\n", tr))
  trial_data <- training_data[training_data$Trial == tr, ]
  cat(sprintf("  %-3s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-8s\n",
              "TS", "US", "S1", "D", "S..1", "M..1", "M..1-D", "dVTA", "d>=disc"))
  for (row in 1:nrow(trial_data)) {
    r <- trial_data[row, ]
    dvta <- as.numeric(r$dVTA)
    cat(sprintf("  %-3s | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %+10.6f | %-8s\n",
                r$TimeStep,
                as.numeric(r$US),
                as.numeric(r$S1),
                as.numeric(r$D),
                as.numeric(r$S..1),
                as.numeric(r$M..1),
                as.numeric(r[[md_col]]),
                dvta,
                ifelse(dvta >= 0.001, "YES", "no")))
  }
}

# --- Part B: First 20 extinction trials ---
cat("\n\n--- Part B: First 20 extinction trials (dVTA trace) ---\n")
ext_data <- result[result$Phase == "extinction", ]
ext_trials <- unique(ext_data$Trial)

for (tr in ext_trials[1:min(20, length(ext_trials))]) {
  cat(sprintf("\nExtinction Trial %s:\n", tr))
  trial_data <- ext_data[ext_data$Trial == tr, ]
  cat(sprintf("  %-3s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-10s | %-8s\n",
              "TS", "US", "S1", "D", "S..1", "M..1", "M..1-D", "dVTA", "d>=disc"))
  for (row in 1:nrow(trial_data)) {
    r <- trial_data[row, ]
    dvta <- as.numeric(r$dVTA)
    cat(sprintf("  %-3s | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %-10.6f | %+10.6f | %-8s\n",
                r$TimeStep,
                as.numeric(r$US),
                as.numeric(r$S1),
                as.numeric(r$D),
                as.numeric(r$S..1),
                as.numeric(r$M..1),
                as.numeric(r[[md_col]]),
                dvta,
                ifelse(dvta >= 0.001, "YES", "no")))
  }
}

# --- Part C: Summary statistics for dVTA in extinction ---
cat("\n\n--- Part C: dVTA summary statistics across ALL 100 extinction trials ---\n")
cat(sprintf("%-10s | %-10s | %-10s | %-10s | %-10s | %-15s\n",
            "Timestep", "Mean_dVTA", "SD_dVTA", "Min_dVTA", "Max_dVTA", "Pct_above_disc"))

for (ts in 1:5) {
  ts_data <- ext_data[ext_data$TimeStep == ts, ]
  dvta_vals <- as.numeric(ts_data$dVTA)
  pct_above <- sum(dvta_vals >= 0.001, na.rm=TRUE) / length(dvta_vals) * 100
  cat(sprintf("%-10d | %+10.6f | %-10.6f | %+10.6f | %+10.6f | %-15.1f\n",
              ts,
              mean(dvta_vals, na.rm=TRUE),
              sd(dvta_vals, na.rm=TRUE),
              min(dvta_vals, na.rm=TRUE),
              max(dvta_vals, na.rm=TRUE),
              pct_above))
}

# --- Part D: What drives dVTA positive during extinction? ---
cat("\n\n--- Part D: D activation dynamics during extinction ---\n")
cat("D(t) - D(t-1) = dVTA. If M''1->D connection is strong and M''1 fires, D activates.\n")
cat("In extinction with no US, D should be driven by M''1->D connection.\n\n")

# Check: when dVTA > disc in extinction, what are D and M''1 doing?
n_above <- 0
n_below <- 0
d_when_above <- c()
m_when_above <- c()
w_when_above <- c()
d_when_below <- c()
m_when_below <- c()

for (row in 1:nrow(ext_data)) {
  dvta <- as.numeric(ext_data[row, "dVTA"])
  if (dvta >= 0.001) {
    n_above <- n_above + 1
    d_when_above <- c(d_when_above, as.numeric(ext_data[row, "D"]))
    m_when_above <- c(m_when_above, as.numeric(ext_data[row, "M..1"]))
    w_when_above <- c(w_when_above, as.numeric(ext_data[row, md_col]))
  } else {
    n_below <- n_below + 1
    d_when_below <- c(d_when_below, as.numeric(ext_data[row, "D"]))
    m_when_below <- c(m_when_below, as.numeric(ext_data[row, "M..1"]))
  }
}

cat(sprintf("Timesteps where dVTA >= disc: %d/%d (%.1f%%)\n", n_above, nrow(ext_data), n_above/nrow(ext_data)*100))
cat(sprintf("Timesteps where dVTA < disc: %d/%d (%.1f%%)\n\n", n_below, nrow(ext_data), n_below/nrow(ext_data)*100))

cat(sprintf("When dVTA >= disc:\n"))
cat(sprintf("  D activation: Mean=%.6f, SD=%.6f\n", mean(d_when_above), sd(d_when_above)))
cat(sprintf("  M''1 activation: Mean=%.6f, SD=%.6f\n", mean(m_when_above), sd(m_when_above)))
cat(sprintf("  M''1->D weight: Mean=%.6f, SD=%.6f\n", mean(w_when_above), sd(w_when_above)))

cat(sprintf("\nWhen dVTA < disc:\n"))
cat(sprintf("  D activation: Mean=%.6f, SD=%.6f\n", mean(d_when_below), sd(d_when_below)))
cat(sprintf("  M''1 activation: Mean=%.6f, SD=%.6f\n", mean(m_when_below), sd(m_when_below)))

# --- Part E: Net increment vs decrement per trial ---
cat("\n\n--- Part E: Per-trial increment vs decrement balance ---\n")
cat("For each trial, count timesteps with increment (dVTA>=disc) vs decrement (dVTA<disc)\n\n")
cat(sprintf("%-8s | %-8s | %-8s | %-12s | %-12s\n",
            "Trial", "Incr_TS", "Decr_TS", "Weight_start", "Weight_end"))

for (tr in ext_trials[1:20]) {
  trial_data <- ext_data[ext_data$Trial == tr, ]
  dvta_vals <- as.numeric(trial_data$dVTA)
  n_incr <- sum(dvta_vals >= 0.001, na.rm=TRUE)
  n_decr <- sum(dvta_vals < 0.001, na.rm=TRUE)
  w_start <- as.numeric(trial_data[1, md_col])
  w_end <- as.numeric(trial_data[nrow(trial_data), md_col])
  cat(sprintf("%-8s | %-8d | %-8d | %-12.6f | %-12.6f\n",
              tr, n_incr, n_decr, w_start, w_end))
}
