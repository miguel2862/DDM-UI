#!/usr/bin/env Rscript
# Experiment 6: Quantify increment vs decrement magnitudes during extinction
# The core question: How much weight is ADDED vs REMOVED per trial?

cat("=== EXPERIMENT 6: INCREMENT vs DECREMENT MAGNITUDE ANALYSIS ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

template <- get_template_data("extinction")
set.seed(42)

# Run a single detailed simulation
timesteps <- Create.Phases(template$contingencies, template$trials)
result <- Simulate.DBP(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = timesteps, HasITI = template$hasITI,
  threshold = "gaussian", disc = 0.001
)

md_col <- "M..1-D"

# For extinction phase, compute per-timestep weight changes
ext_data <- result[result$Phase == "extinction", ]
ext_trials <- unique(ext_data$Trial)

total_increments <- 0
total_decrements <- 0
n_incr_ts <- 0
n_decr_ts <- 0

cat("Per-trial weight change breakdown:\n")
cat(sprintf("%-8s | %-12s | %-12s | %-12s | %-12s | %-12s\n",
            "Trial", "W_start", "W_end", "Net_change", "Sum_incr", "Sum_decr"))

for (tr in ext_trials) {
  trial_data <- ext_data[ext_data$Trial == tr, ]
  weights <- as.numeric(trial_data[[md_col]])

  # Per-timestep changes
  incr_sum <- 0
  decr_sum <- 0
  for (i in 2:length(weights)) {
    dw <- weights[i] - weights[i-1]
    if (dw > 0) {
      incr_sum <- incr_sum + dw
      n_incr_ts <- n_incr_ts + 1
    } else if (dw < 0) {
      decr_sum <- decr_sum + dw
      n_decr_ts <- n_decr_ts + 1
    }
  }
  # Also count first timestep change from previous trial
  total_increments <- total_increments + incr_sum
  total_decrements <- total_decrements + decr_sum

  if (as.numeric(tr) %% 10 == 0 || as.numeric(tr) <= 5) {
    cat(sprintf("%-8s | %-12.6f | %-12.6f | %+12.6f | %+12.6f | %+12.6f\n",
                tr, weights[1], weights[length(weights)],
                weights[length(weights)] - weights[1],
                incr_sum, decr_sum))
  }
}

cat(sprintf("\n\nOverall during extinction:\n"))
cat(sprintf("Total weight increments: %+.6f across %d timesteps\n", total_increments, n_incr_ts))
cat(sprintf("Total weight decrements: %+.6f across %d timesteps\n", total_decrements, n_decr_ts))
cat(sprintf("Increment/Decrement ratio: %.2f\n", abs(total_increments/total_decrements)))

# ---- Analyze the INCREMENT mechanism ----
cat("\n\n=== WHY DOES dVTA EXCEED DISC DURING EXTINCTION? ===\n")
cat("The M''1->D connection is strong (~0.7-0.9 after training).\n")
cat("When S1 is presented, activation cascades: S1->S''1->M''1->D.\n")
cat("D activates because M''1 is active and M''1->D weight is high.\n")
cat("dVTA = D(t) - D(t-1). On early timesteps of a trial:\n")
cat("  - Timestep 1: D resets to ~0 (no HasITI), then M''1 drives D up.\n")
cat("  - dVTA = D(new) - D(reset_to_0) = D(new), which is LARGE.\n")
cat("  => The INCREMENT rule fires because the CS itself predicts D activity!\n\n")

# ---- The self-reinforcing cycle ----
cat("=== THE SELF-REINFORCING CYCLE ===\n")
cat("1. M''1->D is strong from training\n")
cat("2. During extinction, CS activates M''1 (via S1->S''1->M''1)\n")
cat("3. M''1 drives D via the strong M''1->D connection\n")
cat("4. D activation increase => positive dVTA\n")
cat("5. dVTA > disc => INCREMENT rule fires => M''1->D gets STRONGER\n")
cat("6. Stronger M''1->D => more D activation next trial => goto 2\n\n")
cat("The connection that SHOULD be decremented is instead being INCREMENTED\n")
cat("because the network's own prediction of reinforcement creates a positive dVTA signal.\n")

# ---- Test: What if we track all connections ----
cat("\n\n=== ALL CONNECTION WEIGHTS: END OF TRAINING vs END OF EXTINCTION ===\n")
conn_cols <- grep("-", colnames(result), value=TRUE)
training_rows <- which(result$Phase == "training")
ext_rows <- which(result$Phase == "extinction")

cat(sprintf("%-20s | %-12s | %-12s | %-12s\n", "Connection", "Train_end", "Ext_end", "Change"))
for (cc in conn_cols) {
  train_w <- as.numeric(result[max(training_rows), cc])
  ext_w <- as.numeric(result[max(ext_rows), cc])
  cat(sprintf("%-20s | %-12.6f | %-12.6f | %+12.6f\n", cc, train_w, ext_w, ext_w - train_w))
}

# ---- Test: What about the Burgos 2000 extinction template? ----
cat("\n\n=== COMPARISON: BURGOS 2000 EXTINCTION TEMPLATE ===\n")
cat("This template uses different parameters (beta=0.035, mu=0, sigma=1)\n\n")

template_b <- get_template_data("burgos_2000")
set.seed(42)

timesteps_b <- Create.Phases(template_b$contingencies, template_b$trials)
result_b <- Simulate.DBP(
  NPEs = template_b$npes, Connections = template_b$connections,
  TimeSteps = timesteps_b, HasITI = template_b$hasITI,
  threshold = "gaussian", disc = 0.001
)

# Find M''1->D equivalent
conn_cols_b <- grep("-D$", colnames(result_b), value=TRUE)
cat("Connections to D:", paste(conn_cols_b, collapse=", "), "\n\n")

phases_b <- unique(result_b$Phase)
cat("Phases:", paste(phases_b, collapse=", "), "\n\n")

for (cc in conn_cols_b) {
  cat(sprintf("Connection %s:\n", cc))
  for (ph in phases_b) {
    ph_rows <- which(result_b$Phase == ph)
    if (length(ph_rows) > 0) {
      start_w <- as.numeric(result_b[min(ph_rows), cc])
      end_w <- as.numeric(result_b[max(ph_rows), cc])
      cat(sprintf("  %s: start=%.6f, end=%.6f, change=%+.6f\n", ph, start_w, end_w, end_w - start_w))
    }
  }
  cat("\n")
}
