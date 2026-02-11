#!/usr/bin/env Rscript
# Experiment 1: Baseline extinction with current template
# Track M''1->D connection weight across trials in extinction phase

cat("=== EXPERIMENT 1: BASELINE EXTINCTION ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

template <- get_template_data("extinction")

set.seed(42)

n_networks <- 10
results <- list()

for (net in 1:n_networks) {
  cat(sprintf("Running network %d/%d...\n", net, n_networks))

  timesteps <- Create.Phases(template$contingencies, template$trials)

  result <- Simulate.DBP(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = timesteps,
    HasITI = template$hasITI,
    threshold = "gaussian",
    disc = 0.001
  )

  results[[net]] <- result
}

cat("\n--- Results ---\n")
cat("Connection column name for M''1->D: 'M..1-D'\n\n")

# Check available columns
cat("Available columns:\n")
cat(paste(colnames(results[[1]]), collapse=", "), "\n\n")

# For each network, get end-of-trial weights for M''1->D
# Get last timestep per trial in each phase
for (net in 1:n_networks) {
  r <- results[[net]]

  # Get training phase end weight
  training_rows <- which(r$Phase == "training")
  if (length(training_rows) > 0) {
    training_end <- r[max(training_rows), ]
  }

  # Get extinction phase end weight
  ext_rows <- which(r$Phase == "extinction")
  if (length(ext_rows) > 0) {
    ext_end <- r[max(ext_rows), ]
  }

  # M''1->D weight
  md_col <- "M..1-D"
  if (md_col %in% colnames(r)) {
    train_w <- as.numeric(training_end[[md_col]])
    ext_w <- as.numeric(ext_end[[md_col]])
    change <- ext_w - train_w
    pct_change <- ifelse(train_w > 0, (change / train_w) * 100, NA)
    cat(sprintf("Network %2d: Training end = %.6f, Extinction end = %.6f, Change = %+.6f (%.1f%%)\n",
                net, train_w, ext_w, change, pct_change))
  } else {
    cat(sprintf("Network %2d: Column '%s' not found\n", net, md_col))
    cat("  Trying to find connection columns...\n")
    conn_cols <- grep("-", colnames(r), value=TRUE)
    cat("  Connection columns:", paste(conn_cols, collapse=", "), "\n")
  }
}

# Track M''1->D across ALL extinction trials (per-trial trajectory)
cat("\n\n--- Per-trial M''1->D weight trajectory for Network 1 ---\n")
r <- results[[1]]
md_col <- "M..1-D"

if (md_col %in% colnames(r)) {
  # Get last timestep of each trial in extinction phase
  ext_data <- r[r$Phase == "extinction", ]
  trials <- unique(ext_data$Trial)

  cat("Trial, Weight_at_end_of_trial\n")
  for (tr in trials[seq(1, length(trials), by=10)]) {
    trial_rows <- which(ext_data$Trial == tr)
    last_ts <- ext_data[max(trial_rows), ]
    cat(sprintf("  Trial %3d: %.6f\n", as.numeric(tr), as.numeric(last_ts[[md_col]])))
  }
}

# Summary statistics
cat("\n\n--- Summary ---\n")
if (md_col %in% colnames(results[[1]])) {
  train_ends <- numeric(n_networks)
  ext_ends <- numeric(n_networks)

  for (net in 1:n_networks) {
    r <- results[[net]]
    training_rows <- which(r$Phase == "training")
    ext_rows <- which(r$Phase == "extinction")
    train_ends[net] <- as.numeric(r[max(training_rows), md_col])
    ext_ends[net] <- as.numeric(r[max(ext_rows), md_col])
  }

  cat(sprintf("Training end M''1->D: Mean=%.6f, SD=%.6f, Range=[%.6f, %.6f]\n",
              mean(train_ends), sd(train_ends), min(train_ends), max(train_ends)))
  cat(sprintf("Extinction end M''1->D: Mean=%.6f, SD=%.6f, Range=[%.6f, %.6f]\n",
              mean(ext_ends), sd(ext_ends), min(ext_ends), max(ext_ends)))
  cat(sprintf("Mean change: %.6f (%.1f%%)\n",
              mean(ext_ends - train_ends), mean((ext_ends - train_ends) / train_ends * 100)))

  # Count networks showing significant extinction (>20% decrease)
  pct_changes <- (ext_ends - train_ends) / train_ends * 100
  n_extinct <- sum(pct_changes < -20)
  cat(sprintf("\nNetworks showing >20%% extinction: %d/%d\n", n_extinct, n_networks))
  n_extinct10 <- sum(pct_changes < -10)
  cat(sprintf("Networks showing >10%% extinction: %d/%d\n", n_extinct10, n_networks))
  n_extinct5 <- sum(pct_changes < -5)
  cat(sprintf("Networks showing >5%% extinction: %d/%d\n", n_extinct5, n_networks))
  n_any_decrease <- sum(pct_changes < 0)
  cat(sprintf("Networks showing ANY decrease: %d/%d\n", n_any_decrease, n_networks))
}

# Also track dVTA during extinction
cat("\n\n--- dVTA during extinction for Network 1 ---\n")
r <- results[[1]]
ext_data <- r[r$Phase == "extinction", ]
trials <- unique(ext_data$Trial)

cat("Trial, Mean_dVTA, Min_dVTA, Max_dVTA, N_timesteps_above_disc\n")
for (tr in trials[seq(1, length(trials), by=10)]) {
  trial_rows <- ext_data[ext_data$Trial == tr, ]
  dvta_vals <- as.numeric(trial_rows$dVTA)
  n_above <- sum(dvta_vals >= 0.001, na.rm=TRUE)
  cat(sprintf("  Trial %3d: Mean=%.6f, Min=%.6f, Max=%.6f, N_above_disc=%d/%d\n",
              as.numeric(tr), mean(dvta_vals, na.rm=TRUE), min(dvta_vals, na.rm=TRUE),
              max(dvta_vals, na.rm=TRUE), n_above, length(dvta_vals)))
}
