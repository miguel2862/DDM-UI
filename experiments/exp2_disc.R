#!/usr/bin/env Rscript
# Experiment 2: Vary disc criterion
# Test disc = 0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5

cat("=== EXPERIMENT 2: VARY DISC CRITERION ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

template <- get_template_data("extinction")

set.seed(42)

disc_values <- c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
n_networks <- 5
md_col <- "M..1-D"

cat(sprintf("%-10s | %-12s | %-12s | %-12s | %-8s | %-20s\n",
            "disc", "Train_end", "Ext_end", "Change", "Pct", "Extinct(>20%/any)"))
cat(paste(rep("-", 85), collapse=""), "\n")

for (disc_val in disc_values) {
  train_ends <- numeric(n_networks)
  ext_ends <- numeric(n_networks)

  for (net in 1:n_networks) {
    timesteps <- Create.Phases(template$contingencies, template$trials)

    result <- Simulate.DBP(
      NPEs = template$npes,
      Connections = template$connections,
      TimeSteps = timesteps,
      HasITI = template$hasITI,
      threshold = "gaussian",
      disc = disc_val
    )

    training_rows <- which(result$Phase == "training")
    ext_rows <- which(result$Phase == "extinction")
    train_ends[net] <- as.numeric(result[max(training_rows), md_col])
    ext_ends[net] <- as.numeric(result[max(ext_rows), md_col])
  }

  pct_changes <- (ext_ends - train_ends) / train_ends * 100
  n_ext20 <- sum(pct_changes < -20)
  n_any <- sum(pct_changes < 0)

  cat(sprintf("%-10.3f | %-12.6f | %-12.6f | %-12.6f | %-8.1f | %d/%d / %d/%d\n",
              disc_val,
              mean(train_ends), mean(ext_ends),
              mean(ext_ends - train_ends),
              mean(pct_changes),
              n_ext20, n_networks, n_any, n_networks))
}

# Now do a detailed trajectory for disc=0.1 to show the best case
cat("\n\n--- Detailed trajectory for disc=0.1, Network 1 ---\n")
timesteps <- Create.Phases(template$contingencies, template$trials)
set.seed(100)
result <- Simulate.DBP(
  NPEs = template$npes,
  Connections = template$connections,
  TimeSteps = timesteps,
  HasITI = template$hasITI,
  threshold = "gaussian",
  disc = 0.1
)

ext_data <- result[result$Phase == "extinction", ]
trials <- unique(ext_data$Trial)
cat("Trial, Weight_end, dVTA_mean, N_above_disc, N_below_disc\n")
for (tr in trials[seq(1, length(trials), by=5)]) {
  trial_data <- ext_data[ext_data$Trial == tr, ]
  dvta_vals <- as.numeric(trial_data$dVTA)
  last_w <- as.numeric(trial_data[nrow(trial_data), md_col])
  n_above <- sum(dvta_vals >= 0.1, na.rm=TRUE)
  n_below <- sum(dvta_vals < 0.1, na.rm=TRUE)
  cat(sprintf("  Trial %3d: w=%.6f, dVTA_mean=%.4f, above_disc=%d, below_disc=%d\n",
              as.numeric(tr), last_w, mean(dvta_vals, na.rm=TRUE), n_above, n_below))
}
