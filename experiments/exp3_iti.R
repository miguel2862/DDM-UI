#!/usr/bin/env Rscript
# Experiment 3: ITI effect on extinction
# Compare with ITI vs without ITI

cat("=== EXPERIMENT 3: ITI EFFECT ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

template <- get_template_data("extinction")

set.seed(42)

n_networks <- 5
md_col <- "M..1-D"

# Condition 1: No ITI (baseline)
cat("--- Condition 1: No ITI (baseline) ---\n")
no_iti_train <- numeric(n_networks)
no_iti_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template$contingencies, template$trials)

  result <- Simulate.DBP(
    NPEs = template$npes,
    Connections = template$connections,
    TimeSteps = timesteps,
    HasITI = template$hasITI,
    threshold = "gaussian",
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  no_iti_train[net] <- as.numeric(result[max(training_rows), md_col])
  no_iti_ext[net] <- as.numeric(result[max(ext_rows), md_col])

  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, no_iti_train[net], no_iti_ext[net],
              (no_iti_ext[net] - no_iti_train[net]) / no_iti_train[net] * 100))
}

# Condition 2: With ITI (minITI=5, maxITI=10)
cat("\n--- Condition 2: With ITI (5-10 timesteps) ---\n")

# Need to add an ITI trial type and modify contingencies
# ITI trial: just the stimulus at 0 activation
template_iti <- template
template_iti$trials$ITI <- c("US,0.00,S1,0.00,True")

# Modify contingencies to include ITI for extinction phase
template_iti$contingencies <- c(
  "training, Random, Training, 100, False",
  "extinction, Random, Extinction, 100, True, 5, 10, ITI"
)
template_iti$hasITI <- c(FALSE, TRUE)

iti_train <- numeric(n_networks)
iti_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template_iti$contingencies, template_iti$trials)

  result <- Simulate.DBP(
    NPEs = template_iti$npes,
    Connections = template_iti$connections,
    TimeSteps = timesteps,
    HasITI = template_iti$hasITI,
    threshold = "gaussian",
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  iti_train[net] <- as.numeric(result[max(training_rows), md_col])
  iti_ext[net] <- as.numeric(result[max(ext_rows), md_col])

  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, iti_train[net], iti_ext[net],
              (iti_ext[net] - iti_train[net]) / iti_train[net] * 100))
}

# Condition 3: With ITI for BOTH phases
cat("\n--- Condition 3: With ITI for BOTH phases (5-10 timesteps) ---\n")

template_iti2 <- template
template_iti2$trials$ITI <- c("US,0.00,S1,0.00,True")

template_iti2$contingencies <- c(
  "training, Random, Training, 100, True, 5, 10, ITI",
  "extinction, Random, Extinction, 100, True, 5, 10, ITI"
)
template_iti2$hasITI <- c(TRUE, TRUE)

iti2_train <- numeric(n_networks)
iti2_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template_iti2$contingencies, template_iti2$trials)

  result <- Simulate.DBP(
    NPEs = template_iti2$npes,
    Connections = template_iti2$connections,
    TimeSteps = timesteps,
    HasITI = template_iti2$hasITI,
    threshold = "gaussian",
    disc = 0.001
  )

  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  iti2_train[net] <- as.numeric(result[max(training_rows), md_col])
  iti2_ext[net] <- as.numeric(result[max(ext_rows), md_col])

  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, iti2_train[net], iti2_ext[net],
              (iti2_ext[net] - iti2_train[net]) / iti2_train[net] * 100))
}

cat("\n\n--- Summary ---\n")
cat(sprintf("No ITI:        Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(no_iti_train), mean(no_iti_ext), mean((no_iti_ext - no_iti_train) / no_iti_train * 100)))
cat(sprintf("ITI (ext only): Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(iti_train), mean(iti_ext), mean((iti_ext - iti_train) / iti_train * 100)))
cat(sprintf("ITI (both):    Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(iti2_train), mean(iti2_ext), mean((iti2_ext - iti2_train) / iti2_train * 100)))
