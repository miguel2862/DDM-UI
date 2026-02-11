#!/usr/bin/env Rscript
# Experiment 4: Test linear decrement formula
# Original 1993 formula: Δw = -β * a_pre * a_post (linear decay, no w multiplier)
# Current formula:        Δw = -β * w * a_pre * a_post (exponential decay)
#
# We modify Simulate.DBP inline to test linear decrement

cat("=== EXPERIMENT 4: LINEAR vs EXPONENTIAL DECREMENT ===\n\n")

source('/Users/miguel/Documents/DDM-UI/api/simulation.R')
source('/Users/miguel/Documents/DDM-UI/api/helpers.R')
source('/Users/miguel/Documents/DDM-UI/api/templates.R')

# Save original function
Simulate.DBP.original <- Simulate.DBP

# Create modified version with linear decrement (remove w from decay)
# We need to copy the function body and modify the decrement line
# The key change is on the line:
#   network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
#     beta * weight * a_pre * a_post
# becomes:
#   network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
#     beta * a_pre * a_post

# Get the function body as text
sim_code <- readLines('/Users/miguel/Documents/DDM-UI/api/simulation.R')

# Find and replace the decrement line
# Current code (lines ~362-369):
#   network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -
#     ifelse(network[[pre]]@Type == "Excitatory",
#       network[[i]]@InputConnections[[j]]@beta,
#       network[[i]]@InputConnections[[j]]@beta_prime
#     ) *
#       network[[i]]@InputConnections[[j]]@weight *
#       network[[pre]]@Activation *
#       network[[i]]@Activation

# We replace by removing the "* network[[i]]@InputConnections[[j]]@weight" line
modified_code <- gsub(
  'network\\[\\[i\\]\\]@InputConnections\\[\\[j\\]\\]@weight <- network\\[\\[i\\]\\]@InputConnections\\[\\[j\\]\\]@weight -\n                    ifelse\\(network\\[\\[pre\\]\\]@Type == "Excitatory",\n                      network\\[\\[i\\]\\]@InputConnections\\[\\[j\\]\\]@beta,\n                      network\\[\\[i\\]\\]@InputConnections\\[\\[j\\]\\]@beta_prime\n                    \\) \\*\n                      network\\[\\[i\\]\\]@InputConnections\\[\\[j\\]\\]@weight \\*\n                      network\\[\\[pre\\]\\]@Activation \\*\n                      network\\[\\[i\\]\\]@Activation',
  'network[[i]]@InputConnections[[j]]@weight <- network[[i]]@InputConnections[[j]]@weight -\n                    ifelse(network[[pre]]@Type == "Excitatory",\n                      network[[i]]@InputConnections[[j]]@beta,\n                      network[[i]]@InputConnections[[j]]@beta_prime\n                    ) *\n                      network[[pre]]@Activation *\n                      network[[i]]@Activation',
  paste(sim_code, collapse="\n")
)

# Write modified version to temp file
writeLines(modified_code, '/tmp/simulation_linear.R')

# Source the modified version (overwrites Simulate.DBP)
source('/tmp/simulation_linear.R')
Simulate.DBP.linear <- Simulate.DBP

# Restore original
source('/Users/miguel/Documents/DDM-UI/api/simulation.R')

template <- get_template_data("extinction")
n_networks <- 5
md_col <- "M..1-D"

# Run with ORIGINAL (exponential) decrement
cat("--- Exponential decrement (current): Δw = -β * w * a_pre * a_post ---\n")
set.seed(42)
exp_train <- numeric(n_networks)
exp_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template$contingencies, template$trials)
  result <- Simulate.DBP(
    NPEs = template$npes, Connections = template$connections,
    TimeSteps = timesteps, HasITI = template$hasITI,
    threshold = "gaussian", disc = 0.001
  )
  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  exp_train[net] <- as.numeric(result[max(training_rows), md_col])
  exp_ext[net] <- as.numeric(result[max(ext_rows), md_col])
  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, exp_train[net], exp_ext[net],
              (exp_ext[net] - exp_train[net]) / exp_train[net] * 100))
}

# Run with LINEAR decrement
cat("\n--- Linear decrement (modified): Δw = -β * a_pre * a_post ---\n")
set.seed(42)
lin_train <- numeric(n_networks)
lin_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template$contingencies, template$trials)
  result <- Simulate.DBP.linear(
    NPEs = template$npes, Connections = template$connections,
    TimeSteps = timesteps, HasITI = template$hasITI,
    threshold = "gaussian", disc = 0.001
  )
  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  lin_train[net] <- as.numeric(result[max(training_rows), md_col])
  lin_ext[net] <- as.numeric(result[max(ext_rows), md_col])
  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, lin_train[net], lin_ext[net],
              (lin_ext[net] - lin_train[net]) / lin_train[net] * 100))
}

# Run LINEAR with higher disc to also test
cat("\n--- Linear decrement with disc=0.005 ---\n")
set.seed(42)
lin5_train <- numeric(n_networks)
lin5_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template$contingencies, template$trials)
  result <- Simulate.DBP.linear(
    NPEs = template$npes, Connections = template$connections,
    TimeSteps = timesteps, HasITI = template$hasITI,
    threshold = "gaussian", disc = 0.005
  )
  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  lin5_train[net] <- as.numeric(result[max(training_rows), md_col])
  lin5_ext[net] <- as.numeric(result[max(ext_rows), md_col])
  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, lin5_train[net], lin5_ext[net],
              (lin5_ext[net] - lin5_train[net]) / lin5_train[net] * 100))
}

# Run LINEAR with disc=0.01
cat("\n--- Linear decrement with disc=0.01 ---\n")
set.seed(42)
lin01_train <- numeric(n_networks)
lin01_ext <- numeric(n_networks)

for (net in 1:n_networks) {
  timesteps <- Create.Phases(template$contingencies, template$trials)
  result <- Simulate.DBP.linear(
    NPEs = template$npes, Connections = template$connections,
    TimeSteps = timesteps, HasITI = template$hasITI,
    threshold = "gaussian", disc = 0.01
  )
  training_rows <- which(result$Phase == "training")
  ext_rows <- which(result$Phase == "extinction")
  lin01_train[net] <- as.numeric(result[max(training_rows), md_col])
  lin01_ext[net] <- as.numeric(result[max(ext_rows), md_col])
  cat(sprintf("  Network %d: train=%.6f, ext=%.6f, change=%.1f%%\n",
              net, lin01_train[net], lin01_ext[net],
              (lin01_ext[net] - lin01_train[net]) / lin01_train[net] * 100))
}

cat("\n\n--- Summary ---\n")
cat(sprintf("Exponential (disc=0.001): Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(exp_train), mean(exp_ext), mean((exp_ext - exp_train) / exp_train * 100)))
cat(sprintf("Linear (disc=0.001):      Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(lin_train), mean(lin_ext), mean((lin_ext - lin_train) / lin_train * 100)))
cat(sprintf("Linear (disc=0.005):      Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(lin5_train), mean(lin5_ext), mean((lin5_ext - lin5_train) / lin5_train * 100)))
cat(sprintf("Linear (disc=0.01):       Train=%.4f, Ext=%.4f, Change=%.1f%%\n",
            mean(lin01_train), mean(lin01_ext), mean((lin01_ext - lin01_train) / lin01_train * 100)))

# Also track per-trial trajectory for linear decrement
cat("\n\n--- Per-trial trajectory: Linear decrement disc=0.001, Network 1 ---\n")
set.seed(42)
timesteps <- Create.Phases(template$contingencies, template$trials)
result <- Simulate.DBP.linear(
  NPEs = template$npes, Connections = template$connections,
  TimeSteps = timesteps, HasITI = template$hasITI,
  threshold = "gaussian", disc = 0.001
)

# Show every 10 trials in both phases
for (phase_name in c("training", "extinction")) {
  cat(sprintf("\n  Phase: %s\n", phase_name))
  phase_data <- result[result$Phase == phase_name, ]
  trials <- unique(phase_data$Trial)
  for (tr in trials[seq(1, length(trials), by=10)]) {
    trial_rows <- phase_data[phase_data$Trial == tr, ]
    last_w <- as.numeric(trial_rows[nrow(trial_rows), md_col])
    cat(sprintf("    Trial %3d: w=%.6f\n", as.numeric(tr), last_w))
  }
}
